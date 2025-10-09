#ifndef FE_DATA_OBJECT_H
#define FE_DATA_OBJECT_H

#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <unordered_map>
#include <fstream>
#include <optional>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include <vtkTable.h>
#include <vtkAbstractArray.h>
#include <unordered_map>
#include <tuple>
#include <string>
#include <vector>
#include <memory>

using namespace libMesh;

class LoadStep;

typedef std::vector<std::vector<std::variant<int, Real, std::string>>> DataframeMap;
typedef std::unordered_map<std::string, std::variant<int, Real, std::string>> SettingsMap;

// FE data object classes are the basic data repositories of the program.
// They store information into a map (dict) that can be later retrieved when needed
// Each type of data e.g. Materials, Real Constants, Solver Settings, Time Steps, etc. is stored into a unordered map.
// The basic classes can be enriched by peculiar methods (e.g. the ElementLoadsObject has methods to interface with the
// ElementBuilder class that builds the RHS of the elements)
class FEDataObjectAbstract
{
public:
    // Virtual destructor for proper cleanup of derived classes
    virtual ~FEDataObjectAbstract() = 0;
};

inline FEDataObjectAbstract::~FEDataObjectAbstract() {}

class ApecMMesh : public Mesh, public FEDataObjectAbstract
{
public:
    // Default constructor
    using Mesh::Mesh;

    // Virtual destructor for proper cleanup of derived classes
    virtual ~ApecMMesh() = default;
};

template <typename T1>
class FEDataObjectMap : public FEDataObjectAbstract
{
public:
    // Virtual destructor for proper cleanup of derived classes
    virtual ~FEDataObjectMap() = default;

    // virtual method to read data from a binary file stream
    virtual void read_data(std::ifstream &file) {};

    // Getter for the data map
    T1 *get_data_map() const
    {
        return data_map.get();
    }

    T1 *get_data_map()
    {
        return data_map.get();
    }

protected:
    FEDataObjectMap() : data_map(std::make_shared<T1>())
    {
    }
    // Protected constructor to prevent direct instantiation
    FEDataObjectMap(std::ifstream &file) : data_map(std::make_shared<T1>())
    {
        this->read_data(file);
    }

    FEDataObjectMap(std::shared_ptr<T1> external_map)
        : data_map(external_map) {}

    // // NEW: Constructor that receives an external map and a file, sets the map and calls read_data
    // FEDataObjectMap(std::ifstream &file, std::shared_ptr<T1> external_map)
    //     : data_map(external_map)
    // {
    //     this->read_data(file);
    // }

    // Attributes
    std::shared_ptr<T1> data_map; // Stores the data map (type-specific)
};

class DataframeFEDataObject : public FEDataObjectMap<DataframeMap>
{
public:
    // Constructor
    DataframeFEDataObject()
        : FEDataObjectMap<DataframeMap>() // Call base default constructor (does NOT read data)
    {
    }
    DataframeFEDataObject(std::ifstream &file)
        : FEDataObjectMap<DataframeMap>() // Call base default constructor (does NOT read data)
    {
        this->read_data(file); // Now call the child’s read_data
    }
    // Implement the read_data method
    void read_data(std::ifstream &file);
};

template <typename T>
class DataframeMapFEDataObject : public DataframeFEDataObject
{
public:
    // Constructor
    DataframeMapFEDataObject(std::ifstream &file)
        : DataframeFEDataObject(file), // parent constructor
          dataframe_map(std::make_shared<std::unordered_map<int, std::vector<T>>>())
    {
        populate_dataframe_map();
    }

    DataframeMapFEDataObject()
        : DataframeFEDataObject(),
          dataframe_map(std::make_shared<std::unordered_map<int, std::vector<T>>>())
    {
    }

    // Method to populate dataframe_map from data_map
    void populate_dataframe_map();

    // Getter for dataframe_map
    std::unordered_map<int, std::vector<T>> *get_dataframe_map() const
    {
        return dataframe_map.get();
    }

protected:
    std::shared_ptr<std::unordered_map<int, std::vector<T>>> dataframe_map; // Shared pointer to the map
};

class SettingsFEDataObject : public FEDataObjectMap<SettingsMap>
{
public:
    // Constructor
    SettingsFEDataObject(std::ifstream &file)
        : FEDataObjectMap<SettingsMap>() // Call base default constructor (does NOT read data)
    {
        this->read_data(file); // Now call the child’s read_data
    };

    template <typename T>
    T get_setting(const std::string &key) const
    {
        auto it = data_map->find(key);
        if (it == data_map->end())
        {
            libmesh_error_msg("Setting not found: " + key);
        }
        return std::get<T>(it->second);
    }

    template <typename T>
    void set_setting(const std::string &key, T value)
    {
        (*data_map)[key] = value;
    }
    // Implement the read_data method
    void read_data(std::ifstream &file);
};

template class DataframeMapFEDataObject<int>;
template class DataframeMapFEDataObject<Real>;

#include "FEDataObjectMap.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_logging.h"
#include <vtkTable.h>
#include <vtkAbstractArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>
#include <unordered_map>
#include <tuple>
#include <string>
#include <vector>
#include <memory>

template <typename TMap>
class vtkTableFEDataObject : public FEDataObjectMap<TMap>
{
public:
    vtkTableFEDataObject() = default;
    virtual ~vtkTableFEDataObject() = default;

    /**
     * Reads selected arrays from a vtkTable and populates the data_map.
     *
     * Each key in field_names can represent:
     *   - a scalar array (no "_shape" / "_offset")
     *   - an iterable array, accompanied by "<name>_shape" and "<name>_offset"
     *
     * The function supports:
     *   - int, double, string scalars
     *   - iterables (1D–ND) of ints, doubles, or strings
     *
     * Throws a libMesh error message if inconsistencies are found.
     */
    void read_from_vtk_table(vtkTable *table, const std::vector<std::string> &field_names)
    {
        if (!table)
            libmesh_error_msg("vtkTableFEDataObject::read_from_vtk_table: null vtkTable.");

        vtkIdType num_rows = 0;
        if (!field_names.empty())
        {
            auto *arr = table->GetColumnByName(field_names.front().c_str());
            if (arr)
                num_rows = arr->GetNumberOfTuples();
        }

        auto &data_map = *(this->data_map);
        data_map.clear();

        for (vtkIdType irow = 0; irow < num_rows; ++irow)
            data_map[irow] = std::tuple<>();

        for (const auto &base_name : field_names)
        {
            vtkAbstractArray *val_array = table->GetColumnByName(base_name.c_str());
            vtkIntArray *shape_array = vtkIntArray::SafeDownCast(table->GetColumnByName((base_name + "_shape").c_str()));
            vtkIntArray *offset_array = vtkIntArray::SafeDownCast(table->GetColumnByName((base_name + "_offset").c_str()));

            if (!val_array)
            {
                libmesh_warning("Field '" << base_name << "' not found — skipping.");
                continue;
            }

            // ────────────────────────────────
            // SCALAR FIELD
            // ────────────────────────────────
            if (!shape_array || !offset_array)
            {
                if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
                {
                    for (vtkIdType i = 0; i < num_rows; ++i)
                        append_to_tuple(data_map[i], int_arr->GetValue(i));
                }
                else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
                {
                    for (vtkIdType i = 0; i < num_rows; ++i)
                        append_to_tuple(data_map[i], dbl_arr->GetValue(i));
                }
                else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
                {
                    for (vtkIdType i = 0; i < num_rows; ++i)
                        append_to_tuple(data_map[i], std::string(str_arr->GetValue(i)));
                }
                else
                {
                    libmesh_error_msg("Unsupported scalar array type for field '" << base_name << "'.");
                }
                continue;
            }

            // ────────────────────────────────
            // ITERABLE FIELD
            // ────────────────────────────────
            const int *shape_data = shape_array->GetPointer(0);
            const int *offset_data = offset_array->GetPointer(0);

            // integers
            if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
            {
                const int *vals = int_arr->GetPointer(0);
                vtkIdType cursor = 0;

                for (vtkIdType irow = 0; irow < num_rows; ++irow)
                {
                    int shape_start = offset_data[irow];
                    int ndim = shape_data[shape_start];

                    std::vector<int> dims(ndim);
                    vtkIdType total_elems = 1;
                    for (int j = 0; j < ndim; ++j)
                    {
                        dims[j] = shape_data[shape_start + 1 + j];
                        total_elems *= dims[j];
                    }

                    std::vector<int> flattened(vals + cursor, vals + cursor + total_elems);
                    cursor += total_elems;

                    append_to_tuple(data_map[irow], reshape(flattened, dims));
                }
            }
            // doubles
            else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
            {
                const double *vals = dbl_arr->GetPointer(0);
                vtkIdType cursor = 0;

                for (vtkIdType irow = 0; irow < num_rows; ++irow)
                {
                    int shape_start = offset_data[irow];
                    int ndim = shape_data[shape_start];

                    std::vector<int> dims(ndim);
                    vtkIdType total_elems = 1;
                    for (int j = 0; j < ndim; ++j)
                    {
                        dims[j] = shape_data[shape_start + 1 + j];
                        total_elems *= dims[j];
                    }

                    std::vector<double> flattened(vals + cursor, vals + cursor + total_elems);
                    cursor += total_elems;

                    append_to_tuple(data_map[irow], reshape(flattened, dims));
                }
            }
            // strings
            else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
            {
                vtkIdType cursor = 0;
                for (vtkIdType irow = 0; irow < num_rows; ++irow)
                {
                    int shape_start = offset_data[irow];
                    int ndim = shape_data[shape_start];

                    std::vector<int> dims(ndim);
                    vtkIdType total_elems = 1;
                    for (int j = 0; j < ndim; ++j)
                    {
                        dims[j] = shape_data[shape_start + 1 + j];
                        total_elems *= dims[j];
                    }

                    std::vector<std::string> flat;
                    flat.reserve(total_elems);
                    for (vtkIdType k = 0; k < total_elems; ++k)
                        flat.emplace_back(str_arr->GetValue(cursor++));

                    append_to_tuple(data_map[irow], reshape(flat, dims));
                }
            }
            else
            {
                libmesh_error_msg("Unsupported iterable array type for field '" << base_name << "'.");
            }
        }
    }

private:
    // append a new element to an existing tuple
    template <typename Tuple, typename T>
    static void append_to_tuple(Tuple &t, T &&value)
    {
        t = std::tuple_cat(t, std::make_tuple(std::forward<T>(value)));
    }

    // Recursive reshape for flattened vectors (supports ndim ≥ 1)
    template <typename T>
    static std::vector<std::vector<T>> reshape(const std::vector<T> &flat,
                                               const std::vector<int> &dims,
                                               size_t level = 0,
                                               size_t offset = 0)
    {
        if (dims.empty())
            return {};

        if (level == dims.size() - 1)
        {
            std::vector<std::vector<T>> leaf;
            size_t chunk = flat.size() / dims[level];
            for (int i = 0; i < dims[level]; ++i)
                leaf.push_back(std::vector<T>(flat.begin() + offset + i * chunk,
                                              flat.begin() + offset + (i + 1) * chunk));
            return leaf;
        }

        size_t stride = flat.size() / dims[level];
        std::vector<std::vector<T>> result;
        for (int i = 0; i < dims[level]; ++i)
        {
            auto sub = reshape(flat, dims, level + 1, offset + i * stride);
            result.insert(result.end(), sub.begin(), sub.end());
        }
        return result;
    }
};
// class InfiniteElementsOptionsObject : public DataframeMatLikeFEDataObject
// {
// public:
//     InfiniteElementsOptionsObject() = default;

//     // Populate from vtkFieldData (may have multiple tuples)
//     void populate_from_vtk_field_data(vtkFieldData *field_data)
//     {
//         // Clear previous data
//         data_map->clear();

//         if (!field_data)
//             return;

//         // Get arrays
//         vtkIntArray *id_array = vtkIntArray::SafeDownCast(field_data->GetArray("infinite_boundary_ID"));
//         vtkDoubleArray *origin_x_array = vtkDoubleArray::SafeDownCast(field_data->GetArray("origin_x"));
//         vtkDoubleArray *origin_y_array = vtkDoubleArray::SafeDownCast(field_data->GetArray("origin_y"));
//         vtkDoubleArray *origin_z_array = vtkDoubleArray::SafeDownCast(field_data->GetArray("origin_z"));

//         // Check arrays
//         if (!id_array || !origin_x_array || !origin_y_array || !origin_z_array)
//             return;

//         vtkIdType n = id_array->GetNumberOfTuples();
//         if (origin_x_array->GetNumberOfTuples() != n ||
//             origin_y_array->GetNumberOfTuples() != n ||
//             origin_z_array->GetNumberOfTuples() != n)
//             return;

//         // Fill data_map (each row: int, Real, Real, Real)
//         for (vtkIdType i = 0; i < n; ++i)
//         {
//             std::vector<std::variant<int, Real, std::string>> row;
//             row.push_back(static_cast<int>(id_array->GetValue(i)));
//             row.push_back(static_cast<Real>(origin_x_array->GetValue(i)));
//             row.push_back(static_cast<Real>(origin_y_array->GetValue(i)));
//             row.push_back(static_cast<Real>(origin_z_array->GetValue(i)));
//             data_map->push_back(std::move(row));
//         }

//         // Now populate mat_like_map from data_map
//         populate_dataframe_map();
//     }
// };

#endif // FE_DATA_OBJECT_H