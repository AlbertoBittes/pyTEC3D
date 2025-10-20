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

class FEMDataRepository;
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
    // virtual void broadcast(const libMesh::Parallel::Communicator &comm) {};
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

    // virtual void broadcast(const libMesh::Parallel::Communicator &comm) override
    // {
    // }

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

    // NEW: Constructor that copies from an external map
    FEDataObjectMap(const T1 &external_map)
        : data_map(std::make_shared<T1>(external_map)) {}
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

// template <typename TMap>
// class vtkFieldDataFEDataObject : public FEDataObjectMap<TMap>, public vtkFieldDataFEDataObjectBase
// {
// public:
//     using VariantType = std::variant<
//         int,
//         Real,
//         std::string,
//         std::vector<int>,
//         std::vector<Real>,
//         std::vector<std::string>>;

//     vtkFieldDataFEDataObject() = default;
//     virtual ~vtkFieldDataFEDataObject() = default;

//     void read_from_vtk_field_data(
//         vtkFieldData *field_data, const std::vector<std::string> &field_names)
//     {
//         if (!field_data)
//             libmesh_error_msg("vtkFieldDataFEDataObject::read_from_vtk_field_data: null vtkFieldData.");

//         vtkIdType num_rows = 0;
//         if (!field_names.empty())
//         {
//             auto *arr = field_data->GetAbstractArray(field_names.front().c_str());
//             if (arr)
//                 num_rows = arr->GetNumberOfTuples();
//         }

//         auto &data_map = *(this->data_map);
//         data_map.clear();
//         data_map.reserve(static_cast<size_t>(num_rows));

//         for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//             data_map.emplace(irow, typename TMap::mapped_type{});

//         for (const auto &base_name : field_names)
//         {
//             vtkAbstractArray *val_array = field_data->GetAbstractArray(base_name.c_str());
//             vtkIntArray *shape_array = vtkIntArray::SafeDownCast(
//                 field_data->GetAbstractArray((base_name + "_shape").c_str()));
//             vtkIntArray *offset_array = vtkIntArray::SafeDownCast(
//                 field_data->GetAbstractArray((base_name + "_offset").c_str()));

//             if (!val_array)
//             {
//                 libmesh_warning("Field '" << base_name << "' not found — skipping.");
//                 continue;
//             }

//             // --------------------------
//             // SCALAR FIELDS
//             // --------------------------
//             if (!shape_array || !offset_array)
//             {
//                 if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
//                 {
//                     for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//                         data_map[irow].emplace_back(int_arr->GetValue(irow));
//                 }
//                 else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
//                 {
//                     for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//                         data_map[irow].emplace_back(dbl_arr->GetValue(irow));
//                 }
//                 else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
//                 {
//                     for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//                         data_map[irow].emplace_back(std::string(str_arr->GetValue(irow)));
//                 }
//                 else
//                 {
//                     libmesh_error_msg("Unsupported scalar array type for field '" << base_name << "'.");
//                 }
//                 continue;
//             }

//             // --------------------------
//             // ITERABLE FIELDS
//             // --------------------------
//             read_iterable_field(data_map, val_array, shape_array, offset_array, num_rows, base_name);
//         }
//     }

//     virtual void populate_final_map() {};

// private:
//     // Helper: read iterable fields
//     void read_iterable_field(
//         TMap &data_map,
//         vtkAbstractArray *val_array,
//         vtkIntArray *shape_array,
//         vtkIntArray *offset_array,
//         vtkIdType num_rows,
//         const std::string &base_name)
//     {
//         const int *shape_data = shape_array->GetPointer(0);
//         const int *offset_data = offset_array->GetPointer(0);

//         // Integers
//         if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
//         {
//             const int *vals = int_arr->GetPointer(0);
//             vtkIdType cursor = 0;
//             for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//             {
//                 int shape_start = offset_data[irow];
//                 int ndim = shape_data[shape_start];

//                 std::vector<int> dims(ndim);
//                 vtkIdType total_elems = 1;
//                 for (int j = 0; j < ndim; ++j)
//                 {
//                     dims[j] = shape_data[shape_start + 1 + j];
//                     total_elems *= dims[j];
//                 }

//                 std::vector<int> flat(vals + cursor, vals + cursor + total_elems);
//                 cursor += total_elems;
//                 // FIX: wrap flat in a vector of vectors
//                 data_map[irow].emplace_back(std::vector<std::vector<int>>{flat});
//             }
//         }
//         // Doubles
//         else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
//         {
//             const double *vals = dbl_arr->GetPointer(0);
//             vtkIdType cursor = 0;
//             for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//             {
//                 int shape_start = offset_data[irow];
//                 int ndim = shape_data[shape_start];

//                 std::vector<int> dims(ndim);
//                 vtkIdType total_elems = 1;
//                 for (int j = 0; j < ndim; ++j)
//                 {
//                     dims[j] = shape_data[shape_start + 1 + j];
//                     total_elems *= dims[j];
//                 }

//                 std::vector<double> flat(vals + cursor, vals + cursor + total_elems);
//                 cursor += total_elems;
//                 // FIX: wrap flat in a vector of vectors
//                 data_map[irow].emplace_back(std::vector<std::vector<Real>>{flat});
//             }
//         }
//         // Strings
//         else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
//         {
//             vtkIdType cursor = 0;
//             for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
//             {
//                 int shape_start = offset_data[irow];
//                 int ndim = shape_data[shape_start];

//                 std::vector<int> dims(ndim);
//                 vtkIdType total_elems = 1;
//                 for (int j = 0; j < ndim; ++j)
//                 {
//                     dims[j] = shape_data[shape_start + 1 + j];
//                     total_elems *= dims[j];
//                 }

//                 std::vector<std::string> flat;
//                 flat.reserve(total_elems);
//                 for (vtkIdType k = 0; k < total_elems; ++k)
//                     flat.emplace_back(str_arr->GetValue(cursor++));
//                 // FIX: wrap flat in a vector of vectors
//                 data_map[irow].emplace_back(std::vector<std::vector<std::string>>{flat});
//             }
//         }
//         else
//         {
//             libmesh_error_msg("Unsupported iterable array type for field '" + base_name + "'.");
//         }
//     }

//     // Helper: flatten reshaping (currently unused but can be extended)
//     template <typename T>
//     static std::vector<std::vector<T>> reshape(
//         const std::vector<T> &flat,
//         const std::vector<int> &dims,
//         size_t level,
//         size_t offset)
//     {
//         if (dims.empty())
//             return {};

//         if (level == dims.size() - 1)
//         {
//             std::vector<std::vector<T>> leaf;
//             size_t chunk = flat.size() / dims[level];
//             for (int i = 0; i < dims[level]; ++i)
//                 leaf.emplace_back(flat.begin() + offset + i * chunk,
//                                   flat.begin() + offset + (i + 1) * chunk);
//             return leaf;
//         }

//         size_t stride = flat.size() / dims[level];
//         std::vector<std::vector<T>> result;
//         for (int i = 0; i < dims[level]; ++i)
//         {
//             auto sub = reshape(flat, dims, level + 1, offset + i * stride);
//             result.insert(result.end(), sub.begin(), sub.end());
//         }
//         return result;
//     }
// };

#endif // FE_DATA_OBJECT_H