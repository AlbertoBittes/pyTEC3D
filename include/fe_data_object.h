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

using namespace libMesh;

class LoadStep;

typedef std::vector<std::vector<std::variant<int, Real, std::string>>> DataframeMap;
typedef std::unordered_map<std::string, std::variant<int, Real, std::string>> SettingsMap;

// FE data object classes are the basic data repositories of the program.
// They store information into a map (dict) that can be later retrieved when needed
// Each type of data e.g. Materials, Real Constants, Solver Settings, Time Steps, etc. is stored into a unordered map.
// The basic classes can be enriched by peculiar methods (e.g. the LoadsObject has methods to interface with the
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

typedef DataframeMapFEDataObject<Real> MaterialsMap;
typedef DataframeMapFEDataObject<Real> RealConstantsMap;
typedef DataframeMapFEDataObject<int> ETypeKeyoptsMap;
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