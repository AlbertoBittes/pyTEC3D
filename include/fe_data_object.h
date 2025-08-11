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
typedef std::unordered_map<dof_id_type, std::vector<Real>> MaterialsMap;
typedef std::unordered_map<dof_id_type, std::vector<Real>> RealConstantsMap;
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

    // NEW: Constructor that receives an external map and a file, sets the map and calls read_data
    FEDataObjectMap(std::ifstream &file, std::shared_ptr<T1> external_map)
        : data_map(external_map)
    {
        this->read_data(file);
    }

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

class DataframeMatLikeFEDataObject : public DataframeFEDataObject
{
public:
    // Constructor
    DataframeMatLikeFEDataObject(std::ifstream &file)
        : DataframeFEDataObject(file), // parent constructor
          mat_like_map(std::make_shared<MaterialsMap>())
    {
        populate_mat_like_map();
    }

    // Method to populate mat_like_map from data_map
    void populate_mat_like_map();

    // Getter for mat_like_map
    MaterialsMap *get_mat_like_map() const
    {
        return mat_like_map.get();
    }

protected:
    std::shared_ptr<MaterialsMap> mat_like_map; // Shared pointer to the map
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

class InfiniteElementsOptionsObject : public SettingsFEDataObject
{
public:
    InfiniteElementsOptionsObject() = default;

    // Populate the settings map from a vtkFieldData pointer
    void populate_from_vtk_field_data(vtkFieldData *field_data)
    {
        // Set defaults
        this->set_setting("order", std::string("FIRST"));
        this->set_setting("mesh_dimension", 2);
        this->set_setting("edge_elements", 0);
        this->set_setting("boundary_type", std::string("SPHERICAL"));
        this->set_setting("origin_x", Real(0));
        this->set_setting("origin_y", Real(0));
        this->set_setting("origin_z", Real(0));

        if (!field_data)
            return;

        // Helper lambdas
        auto get_string = [&](const char *name) -> std::optional<std::string>
        {
            if (auto arr = vtkStringArray::SafeDownCast(field_data->GetAbstractArray(name)))
                if (arr->GetNumberOfTuples() > 0)
                    return arr->GetValue(0);
            return std::nullopt;
        };
        auto get_int = [&](const char *name) -> std::optional<int>
        {
            if (auto arr = vtkIntArray::SafeDownCast(field_data->GetArray(name)))
                if (arr->GetNumberOfTuples() > 0)
                    return arr->GetValue(0);
            return std::nullopt;
        };
        auto get_real = [&](const char *name) -> std::optional<Real>
        {
            if (auto arr = vtkDoubleArray::SafeDownCast(field_data->GetArray(name)))
                if (arr->GetNumberOfTuples() > 0)
                    return static_cast<Real>(arr->GetValue(0));
            return std::nullopt;
        };

        // Try to read each field, fallback to default if not found
        if (auto v = get_string("order"))
            this->set_setting("order", *v);
        if (auto v = get_int("mesh_dimension"))
            this->set_setting("mesh_dimension", *v);
        if (auto v = get_int("edge_elements"))
            this->set_setting("edge_elements", *v);
        if (auto v = get_string("boundary_type"))
            this->set_setting("boundary_type", *v);
        if (auto v = get_real("origin_x"))
            this->set_setting("origin_x", *v);
        if (auto v = get_real("origin_y"))
            this->set_setting("origin_y", *v);
        if (auto v = get_real("origin_z"))
            this->set_setting("origin_z", *v);
    }
};

#endif // FE_DATA_OBJECT_H