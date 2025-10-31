#ifndef FE_DATA_OBJECT_H
#define FE_DATA_OBJECT_H

#include <cstdint>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>

#include <vtkAbstractArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include "globals.h"

#include "globals.h"

using namespace libMesh;

class LoadStep;

typedef std::vector<std::vector<std::variant<int, Real, std::string>>> DataframeMap;

typedef std::variant<int, Real, std::string, Globals::NDArray> SettingsValue;
typedef std::unordered_map<std::string, SettingsValue> SettingsMap;

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
    T get_scalar_setting(const std::string &key) const
    {
        auto it = data_map->find(key);
        if (it == data_map->end())
        {
            libmesh_error_msg("Setting not found: " + key);
        }

        if (std::holds_alternative<T>(it->second))
            return std::get<T>(it->second);

        if (std::holds_alternative<Globals::NDArray>(it->second))
            return extract_scalar_from_ndarray<T>(key, std::get<Globals::NDArray>(it->second));

        libmesh_error_msg("Setting '" + key + "' cannot be converted to requested scalar type.");
        return T{};
    }

    template <typename T>
    void set_setting(const std::string &key, T value)
    {
        (*data_map)[key] = value;
    }
    // Implement the read_data method
    void read_data(std::ifstream &file);

private:
    template <typename T>
    static T extract_scalar_from_ndarray(const std::string &key, const Globals::NDArray &array)
    {
        std::size_t total_elements = 1;
        for (auto extent : array.shape)
            total_elements *= static_cast<std::size_t>(extent);

        if (total_elements == 0)
            libmesh_error_msg("Setting '" + key + "' is empty; cannot convert to scalar.");
        if (total_elements != 1)
            libmesh_error_msg("Setting '" + key + "' contains multiple values; cannot convert to scalar.");

        return std::visit(
            [&](const auto &vec) -> T
            {
                using VecType = std::decay_t<decltype(vec)>;
                using ElemType = typename VecType::value_type;

                if (vec.size() != total_elements || vec.empty())
                    libmesh_error_msg("NDArray payload for setting '" + key + "' has unexpected size.");

                const auto &value = vec.front();

                if constexpr (std::is_same_v<T, int>)
                {
                    if constexpr (std::is_same_v<ElemType, int>)
                        return value;
                    if constexpr (std::is_same_v<ElemType, Real>)
                        return static_cast<int>(value);
                }
                else if constexpr (std::is_same_v<T, Real>)
                {
                    if constexpr (std::is_same_v<ElemType, Real>)
                        return value;
                    if constexpr (std::is_same_v<ElemType, int>)
                        return static_cast<Real>(value);
                }
                else if constexpr (std::is_same_v<T, std::string>)
                {
                    if constexpr (std::is_same_v<ElemType, std::string>)
                        return value;
                }

                libmesh_error_msg("Unsupported NDArray to scalar conversion for setting '" + key + "'.");
                return T{};
            },
            array.values);
    }
};

template class DataframeMapFEDataObject<int>;
template class DataframeMapFEDataObject<Real>;

#endif // FE_DATA_OBJECT_H
