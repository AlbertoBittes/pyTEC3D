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

#include "globals.h"

using namespace libMesh;

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

template <typename>
struct AlwaysFalse : std::false_type
{
};

struct BinaryColumnarData
{
    std::size_t entry_count = 0;
    std::unordered_map<std::string, Globals::NDArray> columns;
};

class BinaryDataFEDataObject : public FEDataObjectMap<BinaryColumnarData>
{
public:
    BinaryDataFEDataObject()
        : FEDataObjectMap<BinaryColumnarData>(std::make_shared<BinaryColumnarData>())
    {
    }

    explicit BinaryDataFEDataObject(std::shared_ptr<BinaryColumnarData> external_map)
        : FEDataObjectMap<BinaryColumnarData>(external_map ? external_map : std::make_shared<BinaryColumnarData>())
    {
    }

    explicit BinaryDataFEDataObject(std::ifstream &file)
        : FEDataObjectMap<BinaryColumnarData>(std::make_shared<BinaryColumnarData>())
    {
        this->read_data(file);
    }

    template <typename T>
    T get_scalar(const std::string &key, std::size_t index = 0) const;

    int get_int(const std::string &key, std::size_t index = 0) const;
    Real get_real(const std::string &key, std::size_t index = 0) const;
    std::string get_string(const std::string &key, std::size_t index = 0) const;

    const std::vector<int> &get_int_column(const std::string &key) const;
    const std::vector<Real> &get_real_column(const std::string &key) const;
    const std::vector<std::string> &get_string_column(const std::string &key) const;

    bool has_column(const std::string &key) const;
    std::size_t entry_count() const;
    bool empty() const { return entry_count() == 0; }

    void read_data(std::ifstream &file);

protected:
    const Globals::NDArray &get_column(const std::string &key) const;
    static std::size_t flattened_size(const Globals::NDArray &array);
    static void ensure_index_in_bounds(const std::string &key, std::size_t index, std::size_t size);
};

template <typename T>
T BinaryDataFEDataObject::get_scalar(const std::string &key, std::size_t index) const
{
    if constexpr (std::is_same_v<T, int>)
        return static_cast<T>(get_int(key, index));
    else if constexpr (std::is_same_v<T, Real>)
        return static_cast<T>(get_real(key, index));
    else if constexpr (std::is_same_v<T, std::string>)
        return get_string(key, index);
    else
    {
        static_assert(AlwaysFalse<T>::value, "Unsupported scalar type requested from BinaryDataFEDataObject.");
    }
}

#endif // FE_DATA_OBJECT_H
