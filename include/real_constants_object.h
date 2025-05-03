#ifndef REAL_CONSTANTS_OBJECT_H
#define REAL_CONSTANTS_OBJECT_H

#include "../include/fe_data_object.h" // Include the base class
#include <unordered_map>
#include <vector>
#include <memory>
#include <istream>
#include <libmesh/libmesh.h>

using namespace libMesh;

// Define RealConstantsMap as a type alias
typedef std::unordered_map<dof_id_type, std::vector<Real>> RealConstantsMap;

class RealConstantsObject : public CSVFEDataObject
{
public:
    // Default constructor
    explicit RealConstantsObject(std::shared_ptr<LoadStep> parent = nullptr,
                                 std::shared_ptr<RealConstantsMap> map = nullptr)
        : CSVFEDataObject(parent, map) {}

    // Check if a constant exists in the map
    bool has_constant(dof_id_type constant_id) const { return data_map->find(constant_id) != data_map->end(); }

    // Get the values of a constant by ID
    const std::vector<Real> &get_constant(dof_id_type constant_id) const
    {
        auto it = data_map->find(constant_id);
        if (it == data_map->end())
        {
            throw std::runtime_error("Constant not found: " + std::to_string(constant_id));
        }
        return it->second;
    }

    // Set the values of a constant
    void set_constant(dof_id_type constant_id, const std::vector<Real> &values)
    {
        (*data_map)[constant_id] = values;
    }
};

#endif // REAL_CONSTANTS_OBJECT_H