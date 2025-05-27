#ifndef REAL_CONSTANTS_OBJECT_H
#define REAL_CONSTANTS_OBJECT_H

#include "fe_data_object.h" // Include the base class
#include <unordered_map>
#include <vector>
#include <memory>

class RealConstantsObject : public DataframeMatLikeFEDataObject
{
public:
    // Default constructor
    using DataframeMatLikeFEDataObject::DataframeMatLikeFEDataObject;

    // Check if a constant exists in the map
    bool has_constant(dof_id_type constant_id) const { return mat_like_map->find(constant_id) != mat_like_map->end(); }

    // Get the values of a constant by ID
    const std::vector<Real> &get_constant(dof_id_type constant_id) const
    {
        auto it = mat_like_map->find(constant_id);
        if (it == mat_like_map->end())
        {
            libmesh_error_msg("Constant not found: " + std::to_string(constant_id));
        }
        return it->second;
    }

    // // Set the values of a constant
    // void set_constant(dof_id_type constant_id, const std::vector<Real> &values)
    // {
    //     (*mat_like_map)[constant_id] = values;
    // }
};

#endif // REAL_CONSTANTS_OBJECT_H