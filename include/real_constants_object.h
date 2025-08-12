#ifndef REAL_CONSTANTS_OBJECT_H
#define REAL_CONSTANTS_OBJECT_H

#include <unordered_map>
#include <vector>
#include <memory>
#include "fe_data_object.h" // Include the base class

// Holds the RealConstants numbers (R in APDL) and their values to be used by ElementBuilder during
// element matrices' construction
class RealConstantsObject : public DataframeMapFEDataObject<Real>
{
public:
    // Default constructor
    using DataframeMapFEDataObject<Real>::DataframeMapFEDataObject;

    // Check if a constant exists in the map
    bool has_constant(dof_id_type constant_id) const { return dataframe_map->find(constant_id) != dataframe_map->end(); }

    // Get the values of a constant by ID
    const std::vector<Real> &get_constant(dof_id_type constant_id) const
    {
        auto it = dataframe_map->find(constant_id);
        if (it == dataframe_map->end())
        {
            libmesh_error_msg("Constant not found: " + std::to_string(constant_id));
        }
        return it->second;
    }

    // // Set the values of a constant
    // void set_constant(dof_id_type constant_id, const std::vector<Real> &values)
    // {
    //     (*dataframe_map)[constant_id] = values;
    // }
};

#endif // REAL_CONSTANTS_OBJECT_H