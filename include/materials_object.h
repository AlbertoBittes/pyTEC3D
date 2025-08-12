#ifndef MATERIALS_OBJECT_H
#define MATERIALS_OBJECT_H

#include <libmesh/libmesh.h>
#include "fe_data_object.h" // Include the base class
#include "globals.h"

using namespace libMesh;

// Forward declaration of LoadStep
class LoadStep;

// Materials holds the materials numbers and their properties (mu, resistivity, etc.) (MP in APDL)
// And evaluates the
class MaterialsObject : public MaterialsMap
{
public:
    // Default constructor
    explicit MaterialsObject(std::ifstream &file)
        : MaterialsMap(file),
          mu_matrix(3, 0.0),
          rs_matrix(3, 0.0) {}

    // Check if a material ID exists in the map
    bool has_material(dof_id_type mat_id) const { return dataframe_map->find(mat_id) != dataframe_map->end(); }

    // Get material properties by ID
    const std::vector<Real> &get_material_properties(dof_id_type mat_id) const
    {
        return dataframe_map->at(mat_id);
    }

    // Fill the mu_matrix for a given material ID
    inline void fill_mu_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = dataframe_map->at(mat_id);

        mu_matrix[0] = Globals::MU_0 * mat_properties[0];
        mu_matrix[1] = Globals::MU_0 * mat_properties[1];
        mu_matrix[2] = Globals::MU_0 * mat_properties[2];
    }

    // Fill the rs_matrix for a given material ID
    inline void fill_rs_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = dataframe_map->at(mat_id);

        rs_matrix[0] = mat_properties[3];
        rs_matrix[1] = mat_properties[4];
        rs_matrix[2] = mat_properties[5];
    }

    // Compute material properties for a given material ID
    inline void compute_material_properties(dof_id_type mat_id, bool conductive)
    {
        fill_mu_matrix(mat_id);
        if (conductive)
            fill_rs_matrix(mat_id);
    }

    // Get the mu_matrix
    const std::vector<Real> &get_mu_matrix() const { return mu_matrix; }

    // Get the rs_matrix
    const std::vector<Real> &get_rs_matrix() const { return rs_matrix; }

private:
    std::vector<Real> mu_matrix; // Magnetic permeability matrix
    std::vector<Real> rs_matrix; // Resistivity matrix
};

#endif // MATERIALS_OBJECT_H