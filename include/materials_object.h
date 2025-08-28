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
class MaterialsObject : public DataframeMapFEDataObject<Real>
{
public:
    // Default constructor
    explicit MaterialsObject(std::ifstream &file)
        : DataframeMapFEDataObject<Real>(file),
          nu_matrix(3, 0.0),
          sigma_matrix(3, 0.0) {}

    // Check if a material ID exists in the map
    bool has_material(dof_id_type mat_id) const { return dataframe_map->find(mat_id) != dataframe_map->end(); }

    // Get material properties by ID
    const std::vector<Real> &get_material_properties(dof_id_type mat_id) const
    {
        return dataframe_map->at(mat_id);
    }

    // Fill the nu_matrix for a given material ID
    inline void fill_nu_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = dataframe_map->at(mat_id);

        nu_matrix[0] = 1 / Globals::MU_0 * mat_properties[0];
        nu_matrix[1] = 1 / Globals::MU_0 * mat_properties[1];
        nu_matrix[2] = 1 / Globals::MU_0 * mat_properties[2];
    }

    // Fill the sigma_matrix for a given material ID
    inline void fill_sigma_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = dataframe_map->at(mat_id);

        sigma_matrix[0] = 1 / mat_properties[3];
        sigma_matrix[1] = 1 / mat_properties[4];
        sigma_matrix[2] = 1 / mat_properties[5];
    }

    // Compute material properties for a given material ID
    inline void compute_material_properties(dof_id_type mat_id, bool conductive)
    {
        fill_nu_matrix(mat_id);
        if (conductive)
            fill_sigma_matrix(mat_id);
    }

    // Get the nu_matrix
    const std::vector<Real> &get_nu_matrix() const { return nu_matrix; }

    // Get the sigma_matrix
    const std::vector<Real> &get_sigma_matrix() const { return sigma_matrix; }

private:
    std::vector<Real> nu_matrix;    // Magnetic permeability matrix
    std::vector<Real> sigma_matrix; // Resistivity matrix
};

#endif // MATERIALS_OBJECT_H