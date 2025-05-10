#ifndef MATERIALS_OBJECT_H
#define MATERIALS_OBJECT_H

#include "../include/fe_data_object.h" // Include the base class
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <istream>
#include <libmesh/libmesh.h>

using namespace libMesh;

// Forward declaration of LoadStep
class LoadStep;

// Define MaterialsMap as a type alias
typedef std::unordered_map<dof_id_type, std::vector<Real>> MaterialsMap;

class MaterialsObject : public DataframeMatLikeFEDataObject
{
public:
    // Default constructor
    explicit MaterialsObject(std::shared_ptr<LoadStep> parent = nullptr,
                             std::shared_ptr<MaterialsMap> map = nullptr)
        : DataframeMatLikeFEDataObject(parent, map),
          mu_matrix(3, 0.0),
          rs_matrix(3, 0.0) {}

    // Check if a material ID exists in the map
    bool has_material(dof_id_type mat_id) const { return mat_like_map->find(mat_id) != mat_like_map->end(); }

    // Get material properties by ID
    const std::vector<Real> &get_material_properties(dof_id_type mat_id) const
    {
        return mat_like_map->at(mat_id);
    }

    // Fill the mu_matrix for a given material ID
    inline void fill_mu_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = mat_like_map->at(mat_id);

        mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties[0];
        mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties[1];
        mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties[2];
    }

    // Fill the rs_matrix for a given material ID
    inline void fill_rs_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = mat_like_map->at(mat_id);

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