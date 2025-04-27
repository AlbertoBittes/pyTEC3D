#ifndef MATERIALS_OBJECT_H
#define MATERIALS_OBJECT_H

#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <istream>
#include <libmesh/libmesh.h>

using namespace libMesh;

// Forward declaration of LoadStep
class LoadStep;

// // Structure to hold material properties
// struct MaterialsProperties
// {
//     Real murx, mury, murz;
//     Real rsvx, rsvy, rsvz;
//     Real perx, pery, perz;
// };

typedef std::unordered_map<dof_id_type, std::vector<Real>> MaterialsMap;

class MaterialsObject
{
public:
    // Default constructor
    explicit MaterialsObject(LoadStep *parent = nullptr)
        : parent_load_step(parent),
          materials_map(std::make_shared<MaterialsMap>()),
          mu_matrix(3, 0.0),
          rs_matrix(3, 0.0) {}

    // Constructor with an existing MaterialsMap
    explicit MaterialsObject(std::shared_ptr<MaterialsMap> map, LoadStep *parent = nullptr)
        : parent_load_step(parent),
          materials_map(std::move(map)),
          mu_matrix(3, 0.0),
          rs_matrix(3, 0.0) {}

    // Getter for the MaterialsMap
    const MaterialsMap &getMaterialsMap() const { return *materials_map; }

    // Setter for the MaterialsMap
    void setMaterialsMap(std::shared_ptr<MaterialsMap> new_map) { materials_map = std::move(new_map); }

    // Check if a material ID exists in the map
    bool hasMaterial(dof_id_type mat_id) const { return materials_map->find(mat_id) != materials_map->end(); }

    // Get material properties by ID
    const std::vector<Real> &getMaterialProperties(dof_id_type mat_id) const
    {
        return materials_map->at(mat_id);
    }

    // Method to read material properties from an already open binary file stream
    void readMaterialsProperties(std::istream &file_stream);

    // Fill the mu_matrix for a given material ID
    inline void fill_mu_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = materials_map->at(mat_id);

        mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties[0];
        mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties[1];
        mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties[2];
    }

    // Fill the rs_matrix for a given material ID
    inline void fill_rs_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = materials_map->at(mat_id);

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
    const std::vector<Real> &getMuMatrix() const { return mu_matrix; }

    // Get the rs_matrix
    const std::vector<Real> &getRsMatrix() const { return rs_matrix; }

    // Set the parent LoadStep
    void setParentLoadStep(LoadStep *parent) { parent_load_step = parent; }

    // Get the parent LoadStep
    LoadStep *getParentLoadStep() const { return parent_load_step; }

private:
    LoadStep *parent_load_step;                  // Pointer to the parent LoadStep
    std::shared_ptr<MaterialsMap> materials_map; // Holds the MaterialsMap
    std::vector<Real> mu_matrix;                 // Magnetic permeability matrix
    std::vector<Real> rs_matrix;                 // Resistivity matrix
};

#endif // MATERIALS_OBJECT_H