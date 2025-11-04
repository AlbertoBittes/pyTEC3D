#ifndef MATERIALS_OBJECT_H
#define MATERIALS_OBJECT_H

#include <array>
#include <unordered_map>
#include <vector>
#include <libmesh/libmesh.h>
#include "fe_data_object.h" // Include the base class
#include "globals.h"

using namespace libMesh;

// Materials holds the materials numbers and their properties (mu, resistivity, etc.) (MP in APDL)
// And evaluates the
class MaterialsObject : public BinaryDataFEDataObject
{
public:
    // Default constructor
    explicit MaterialsObject(std::ifstream &file)
        : BinaryDataFEDataObject(file),
          nu_matrix(3, 0.0),
          sigma_matrix(3, 0.0)
    {
        build_materials_map();
    }

    void read_data(std::ifstream &file) override
    {
        BinaryDataFEDataObject::read_data(file);
        build_materials_map();
    }

    // Check if a material ID exists in the map
    bool has_material(dof_id_type mat_id) const { return materials_map.find(mat_id) != materials_map.end(); }

    // Get material properties by ID
    const std::vector<Real> &get_material_properties(dof_id_type mat_id) const
    {
        auto it = materials_map.find(mat_id);
        if (it == materials_map.end())
            libmesh_error_msg("MaterialsObject: material not found with id " + std::to_string(mat_id));
        return it->second;
    }

    // Fill the nu_matrix for a given material ID
    inline void fill_nu_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = get_material_properties(mat_id);

        nu_matrix[0] = 1 / (Globals::MU_0 * mat_properties[0]);
        nu_matrix[1] = 1 / (Globals::MU_0 * mat_properties[1]);
        nu_matrix[2] = 1 / (Globals::MU_0 * mat_properties[2]);
    }

    // Fill the sigma_matrix for a given material ID
    inline void fill_sigma_matrix(dof_id_type mat_id)
    {
        const std::vector<Real> &mat_properties = get_material_properties(mat_id);

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

    const std::unordered_map<dof_id_type, std::vector<Real>> &get_materials_map() const { return materials_map; }

private:
    void build_materials_map()
    {
        materials_map.clear();
        if (empty())
            return;

        static constexpr std::array<const char *, 12> PROPERTY_COLUMNS = {
            "MURX", "MURY", "MURZ", "RSVX", "RSVY", "RSVZ",
            "PERX", "PERY", "PERZ", "MGXX", "MGYY", "MGZZ"};

        if (!has_column("ansys_material_type"))
            libmesh_error_msg("MaterialsObject: missing 'ansys_material_type' column.");

        const auto &ids = get_int_column("ansys_material_type");
        const std::size_t count = entry_count();

        if (ids.size() < count)
            libmesh_error_msg("MaterialsObject: inconsistent 'ansys_material_type' column length.");

        for (std::size_t row = 0; row < count; ++row)
        {
            std::vector<Real> properties;
            properties.reserve(PROPERTY_COLUMNS.size());
            for (const auto *name : PROPERTY_COLUMNS)
            {
                if (!has_column(name))
                    libmesh_error_msg(std::string("MaterialsObject: missing column '") + name + "'.");
                const auto &col = get_real_column(name);
                if (row >= col.size())
                    libmesh_error_msg(std::string("MaterialsObject: column '") + name + "' has insufficient entries.");
                properties.push_back(col[row]);
            }
            materials_map[ids[row]] = std::move(properties);
        }
    }

    std::unordered_map<dof_id_type, std::vector<Real>> materials_map;
    std::vector<Real> nu_matrix;    // Magnetic permeability matrix
    std::vector<Real> sigma_matrix; // Resistivity matrix
};

#endif // MATERIALS_OBJECT_H
