#include "../include/materials_object.h"
#include <iostream>
#include <stdexcept>

void MaterialsObject::readMaterialsProperties(std::istream &file_stream)
{
    libMesh::out << "Reading materials properties from binary stream..." << std::endl;
    auto new_map = std::make_shared<MaterialsMap>();

    // Read the number of materials
    unsigned int num_materials;
    file_stream.read(reinterpret_cast<char *>(&num_materials), sizeof(num_materials));
    if (!file_stream)
    {
        libmesh_error_msg("Failed to read the number of materials from the binary stream.");
    }

    // Read the number of values per material
    unsigned int num_values_per_material;
    file_stream.read(reinterpret_cast<char *>(&num_values_per_material), sizeof(num_values_per_material));
    if (!file_stream)
    {
        libmesh_error_msg("Failed to read the number of values per material from the binary stream.");
    }

    // Read each material
    for (unsigned int i = 0; i < num_materials; ++i)
    {
        unsigned int material_id;
        std::vector<Real> material_properties(num_values_per_material);

        // Read the material ID
        file_stream.read(reinterpret_cast<char *>(&material_id), sizeof(material_id));
        if (!file_stream)
        {
            libmesh_error_msg("Failed to read material ID from the binary stream.");
        }

        // Read the material properties
        file_stream.read(reinterpret_cast<char *>(material_properties.data()), num_values_per_material * sizeof(Real));
        if (!file_stream)
        {
            libmesh_error_msg("Failed to read material properties from the binary stream.");
        }

        // Add the material to the map
        (*new_map)[material_id] = std::move(material_properties);
    }

    // Set the new MaterialsMap
    setMaterialsMap(new_map);
}