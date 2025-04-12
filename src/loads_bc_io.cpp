#include "../include/loads_bc_io.h"
#include <fstream>
#include <memory>

void read_boundary_conditions(std::string &file_name, BoundaryConditionsMap &node_conditions)
{
    libMesh::out << "Reading boundary conditions data..." << std::endl;
    // Open the file for reading in binary mode
    std::ifstream file(file_name, std::ios::binary);
    if (!file.is_open())
    {
        libmesh_error_msg("Error: Could not open file " + file_name);
    }

    while (file)
    {
        dof_id_type node;
        uint8_t label_length;
        Real value;

        // Read node number
        file.read(reinterpret_cast<char *>(&node), sizeof(node));
        if (!file)
            break;

        // Read label length
        file.read(reinterpret_cast<char *>(&label_length), sizeof(label_length));
        if (!file)
            break;

        // Read label
        std::vector<char> label(label_length);
        file.read(label.data(), label_length);
        if (!file)
            break;

        // Read value
        file.read(reinterpret_cast<char *>(&value), sizeof(value));
        if (!file)
            break;

        // Store the boundary condition if the node exists in the map
        auto it = node_conditions.find(node);
        if (it != node_conditions.end())
        {
            std::string label_str(label.data(), label_length); // Convert char array to string
            it->second.emplace(std::move(label_str), value);
        }
    }
    // Close the file when done
    file.close();
}

void read_loads_arrays(std::string &file_name, std::unique_ptr<LoadArraysMap> loads_arrays)
{
    libMesh::out << "Reading loads data..." << std::endl;
    // Open the file for reading in binary mode
    std::ifstream file(file_name, std::ios::binary);
    if (!file.is_open())
    {
        libmesh_error_msg("Error: Could not open file " + file_name);
    }
    while (file)
    {
        uint8_t name_length;
        unsigned int num_components;
        unsigned int num_elements;

        // Read name length
        file.read(reinterpret_cast<char *>(&name_length), sizeof(name_length));
        if (!file)
            break;

        // Read name
        std::vector<char> name(name_length);
        file.read(name.data(), name_length);
        if (!file)
            break;

        // Read number of components
        file.read(reinterpret_cast<char *>(&num_components), sizeof(num_components));
        if (!file)
            break;

        // Read number of elements
        file.read(reinterpret_cast<char *>(&num_elements), sizeof(num_elements));
        if (!file)
            break;

        std::string variable_name(name.data(), name_length);
        ElementLoad &load_arrays = (*loads_arrays)[variable_name];
        // load_arrays.reserve(num_elements); // Preallocate memory

        for (unsigned int i = 0; i < num_elements; ++i)
        {
            dof_id_type element;
            std::vector<Real> values(num_components);

            // Read element number
            file.read(reinterpret_cast<char *>(&element), sizeof(element));
            if (!file)
                break;

            // Read values
            file.read(reinterpret_cast<char *>(values.data()), num_components * sizeof(Real));
            if (!file)
                break;

            // Store the load array if the element exists in the map
            auto it = load_arrays.find(element);
            if (it != load_arrays.end())
            {
                it->second = std::move(values);
            }
        }
    }
    // Close the file when done
    file.close();
}