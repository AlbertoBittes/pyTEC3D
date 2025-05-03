#include "../include/loads_object.h"

void LoadsObject::read_data(std::ifstream &file)
{
    libMesh::out << "Reading loads data from an open file stream..." << std::endl;

    // Initialize the data_map if it is a nullptr
    if (!data_map)
    {
        data_map = std::make_shared<LoadArraysMap>();
    }

    // Read the number of load instances
    unsigned int num_loads;
    file.read(reinterpret_cast<char *>(&num_loads), sizeof(num_loads));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of loads from the binary stream.");
    }

    // Define the vector of load labels for mapping
    std::vector<std::string> loads_labels = {"Js"}; // Example load labels
    std::string load_label;

    // Read the loads data
    for (unsigned int i = 0; i < num_loads; ++i)
    {
        dof_id_type element_number;
        unsigned int load_label_number;
        unsigned int component_number;
        Real load_value;

        // Read element number
        file.read(reinterpret_cast<char *>(&element_number), sizeof(element_number));

        // Read load label number
        file.read(reinterpret_cast<char *>(&load_label_number), sizeof(load_label_number));

        // Read component number
        file.read(reinterpret_cast<char *>(&component_number), sizeof(component_number));

        // Read load value
        file.read(reinterpret_cast<char *>(&load_value), sizeof(load_value));

        // Validate the load label number
        if (load_label_number >= loads_labels.size())
        {
            libmesh_error_msg("Error: Invalid load label number " + std::to_string(load_label_number));
        }

        // Map the load label number to the corresponding string
        load_label = loads_labels[load_label_number];

        // Access the ElementLoad map for the given load label
        auto &element_load = (*data_map)[load_label];

        // Check if the element exists in the map
        auto element_it = element_load.find(element_number);
        if (element_it != element_load.end())
        {
            // Ensure the component number is valid
            if (component_number < element_it->second.size())
            {
                // Update the specific component with the load value
                element_it->second[component_number] = load_value;
            }
            else
            {
                libmesh_error_msg("Warning: Invalid component number " +
                                  std::to_string(component_number) + " for element " +
                                  std::to_string(element_number) + " in load label " +
                                  load_label);
            }
        }
    }
}

void LoadsObject::zero_loads_array()
{
    if (!data_map)
    {
        return;
    }
    // Iterate over each load in the LoadArraysMap
    for (auto &load_array_pair : *data_map)
    {
        ElementLoad &element_load = load_array_pair.second; // Access the ElementLoad map

        // Iterate over each element in the ElementLoad map
        for (auto &element_pair : element_load)
        {
            std::vector<Real> &values = element_pair.second; // Access the vector of Real values

            // Zero out the vector using std::fill for efficiency
            if (!values.empty())
            {
                std::fill(values.begin(), values.end(), 0.0);
            }
        }
    }
}

void LoadsObject::build_load_data_vectors(const MeshBase &mesh, std::set<std::string> &loads_labels)
{
    // Ensure data_map is initialized
    if (!data_map)
    {
        data_map = std::make_shared<LoadArraysMap>();
    }

    // Get the load quantity from the cell data
    for (const auto &label : loads_labels)
    {
        // If not, create a new load vector for this label
        (*data_map)[label] = ElementLoad();
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            // If not, create a new load for this element ID
            (*data_map)[label][elem->id()] = std::vector<Real>(3, 0.0);
        }
    }

    if (!oldLoadsVector)
    {
        oldLoadsVector = std::make_shared<LoadArraysMap>();
    }
    *oldLoadsVector = *data_map;
}
