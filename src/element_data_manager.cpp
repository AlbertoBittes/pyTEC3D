#include "../include/element_data_manager.h"

void ElementDataManager::read_loads_arrays(std::string &file_name)
{
    libMesh::out << "Reading loads data..." << std::endl;

    // Open the file for reading in binary mode
    std::ifstream file(file_name, std::ios::binary);
    if (!file.is_open())
    {
        libmesh_error_msg("Error: Could not open file " + file_name);
    }

    // Define the vector of load labels for mapping
    std::vector<std::string> loads_labels = {"Js", "Jm", "Jc"}; // Example load labels
    std::string load_label;
    // Read the loads data
    while (file.peek() != EOF)
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
        auto &element_load = (*activeLoadsVector)[load_label];

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
                libMesh::out << "Warning: Invalid component number " << component_number
                             << " for element " << element_number << std::endl;
            }
        }
        else
        {
            // Element ID not found, skip this entry
            libMesh::out << "Warning: Element ID " << element_number
                         << " not found in activeLoadsVector for load label " << load_label << std::endl;
        }
    }

    file.close();
}

void ElementDataManager::read_update_loads_vectors(std::string file_name)
{
    this->stashLoadsVector();
    // Read the loads data from the file
    this->read_loads_arrays(file_name);
    // Update the old load vector with the new load vector
}

void ElementDataManager::zero_load_arrays()
{
    if (!activeLoadsVector)
    {
        return;
    }
    // Iterate over each load in the LoadArraysMap
    for (auto &load_array_pair : *activeLoadsVector)
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

void ElementDataManager::build_load_data_vectors(const MeshBase &mesh, std::set<std::string> loads_labels)
{
    // Get the load quantity from the cell data
    for (const auto &label : loads_labels)
    {
        // If not, create a new load vector for this label
        (*activeLoadsVector)[label] = ElementLoad();
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            // If not, create a new load for this element ID
            (*activeLoadsVector)[label][elem->id()] = std::vector<Real>(3, 0.0);
        }
    }
    oldLoadsVector = std::make_unique<LoadArraysMap>(*activeLoadsVector);
}

void ElementDataManager::set_load_quantity(std::string load_name)
{
    // Get the load quantity from the cell data
    load_array = &(*activeLoadsVector)[load_name];
    load_array = &(*oldLoadsVector)[load_name];
}

void ElementDataManager::stashLoadsVector()
{
    // Swap the old load with the new load vector
    std::swap(oldLoadsVector, activeLoadsVector);
}