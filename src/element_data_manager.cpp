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

    // Read the length of the load name
    // size_t load_name_length;
    // file.read(reinterpret_cast<char *>(&load_name_length), sizeof(load_name_length));

    // // Read the load name
    // std::string load_name(load_name_length, '\0');
    // file.read(&load_name[0], load_name_length);

    // // Check if the load name exists in the activeLoadsVector
    // if (!activeLoadsVector || activeLoadsVector->find(load_name) == activeLoadsVector->end())
    // {
    //     libmesh_error_msg("Error: Load name " + load_name + " not found in activeLoadsVector.");
    // }

    // // Read the number of components
    // unsigned int num_components;
    // file.read(reinterpret_cast<char *>(&num_components), sizeof(num_components));

    // Get the ElementLoad map for the given load name
    // ElementLoad &element_load = (*activeLoadsVector)[load_name];

    // Read the element data and modify the existing pairs in the map
    while (file.peek() != EOF)
    {
        dof_id_type element_number;
        unsigned int component_number;
        Real load_value;

        // Read element number
        file.read(reinterpret_cast<char *>(&element_number), sizeof(element_number));

        // Read component number
        file.read(reinterpret_cast<char *>(&component_number), sizeof(component_number));

        // Read load value
        file.read(reinterpret_cast<char *>(&load_value), sizeof(load_value));

        // Modify the existing pair in the map
        auto it = load_array->find(element_number);
        if (it != load_array->end())
        {
            it->second[component_number] = load_value; // Update the component number
        } // Update the load value
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

inline void ElementDataManager::fill_mu_matrix(dof_id_type elem_id)
{
    // get ansys_material_type from the materials map
    const MaterialsProperties &mat_properties = materials->at(elem_id);

    mu_matrix[0] = 4 * libMesh::pi * 1e-7 * mat_properties.murx;
    mu_matrix[1] = 4 * libMesh::pi * 1e-7 * mat_properties.mury;
    mu_matrix[2] = 4 * libMesh::pi * 1e-7 * mat_properties.murz;
}

inline void ElementDataManager::fill_rs_matrix(dof_id_type elem_id)
{
    const MaterialsProperties &mat_properties = materials->at(elem_id);

    rs_matrix[0] = mat_properties.rsvx;
    rs_matrix[1] = mat_properties.rsvy;
    rs_matrix[2] = mat_properties.rsvz;
}

inline void ElementDataManager::compute_material_properties(dof_id_type elem_id, bool conductive)
{
    // Call fill_mu_matrix and fill_rs_matrix to compute the material properties
    this->fill_mu_matrix(elem_id);
    if (conductive)
        this->fill_rs_matrix(elem_id);
}

inline Real ElementDataManager::compute_js(dof_id_type elem_id, unsigned int comp)
{
    // Get the value from the "Js" array using the element ID
    if (time_substep->get_parent().get_ramped_load())
        return this->ramp_current(elem_id, comp);
    else
        return (*load_array)[elem_id][comp];
}

inline Real ElementDataManager::ramp_current(dof_id_type elem_id, unsigned int comp)
{
    Real y1;
    Real y2;
    if (oldLoadsVector != nullptr)
    {
        y1 = (*old_load_array)[elem_id][comp];
        y2 = (*load_array)[elem_id][comp];
    }
    else
    {
        y1 = 0;
        y2 = (*load_array)[elem_id][comp];
    }
    return y1 + ((y2 - y1) / (time_substep->get_parent().get_time() - time_substep->get_parent().get_time_beginning())) * (time_substep->get_time() - time_substep->get_parent().get_time_beginning());
}

inline void ElementDataManager::set_load_quantity(std::string load_name)
{
    // Get the load quantity from the cell data
    load_array = &(*activeLoadsVector)[load_name];
    load_array = &(*oldLoadsVector)[load_name];
}

inline void ElementDataManager::stashLoadsVector()
{
    // Swap the old load with the new load vector
    std::swap(oldLoadsVector, activeLoadsVector);
}