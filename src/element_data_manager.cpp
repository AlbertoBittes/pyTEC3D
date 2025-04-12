#include <../include/element_data_manager.h>

void ElementDataManager::read_loads_arrays(std::string &file_name)
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
        ElementLoad &load_arrays = (*activeLoadsVector)[variable_name];
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

void ElementDataManager::read_update_loads_vectors(std::string file_name)
{
    this->stashLoadsVector();
    // Read the loads data from the file
    this->read_loads_arrays(file_name);
    // Update the old load vector with the new load vector
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
            (*activeLoadsVector)[label][elem->id()] = std::vector<Real>();
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

inline void ElementDataManager::update_time_step(const Substep &time_substep)
{
    this->time_substep = &time_substep;
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