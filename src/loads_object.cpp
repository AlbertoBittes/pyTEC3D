#include "loads_object.h"
#include "globals.h"
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>

void LoadsObject::read_data(std::ifstream &file)
{
    libMesh::out << "Reading loads data from an open file stream..." << std::endl;

    // Use dof_id_type (uint32_t) for binary I/O
    dof_id_type num_loads = 0;
    file.read(reinterpret_cast<char *>(&num_loads), sizeof(dof_id_type));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of loads from the binary stream.");
    }

    std::string load_label_api;
    std::string load_label_global;
    unsigned int sys_number = data_map->number();
    const MeshBase &mesh = data_map->get_mesh();
    dof_id_type dof_id;

    for (dof_id_type i = 0; i < num_loads; ++i)
    {
        dof_id_type element_number = 0;
        dof_id_type load_label_number = 0;
        dof_id_type component_number = 0;
        Real load_value = 0.0;

        file.read(reinterpret_cast<char *>(&element_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&load_label_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&component_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&load_value), sizeof(Real));

        // Map the load label number to the corresponding API enum
        Globals::APIVariableName api_var_enum = static_cast<Globals::APIVariableName>(static_cast<unsigned int>(load_label_number));

        // Use the map to get the global variable string
        const auto &global_vars = Globals::API_TO_GLOBAL_VARIABLE_MAP.at(api_var_enum);
        load_label_global = global_vars[static_cast<unsigned int>(component_number)];

        const Elem *elem_ptr = mesh.query_elem_ptr(static_cast<unsigned int>(element_number));

        if (!elem_ptr)
        {
            continue;
        }
        else
        {
            dof_id = elem_ptr->dof_number(sys_number,
                                          static_cast<unsigned int>(load_label_number),
                                          static_cast<unsigned int>(component_number));
            data_map->solution->set(dof_id, load_value);
        }
    }
}

// void LoadsObject::zero_loads_array()
// {
//     if (!data_map)
//     {
//         return;
//     }
//     // Iterate over each load in the LoadArraysMap
//     for (auto &load_array_pair : *data_map)
//     {
//         ElementLoad &element_load = load_array_pair.second; // Access the ElementLoad map

//         // Iterate over each element in the ElementLoad map
//         for (auto &element_pair : element_load)
//         {
//             std::vector<Real> &values = element_pair.second; // Access the vector of Real values

//             // Zero out the vector using std::fill for efficiency
//             if (!values.empty())
//             {
//                 std::fill(values.begin(), values.end(), 0.0);
//             }
//         }
//     }
// }

// void LoadsObject::build_load_data_vectors(const MeshBase &mesh)
// {

//     // Get the load quantity from the cell data
//     for (const auto &label : Globals::API_LOADS_LABELS)
//     {
//         // If not, create a new load vector for this label
//         (*data_map)[label] = ElementLoad();
//         for (const auto &elem : mesh.active_local_element_ptr_range())
//         {
//             // If not, create a new load for this element ID
//             (*data_map)[label][elem->id()] = std::vector<Real>(3, 0.0);
//         }
//     }

//     if (!oldLoadsVector)
//     {
//         oldLoadsVector = std::make_unique<LoadArraysMap>();
//     }
//     *oldLoadsVector = *data_map;
// }
