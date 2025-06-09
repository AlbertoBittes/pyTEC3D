#include "node_bc_object.h"
#include "globals.h"

using namespace libMesh;

void NodeBCObject::read_data(std::ifstream &file)
{
    libMesh::out << "Reading boundary conditions data..." << std::endl;

    // Use dof_id_type for binary I/O
    dof_id_type num_bcs = 0;
    file.read(reinterpret_cast<char *>(&num_bcs), sizeof(dof_id_type));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of boundary conditions from the binary stream.");
    }

    // Read each boundary condition
    for (dof_id_type i = 0; i < num_bcs; ++i)
    {
        dof_id_type node = 0;
        dof_id_type label_number = 0;
        dof_id_type component_number = 0;
        libMesh::Real value = 0.0;

        // Read node number
        file.read(reinterpret_cast<char *>(&node), sizeof(dof_id_type));

        // Read label number
        file.read(reinterpret_cast<char *>(&label_number), sizeof(dof_id_type));

        // Read component number
        file.read(reinterpret_cast<char *>(&component_number), sizeof(dof_id_type));

        // Read value
        file.read(reinterpret_cast<char *>(&value), sizeof(libMesh::Real));

        Globals::APIVariableName api_var_enum = static_cast<Globals::APIVariableName>(static_cast<unsigned int>(label_number));
        const auto &global_vars = Globals::API_TO_GLOBAL_VARIABLE_MAP.at(api_var_enum);
        const std::string &label_str = global_vars[static_cast<unsigned int>(component_number)];
        Globals::GlobalVariableName label_enum =
            Globals::string_to_enum(label_str, Globals::STRING_TO_ENUM, Globals::GlobalVariableName::UNKNOWN);

        // Populate the map
        (*data_map)[node][label_enum] = value;
    }

    // Do not close the file, as another function will continue reading
}
