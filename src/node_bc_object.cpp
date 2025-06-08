#include "node_bc_object.h"
#include "globals.h"

using namespace libMesh;

void NodeBCObject::read_data(std::ifstream &file)
{
    libMesh::out << "Reading boundary conditions data..." << std::endl;

    // Read the number of boundary conditions
    unsigned int num_bcs;
    file.read(reinterpret_cast<char *>(&num_bcs), sizeof(num_bcs));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of boundary conditions from the binary stream.");
    }

    // Read each boundary condition
    for (unsigned int i = 0; i < num_bcs; ++i)
    {
        dof_id_type node;
        unsigned int label_number;
        unsigned int component_number;
        libMesh::Real value;

        // Read node number
        file.read(reinterpret_cast<char *>(&node), sizeof(node));

        // Read label number
        file.read(reinterpret_cast<char *>(&label_number), sizeof(label_number));

        // Read label
        file.read(reinterpret_cast<char *>(&component_number), sizeof(component_number));

        // Read value
        file.read(reinterpret_cast<char *>(&value), sizeof(value));

        Globals::APIVariableName api_var_enum = static_cast<Globals::APIVariableName>(label_number);
        const auto &global_vars = Globals::API_TO_GLOBAL_VARIABLE_MAP.at(api_var_enum);
        const std::string &label_str = global_vars[component_number];
        Globals::GlobalVariableName label_enum =
            Globals::string_to_enum(label_str, Globals::STRING_TO_ENUM, Globals::GlobalVariableName::UNKNOWN);

        // Populate the map
        (*data_map)[node][label_enum] = value;
    }

    // Do not close the file, as another function will continue reading
}
