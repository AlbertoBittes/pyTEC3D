#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h" // Include the header file for LinearImplicitSystem
#include "../include/node_bc_object.h"

using namespace libMesh;
// BCs are probably broken now

void NodeBCObject::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.
    // Retrieve the system by name
    MeshBase &mesh = system->get_mesh();
    DofMap &dof_map = system->get_dof_map();

    const Node *node_ptr = nullptr;
    for (const auto &bc_pair : *(active_bc_map))
    {
        const auto &bc_values = bc_pair.second;

        // This will not work for edge elements
        // Access the node pointer using the node_id, you will have to replace this with something taht can apply bcs for edge elements
        node_ptr = mesh.node_ptr(bc_pair.first);
        if (node_ptr)
        {
            for (const auto &bc_value : bc_values)
            {
                if (system->has_variable(bc_value.first))
                {
                    dof_id_type dof_id = node_ptr->dof_number(system->number(), system->variable_number(bc_value.first), 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, bc_value.second, true);
                }
                else
                    libmesh_error_msg("Variable " + bc_value.first + " not found in system " + system->name());
            }
        }
    }

    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    if (!this->parent_load_step->get_time_step()->get_transient_effect() && system->has_variable("VOLT"))
    {
        // Apply the constraint for the time derivative of the magnetic vector potential
        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, system->get_mesh(), system->variable_number("VOLT"));
        for (const auto &dof : dof_indices)
        {
            // Apply the constraint
            DofConstraintRow constraint_row;
            dof_map.add_constraint_row(dof, constraint_row, true);
        }
    }
}

void NodeBCObject::read_boundary_conditions(std::istream &file)
{
    libMesh::out << "Reading boundary conditions data..." << std::endl;

    // Initialize the map if it is a nullptr
    if (!active_bc_map)
    {
        active_bc_map = std::make_shared<BoundaryConditionsMap>();
    }

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
        uint8_t label_length;
        Real value;

        // Read node number
        file.read(reinterpret_cast<char *>(&node), sizeof(node));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read node number from the binary stream.");
        }

        // Read label length
        file.read(reinterpret_cast<char *>(&label_length), sizeof(label_length));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read label length from the binary stream.");
        }

        // Read label
        std::vector<char> label(label_length);
        file.read(label.data(), label_length);
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read label from the binary stream.");
        }

        // Read value
        file.read(reinterpret_cast<char *>(&value), sizeof(value));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read value from the binary stream.");
        }

        // Populate the map
        std::string label_str(label.data(), label_length); // Convert char array to string
        (*active_bc_map)[node].emplace(std::move(label_str), value);
    }

    // Do not close the file, as another function will continue reading
}

void NodeBCObject::check_boundary_conditions(LoadStep *prev_load_step)
{
    // Example implementation
    // This function should implement the logic to check the boundary conditions applied to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be checked.

    // Retrieve the system by name
    this->update_boundary_conditions = false;

    if (prev_load_step == nullptr)
        this->update_boundary_conditions = true;
    else
    {
        if (*(prev_load_step->get_node_bc()->get_active_bc_map()) == *active_bc_map)
            this->update_boundary_conditions = false;
        else
            this->update_boundary_conditions = true;
    }
}

void NodeBCObject::set_boundary_conditions()
{
    this->check_boundary_conditions();

    if (this->update_boundary_conditions)
    {
        system->reinit_constraints();
    }
}