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

    std::string variable_label;
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
                if (bc_value.second.first == 0)
                {
                    if (system->has_variable(bc_value.first))
                    {
                        variable_label = bc_value.first;
                    }
                    else
                    {
                        variable_label = bc_value.first + "_x";
                    }
                }
                else if (bc_value.second.first == 1)
                {
                    variable_label = bc_value.first + "_y";
                }
                else if (bc_value.second.first == 2)
                {
                    variable_label = bc_value.first + "_z";
                }
                if (system->has_variable(variable_label))
                {
                    dof_id_type dof_id = node_ptr->dof_number(system->number(), system->variable_number(variable_label), 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, bc_value.second.second, true);
                }
                else
                    libmesh_error_msg("Variable " + variable_label + " not found");
            }
        }
        else
        {
            libmesh_error_msg("Node " + std::to_string(bc_pair.first) + " not found in mesh.");
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

    std::vector<std::string> labels = {"A", "VOLT", "Az"};
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
        unsigned int label_number;
        unsigned int component_number;
        Real value;

        // Read node number
        file.read(reinterpret_cast<char *>(&node), sizeof(node));

        // Read label number
        file.read(reinterpret_cast<char *>(&label_number), sizeof(label_number));

        // Read label
        file.read(reinterpret_cast<char *>(&component_number), sizeof(component_number));

        // Read value
        file.read(reinterpret_cast<char *>(&value), sizeof(value));

        // Populate the map
        (*active_bc_map)[node][labels[label_number]] = std::make_pair(component_number, value);
    }

    // Do not close the file, as another function will continue reading
}

void NodeBCObject::check_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step)
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

void NodeBCObject::set_boundary_conditions(std::shared_ptr<LoadStep> prev_load_step)
{
    this->check_boundary_conditions(prev_load_step);

    if (this->update_boundary_conditions)
    {
        system->reinit_constraints();
    }
}