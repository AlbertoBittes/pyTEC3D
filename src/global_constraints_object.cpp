#include "fe_parent_object.h"
#include "global_constraints_object.h"
#include "apec_m_system_objects.h"
#include "globals.h"
#include "fem_data_repository.h"
#include "fe_data_object.h"
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>

void GlobalConstraintsObject::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.
    // Retrieve the system by name
    // at system init
    DofMap &dof_map = this->parent_object->get_dof_map();
    MeshBase &mesh = this->parent_object->get_mesh();
    EquationSystems &equation_systems = this->parent_object->get_equation_systems();
    if (!this->ctx_view)
    {
        return;
    }
    libMesh::out << "Applying constraints to the system." << std::endl;

    NumericVector<Number> &node_dof_constraints = *this->get_parent_object()->request_vector("ramped_node_dof_constraints");
    // std::string string_variable_label;
    for (numeric_index_type k = node_dof_constraints.first_local_index(); k < node_dof_constraints.last_local_index(); ++k)
    {
        Real value = node_dof_constraints(k);
        if (!std::isnan(value))
        {
            DofConstraintRow constraint_row;
            dof_map.add_constraint_row(k, constraint_row, value, true);
        }
    }

    // apply tree-cotree gauge
    if (equation_systems.has_system("SupportSystem"))
    {
        System &tree_cotree_gauge_system = equation_systems.get_system<System>("SupportSystem");
        if (tree_cotree_gauge_system.has_variable("TREE"))
        {
            unsigned int az_var_num = this->get_parent_object()->variable_number("Az");
            unsigned int main_sys_num = this->get_parent_object()->number();
            unsigned int tree_sys_num = tree_cotree_gauge_system.number();
            unsigned int tree_var_num = tree_cotree_gauge_system.variable_number("TREE");
            NumericVector<Number> &tree_cotree_gauge = *tree_cotree_gauge_system.current_local_solution;
            for (auto node_ptr : mesh.local_node_ptr_range())
            {
                dof_id_type tree_dof_id = node_ptr->dof_number(tree_sys_num, tree_var_num, 0);
                Real value = tree_cotree_gauge(tree_dof_id);
                if (value > 0.5)
                {
                    dof_id_type az_dof_id = node_ptr->dof_number(main_sys_num, az_var_num, 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(az_dof_id, constraint_row, true);
                }
            }
        }
    }
}

void GlobalConstraintsObject::check_update_boundary_conditions()
{
    DofMap &dof_map = this->parent_object->get_dof_map();

    // dof_map.swap_dof_constraints();

    this->update_boundary_conditions = false;

    // if (!this->ctx_view)
    // {
    //     libmesh_error_msg("Context view is not set in GlobalConstraintsObject.");
    // }
    // if (!old_ctx_view)
    // {
    //     this->update_boundary_conditions = true;
    //     libMesh::out << "Boundary conditions update required" << std::endl;
    //     return;
    // }

    NumericVector<Number> &node_dof_constraints = *this->get_parent_object()->request_vector("node_dof_constraints");
    NumericVector<Number> &old_node_dof_constraints = *this->get_parent_object()->request_vector("old_node_dof_constraints");

    for (numeric_index_type k = node_dof_constraints.first_local_index(); k < node_dof_constraints.last_local_index(); ++k)
    {
        Real new_val = node_dof_constraints(k);
        Real old_val = old_node_dof_constraints(k);
        // If one is NaN and the other is not,
        if ((std::isnan(new_val) != std::isnan(old_val))) //|| (!std::isnan(new_val) && !std::isnan(old_val) && new_val != old_val))
        {
            this->update_boundary_conditions = true;
            libMesh::out << "New DOF constraints: Matrix update required" << std::endl;
            break;
        }
    }
    return;
}