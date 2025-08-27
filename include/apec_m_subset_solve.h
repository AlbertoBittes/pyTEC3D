

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include <libmesh/system_subset_by_subdomain.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include "apec_m_system_objects.h" // <-- Add this line

#include "fem_service_abstract.h"
#include "time_step.h"

using namespace libMesh;

class ApecMSubset : public SystemSubsetBySubdomain
{
public:
    using SystemSubsetBySubdomain::SystemSubsetBySubdomain; // Inherit constructors

    // void set_tree_gauge()
    // {
    //     // apply tree-cotree gauge
    //     const EquationSystems &equation_systems = this->get_system().get_equation_systems();
    //     const MeshBase &mesh = equation_systems.get_mesh();

    //     if (equation_systems.has_system("SupportSystem"))
    //     {
    //         const System &tree_cotree_gauge_system = equation_systems.get_system<System>("SupportSystem");
    //         if (tree_cotree_gauge_system.has_variable("TREE"))
    //         {
    //             unsigned int az_var_num = this->get_system().variable_number("Az");
    //             unsigned int main_sys_num = this->get_system().number();
    //             unsigned int tree_sys_num = tree_cotree_gauge_system.number();
    //             unsigned int tree_var_num = tree_cotree_gauge_system.variable_number("TREE");
    //             NumericVector<Number> &tree_cotree_gauge = *tree_cotree_gauge_system.current_local_solution;

    //             // Collect DOF ids to remove
    //             std::vector<unsigned int> to_remove;
    //             for (auto node_ptr : mesh.local_node_ptr_range())
    //             {
    //                 if (node_ptr->n_dofs(tree_sys_num, tree_var_num) > 0)
    //                 {
    //                     dof_id_type tree_dof_id = node_ptr->dof_number(tree_sys_num, tree_var_num, 0);
    //                     Real value = tree_cotree_gauge(tree_dof_id);
    //                     if (value > 0.5)
    //                     {
    //                         dof_id_type az_dof_id = node_ptr->dof_number(main_sys_num, az_var_num, 0);
    //                         to_remove.push_back(az_dof_id);
    //                     }
    //                 }
    //             }
    //             // Remove the constrained ids from _dof_ids
    //             auto &dof_ids_vec = this->_dof_ids;
    //             libMesh::out << "Setting tree gauge on " << to_remove.size() << " DOFs" << std::endl;
    //             dof_ids_vec.erase(
    //                 std::remove_if(
    //                     dof_ids_vec.begin(),
    //                     dof_ids_vec.end(),
    //                     [&to_remove](unsigned int id)
    //                     {
    //                         return std::find(to_remove.begin(), to_remove.end(), id) != to_remove.end();
    //                     }),
    //                 dof_ids_vec.end());
    //         }
    //     }
    // }
    // Static factory method to create the subset if needed
    static ApecMSubset *create(ApecMLinearSystem &system)
    {
        // Get mesh and equation systems from system
        EquationSystems &equation_systems = system.get_equation_systems();

        bool has_volt = system.has_variable("VOLT");
        bool transient_effect = system.get_assembler()->context<TimeStep>()->get_transient_effect();

        // Retrieve ElementTypesFactory from system context
        auto element_type_factory = system.context<ElementTypesFactory>();
        if (!element_type_factory)
            return nullptr;

        // Collect subdomains from element_type_factory
        std::set<subdomain_id_type> subdomains;
        for (const auto &pair : *element_type_factory->get_element_types_map())
        {
            subdomains.insert(pair.first);
        }
        SubdomainSelectionByList selection(subdomains);

        // Collect all variable numbers except VOLT
        std::vector<unsigned int> all_var_nums;
        system.get_all_variable_numbers(all_var_nums);
        std::set<unsigned int> all_var_nums_set(all_var_nums.begin(), all_var_nums.end());
        if (has_volt)
            all_var_nums_set.erase(system.variable_number("VOLT"));

        // Check if we need to restrict solve (non-transient and has VOLT variable)
        bool restrict_solve = !transient_effect && has_volt;

        // // Check if tree gauge should be applied
        // bool apply_tree_gauge = false;
        // if (equation_systems.has_system("SupportSystem"))
        // {
        //     const System &tree_cotree_gauge_system = equation_systems.get_system<System>("SupportSystem");
        //     if (tree_cotree_gauge_system.has_variable("TREE"))
        //         apply_tree_gauge = true;
        // }

        // if (restrict_solve || apply_tree_gauge)
        if (restrict_solve)
        {
            auto *subset_ptr = new ApecMSubset(system, selection, &all_var_nums_set);
            // if (apply_tree_gauge)
            //     subset_ptr->set_tree_gauge();
            return subset_ptr;
        }
        return nullptr;
    }
};