#ifndef APECM_EQUILIBRIUM_ANALYSIS_H
#define APECM_EQUILIBRIUM_ANALYSIS_H

#include "transient_analysis.h"
#include "apec_m_system_objects.h"

#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <limits>

class APECMEquilibriumAnalysis : public TransientAnalysis
{
public:
    explicit APECMEquilibriumAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    void compute_jplasma_gradshafranov();
    void find_separatrix();
    void check_poloidal_flux_convergence();
    void get_b_central_axis();
    void apply_control_currents();

    Real check_limiter();
    Real check_diverted();
    Real get_psi_magnetic_axis();
    void identify_plasma_region();
    void reset_plasma_components();
    void remove_plasma_out_of_X_point();

    inline bool check_distance_magnetic_axis(const Node *node)
    {
        return (node->operator-(magnetic_axis)).norm() < min_dist_magnetic_axis;
    };

private:
    Real psi_min;
    Real psi_magnetic_axis;
    std::vector<Real> b_magnetic_axis;
    Point magnetic_axis;
    Real min_dist_magnetic_axis = 0.3;
    Real plasma_current_Ip;
};

void APECMEquilibriumAnalysis::reset_plasma_components()
{
    auto eq_sys = this->get_parent_object();
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_id = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_id = mesh.get_elem_integer_index("vacuum_chamber");
    dof_id_type vacuum_boundary_id = mesh.get_node_integer_index("vacuum_chamber_boundary");
    dof_id_type new_vacuum_boundary_id = mesh.get_node_integer_index("new_vacuum_chamber_boundary");

    for (auto &elem : mesh.active_local_element_ptr_range())
        if (elem->get_extra_integer(vacuum_chamber_id))
            elem->set_extra_integer(plasma_region_id, 1);

    for (auto &node : mesh.local_node_ptr_range())
        if (node->get_extra_integer(vacuum_boundary_id))
            node->set_extra_integer(new_vacuum_boundary_id, 1);
}

void APECMEquilibriumAnalysis::remove_plasma_out_of_X_point()
{
    auto eq_sys = this->get_parent_object();
    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    // unsigned int psi_var_number = element_solution_system.variable_number("PSI");
    unsigned int bx_var_number = element_solution_system.variable_number("Be_x");
    unsigned int by_var_number = element_solution_system.variable_number("Be_y");
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_id = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_id = mesh.get_elem_integer_index("vacuum_chamber");
    // dof_id_type vacuum_boundary_id = mesh.get_node_integer_index("vacuum_chamber_boundary");
    dof_id_type new_vacuum_boundary_id = mesh.get_node_integer_index("new_vacuum_chamber_boundary");

    int current_sign = (plasma_current_Ip >= 0) ? 1 : -1;

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(vacuum_chamber_id))
            continue;

        Real bx_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, bx_var_number, 0);
        Real by_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, by_var_number, 0);
        Real bx_value = element_solution_system.current_local_solution->el(bx_dof);
        Real by_value = element_solution_system.current_local_solution->el(by_dof);

        if (current_sign * (by_value * (elem->centroid()(0) - magnetic_axis(0)) - bx_value * (elem->centroid()(1) - magnetic_axis(1))) < 0)
        {
            elem->set_extra_integer(plasma_region_id, 0);
            for (unsigned int n = 0; n < elem->n_nodes(); n++)
            {
                Node *node_ptr = elem->node_ptr(n);
                if (node_ptr->get_extra_integer(new_vacuum_boundary_id))
                    node_ptr->set_extra_integer(new_vacuum_boundary_id, 0);
            }
        }
    }
}

Real APECMEquilibriumAnalysis::get_psi_magnetic_axis()
{
    auto eq_sys = this->get_parent_object();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();

    dof_id_type vacuum_chamber_region_id = mesh.get_node_integer_index("vacuum_chamber");

    std::pair<dof_id_type, Real> magnetic_axis_data;

    auto &comm = mesh.comm();

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_chamber_region_id);
        if (!is_vacuum_chamber)
            continue;

        if (!check_distance_magnetic_axis(node))
            continue;
        // Process the element as needed
        Real psi_dof = node->dof_number(Globals::SystemType::NODAL_SOLUTION, psi_var_number, 0);
        Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
        if (psi_value > magnetic_axis_data.second)
            magnetic_axis_data = std::make_pair(node->id(), psi_value);
    }
    std::set<std::pair<dof_id_type, Real>> magnetic_axis_data_all;
    magnetic_axis_data_all.insert(magnetic_axis_data);
    comm.set_union(magnetic_axis_data_all);
    // Find the global maximum
    dof_id_type max_node_id = 0;
    Real max_psi_value = -std::numeric_limits<Real>::max();
    for (const auto &data : magnetic_axis_data_all)
    {
        if (data.second > max_psi_value)
        {
            max_psi_value = data.second;
            max_node_id = data.first;
        }
    }

    psi_magnetic_axis = max_psi_value;

    auto *node_ptr = mesh.node_ptr(max_node_id);
    Real axis_x = -std::numeric_limits<Real>::max();
    Real axis_y = -std::numeric_limits<Real>::max();
    Real axis_z = -std::numeric_limits<Real>::max();
    if (node_ptr)
    {
        axis_x = (*node_ptr)(0);
        axis_y = (*node_ptr)(1);
        axis_z = (*node_ptr)(2);
    }
    comm.max(axis_x);
    comm.max(axis_y);
    comm.max(axis_z);
    magnetic_axis(0) = axis_x;
    magnetic_axis(1) = axis_y;
    magnetic_axis(2) = axis_z;

    return psi_magnetic_axis;
}

Real APECMEquilibriumAnalysis::check_limiter()
{
    auto eq_sys = this->get_parent_object();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    unsigned int bx_var_number = nodal_solution_system.variable_number("B_x");
    unsigned int by_var_number = nodal_solution_system.variable_number("B_y");
    const auto &mesh = eq_sys->get_mesh();
    Real psi_limiter = -std::numeric_limits<Real>::max();

    dof_id_type vacuum_boundary_region_id = mesh.get_node_integer_index("vacuum_chamber_boundary");

    int current_sign = (plasma_current_Ip >= 0) ? 1 : -1;

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_boundary_region_id);
        if (!is_vacuum_chamber)
            continue;

        Real bx_dof = node->dof_number(Globals::SystemType::NODAL_SOLUTION, bx_var_number, 0);
        Real by_dof = node->dof_number(Globals::SystemType::NODAL_SOLUTION, by_var_number, 0);
        Real bx_value = nodal_solution_system.current_local_solution->el(bx_dof);
        Real by_value = nodal_solution_system.current_local_solution->el(by_dof);

        if (current_sign * (by_value * ((*node)(0) - magnetic_axis(0)) - bx_value * ((*node)(1) - magnetic_axis(1))) < 0)
            continue;
        // Process the element as needed
        dof_id_type psi_dof = node->dof_number(Globals::SystemType::NODAL_SOLUTION, psi_var_number, 0);
        Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
        if (psi_value > psi_limiter)
            psi_limiter = psi_value;
    }
    mesh.comm().max(psi_limiter);
    return psi_limiter;
}

Real APECMEquilibriumAnalysis::check_diverted()
{
    auto eq_sys = this->get_parent_object();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();
    Real psi_X_point = std::numeric_limits<Real>::max();

    dof_id_type vacuum_chamber_region_id = mesh.get_node_integer_index("vacuum_chamber");

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_chamber_region_id);
        if (!is_vacuum_chamber)
            continue;

        if (!check_distance_magnetic_axis(node))
            continue;
        // Process the element as needed
        Real psi_dof = node->dof_number(Globals::SystemType::NODAL_SOLUTION, psi_var_number, 0);
        Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
        if (psi_value < psi_X_point)
            psi_X_point = psi_value;
    }
    mesh.comm().min(psi_X_point);
    return psi_X_point;
}

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(ApecMEquationSystems *eq_systems = nullptr);
};

#endif // APECM_EQUILIBRIUM_ANALYSIS_H
