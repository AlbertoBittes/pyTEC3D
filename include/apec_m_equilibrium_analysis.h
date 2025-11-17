#ifndef APECM_EQUILIBRIUM_ANALYSIS_H
#define APECM_EQUILIBRIUM_ANALYSIS_H

#include "transient_analysis.h"
#include "apec_m_system_objects.h"
#include "element_builder.h"

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
    void control_currents();
    void apply_control_currents(std::pair<Real, Real> &correction_currents);
    std::pair<Real, Real> compute_correction_currents();

    Real check_limiter();
    Real check_X_point();

    Real get_psi_magnetic_axis() const { return psi_magnetic_axis; }
    void set_psi_magnetic_axis();

    DenseVector<Real> &get_b_target_axis() { return b_target_axis; };
    void set_b_target_axis();

    void identify_plasma_region();
    void reset_plasma_components();
    void remove_plasma_out_of_X_point();

    inline bool check_distance_magnetic_axis(const Node *node)
    {
        return (node->operator-(magnetic_axis)).norm() < min_dist_magnetic_axis;
    };

private:
    Real psi_magnetic_axis;
    Real psi_min;

    Point magnetic_axis;
    Real min_dist_magnetic_axis = 0.3;

    Point target_axis;
    DenseVector<Real> b_target_axis;

    Real plasma_current_Ip;

    // Real bb;
    Real bradx;
    Real brady;
    Real bverx;
    Real bvery;
    Real cpgain;
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
        {
            elem->set_extra_integer(plasma_region_id, 1);
            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(vacuum_boundary_id))
                    node.set_extra_integer(new_vacuum_boundary_id, 1);
        }

    // for (auto &node : mesh.local_node_ptr_range())
    //     if (node->get_extra_integer(vacuum_boundary_id))
    //         node->set_extra_integer(new_vacuum_boundary_id, 1);
}

void APECMEquilibriumAnalysis::set_b_target_axis()
{
    auto eq_sys = this->get_parent_object();
    auto session = this->context<Session>();
    auto &main_sys = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int bx_var_number = nodal_solution_system.variable_number("B_x");
    unsigned int by_var_number = nodal_solution_system.variable_number("B_y");
    const auto &mesh = eq_sys->get_mesh();
    auto builder = ElementBuilder::make_element_builder(eq_sys);
    auto point_locator = mesh.sub_point_locator();
    // Here you will be able to restrict to the subdomain of stranded conductors for speedup
    const Elem *target_elem = point_locator->operator()(target_axis);
    auto &comm = mesh.comm();
    // target_axis(0) = -std::numeric_limits<Real>::max();
    // target_axis(1) = -std::numeric_limits<Real>::max();
    // target_axis(2) = -std::numeric_limits<Real>::max();
    // if (target_elem)
    // {
    builder->build_solution_vectors_from_emat_files(target_elem->id(), session->get_path("emat_files_path"));
    // }
    b_target_axis.resize(2);
    b_target_axis(0) = nodal_solution_system.point_value(bx_var_number, target_axis, target_elem, nodal_solution_system.solution.get());
    b_target_axis(1) = nodal_solution_system.point_value(by_var_number, target_axis, target_elem, nodal_solution_system.solution.get());
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
            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(new_vacuum_boundary_id))
                    node.set_extra_integer(new_vacuum_boundary_id, 0);
        }
    }
}

void APECMEquilibriumAnalysis::set_psi_magnetic_axis()
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

    auto *node_ptr = mesh.query_node_ptr(max_node_id);
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

    return;
}

void APECMEquilibriumAnalysis::identify_plasma_region()
{
    auto eq_sys = this->get_parent_object();
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_id = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_id = mesh.get_elem_integer_index("vacuum_chamber");

    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(vacuum_chamber_id))
            continue;

        for (auto &node : elem->node_ref_range())
        {
            Real psi_dof = node.dof_number(Globals::SystemType::NODAL_SOLUTION, psi_var_number, 0);
            Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
            if (psi_value < psi_min)
            {
                elem->set_extra_integer(plasma_region_id, 0);
                break;
            }
        }
    }
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

Real APECMEquilibriumAnalysis::check_X_point()
{
    auto eq_sys = this->get_parent_object();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();
    Real psi_X_point = std::numeric_limits<Real>::max();

    dof_id_type vacuum_chamber_region_id = mesh.get_node_integer_index("vacuum_chamber_nodes");

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

void APECMEquilibriumAnalysis::find_separatrix()
{
    this->reset_plasma_components();
    this->set_psi_magnetic_axis();

    Real psi_limiter = this->check_limiter();
    Real psi_X_point = this->check_X_point();
    if (psi_X_point < psi_limiter)
        psi_min = psi_limiter;
    else
        psi_min = psi_X_point;

    this->remove_plasma_out_of_X_point();
    this->identify_plasma_region();
}

std::pair<Real, Real> APECMEquilibriumAnalysis::compute_correction_currents()
{
    Real bb = bradx * bvery - bverx * brady;
    Real dira = -(bverx * b_target_axis(1) - bvery * b_target_axis(0)) / bb;
    // delta vertical current
    Real dive = -(brady * b_target_axis(0) - bradx * b_target_axis(1)) / bb;

    // Integral component
    // xiradi = xiradi + cigain*dira
    // Proportional+Integral component
    Real xirad = -cpgain * dira;

    // Integral component
    // xiveri = xiveri + cigain*dive
    // Proportional+Integral component
    Real xiver = -cpgain * dive;
    return std::make_pair(xirad, xiver);
}

void apply_control_currents(std::pair<Real, Real> &correction_currents)
{
    auto eq_sys = this->get_parent_object();
    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int vertical_coil_current_var_number = element_loads_system.variable_number("Js_z");

    for (auto &elem : eq_sys->get_mesh().active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(eq_sys->get_mesh().get_elem_integer_index("radial_coil")))
            continue;

        for (auto &node : elem->node_ref_range())
        {
            Real iradial_dof = node.dof_number(Globals::SystemType::EXPLICIT_SYSTEM, radial_coil_current_var_number, 0);
            Real iverical_dof = node.dof_number(Globals::SystemType::EXPLICIT_SYSTEM, vertical_coil_current_var_number, 0);

            Real iradial_value = coil_system.current_local_solution->el(iradial_dof);
            Real iverical_value = coil_system.current_local_solution->el(iverical_dof);

            coil_system.current_local_solution->set(iradial_dof, iradial_value + correction_currents.first);
            coil_system.current_local_solution->set(iverical_dof, iverical_value + correction_currents.second);
        }
    }
}
class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(ApecMEquationSystems *eq_systems = nullptr);
};

#endif // APECM_EQUILIBRIUM_ANALYSIS_H
