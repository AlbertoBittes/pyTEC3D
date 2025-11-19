#ifndef APECM_EQUILIBRIUM_ANALYSIS_H
#define APECM_EQUILIBRIUM_ANALYSIS_H

#include "transient_analysis.h"
#include "apec_m_system_objects.h"
#include "element_builder.h"

#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <limits>
#include <utility>
#include <vector>

class APECMEquilibriumAnalysis : public TransientAnalysis
{
public:
    explicit APECMEquilibriumAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    void apply_jplasma_gradshafranov();
    void constant_plasma_current_density();
    void find_separatrix();
    bool check_poloidal_flux_convergence();
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

    bool picard_loop(Substep &substep);

    inline bool check_distance_magnetic_axis(const Node *node)
    {
        return (node->operator-(magnetic_axis)).norm() < min_dist_magnetic_axis;
    };

protected:
    bool initialize_model_context(const std::shared_ptr<FEMDataRepository> &model_repo) override;
    bool execute_solve() override;

    virtual void configure_from_settings();
    void set_plasma_state_from_profile(std::size_t index);
    bool solve_equilibrium_step(TimeStep &time_step);
    bool run_picard_iterations(TimeStep &time_step);
    void prepare_time_step(TimeStep &time_step, bool enforce_static);
    void reset_iteration_state();
    std::vector<Real> plasma_current_profile;
    std::vector<Real> beta_p_profile;
    std::size_t active_profile_index = 0;
    bool settings_loaded = false;

private:
    Real psi_magnetic_axis = 0.0;
    Real psi_min = 0.0;

    Point magnetic_axis = Point();
    Real min_dist_magnetic_axis = 0.3;

    Point target_axis = Point();
    DenseVector<Real> b_target_axis;

    Real plasma_current_Ip = 0.0;
    Real alpha = 1.0;
    Real beta = 0.0;

    // Real bb;
    Real bradx = 0.0;
    Real brady = 0.0;
    Real bverx = 0.0;
    Real bvery = 0.0;
    Real cpgain = 0.0;

    Real R0 = 0.0;
    bool first_absolute_iteration = true;
    unsigned char max_picard_iterations = 40;
    Real tol = 1e-3;
    std::vector<dof_id_type> vacuum_chamber_psi_dofs;
    std::vector<Number> prev_it_vacuum_chamber_psi_vector;
    std::vector<Number> curr_it_vacuum_chamber_psi_vector;
};

void APECMEquilibriumAnalysis::reset_plasma_components()
{
    auto eq_sys = this->get_parent_object();
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("vacuum_chamber");
    dof_id_type vacuum_boundary_component = mesh.get_node_integer_index("vacuum_chamber_boundary");
    dof_id_type new_vacuum_boundary_component = mesh.get_node_integer_index("new_vacuum_chamber_boundary");

    for (auto &elem : mesh.active_local_element_ptr_range())
        if (elem->get_extra_integer(vacuum_chamber_component))
        {
            elem->set_extra_integer(plasma_region_component, 1);
            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(vacuum_boundary_component))
                    node.set_extra_integer(new_vacuum_boundary_component, 1);
        }

    // for (auto &node : mesh.local_node_ptr_range())
    //     if (node->get_extra_integer(vacuum_boundary_component))
    //         node->set_extra_integer(new_vacuum_boundary_component, 1);
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

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("vacuum_chamber");
    // dof_id_type vacuum_boundary_component = mesh.get_node_integer_index("vacuum_chamber_boundary");
    dof_id_type new_vacuum_boundary_component = mesh.get_node_integer_index("new_vacuum_chamber_boundary");

    int current_sign = (plasma_current_Ip >= 0) ? 1 : -1;

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(vacuum_chamber_component))
            continue;

        Real bx_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, bx_var_number, 0);
        Real by_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, by_var_number, 0);
        Real bx_value = element_solution_system.current_local_solution->el(bx_dof);
        Real by_value = element_solution_system.current_local_solution->el(by_dof);

        if (current_sign * (by_value * (elem->centroid()(0) - magnetic_axis(0)) - bx_value * (elem->centroid()(1) - magnetic_axis(1))) < 0)
        {
            elem->set_extra_integer(plasma_region_component, 0);
            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(new_vacuum_boundary_component))
                    node.set_extra_integer(new_vacuum_boundary_component, 0);
        }
    }
}

void APECMEquilibriumAnalysis::set_psi_magnetic_axis()
{
    auto eq_sys = this->get_parent_object();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();

    dof_id_type vacuum_chamber_region_component = mesh.get_node_integer_index("vacuum_chamber");

    std::pair<dof_id_type, Real> magnetic_axis_data;

    auto &comm = mesh.comm();

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_chamber_region_component);
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

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("plasma");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("vacuum_chamber");

    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(vacuum_chamber_component))
            continue;

        for (auto &node : elem->node_ref_range())
        {
            Real psi_dof = node.dof_number(Globals::SystemType::NODAL_SOLUTION, psi_var_number, 0);
            Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
            if (psi_value < psi_min)
            {
                elem->set_extra_integer(plasma_region_component, 0);
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

    dof_id_type vacuum_boundary_region_component = mesh.get_node_integer_index("vacuum_chamber_boundary");

    int current_sign = (plasma_current_Ip >= 0) ? 1 : -1;

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_boundary_region_component);
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

    dof_id_type vacuum_chamber_region_component = mesh.get_node_integer_index("vacuum_chamber_nodes");
    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_chamber_region_component);
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

void APECMEquilibriumAnalysis::apply_control_currents(std::pair<Real, Real> &correction_currents)
{
    auto eq_sys = this->get_parent_object();
    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int vertical_coil_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &coils_current_vector = element_loads_system.get_vector("ramped_loads");

    const auto &mesh = eq_sys->get_mesh();
    dof_id_type vertical_control_coils_component = mesh.get_elem_integer_index("vertical_control_coils");
    dof_id_type radial_control_coils_component = mesh.get_elem_integer_index("radial_control_coils");
    dof_id_type real_constant_mesh_id = mesh.get_elem_integer_index("ansys_real_constant");

    auto real_constants_map = this->context<RealConstantsObject>();

    Real correction_current;

    for (auto &elem : eq_sys->get_mesh().active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(vertical_control_coils_component))
            correction_current = correction_currents.second;
        else if (elem->get_extra_integer(radial_control_coils_component))
            correction_current = correction_currents.first;
        else
            continue;

        dof_id_type real_constant_elem_id = elem->get_extra_integer(real_constant_mesh_id);
        Real coil_area = real_constants_map->get_constant(real_constant_elem_id)[0];
        Real correction_current_density = correction_current / coil_area;

        char correction_current_sign = real_constants_map->get_constant(real_constant_elem_id)[3];
        dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, vertical_coil_current_var_number, 0);
        Real existing_current_density = element_loads_system.current_local_solution->el(elem_current_density_dof);
        Real new_current_density = existing_current_density + correction_current_sign * correction_current_density;
        coils_current_vector.set(elem_current_density_dof, new_current_density);
    }

    first_absolute_iteration = false;
}

void APECMEquilibriumAnalysis::apply_jplasma_gradshafranov()
{
    auto eq_sys = this->get_parent_object();

    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");

    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int plasma_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &plasma_current_vector = element_loads_system.get_vector("ramped_loads");

    const auto &mesh = eq_sys->get_mesh();
    auto &comm = mesh.comm();
    dof_id_type plasma_component = mesh.get_elem_integer_index("plasma");

    Real correction_current;
    Real total_plasma_current = 0.0;

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(plasma_component))
            continue;

        dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, plasma_current_var_number, 0);
        dof_id_type psi_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, psi_var_number, 0);

        Real psi_elem = element_solution_system.current_local_solution->el(psi_dof);
        Real elem_centroid_rad = elem->centroid()(0);

        Real psi_segn = (psi_elem - psi_magnetic_axis) / (psi_min - psi_magnetic_axis);
        Real grad_shafranov_current_density = (1 - std::pow(psi_segn, alpha)) *
                                              (beta * elem_centroid_rad / R0 + (1 - beta) * R0 / elem_centroid_rad);
        total_plasma_current += grad_shafranov_current_density * elem->volume();
        plasma_current_vector.set(elem_current_density_dof, grad_shafranov_current_density);
    }

    // Take total current across all processors
    comm.sum(total_plasma_current);

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(plasma_component))
            continue;

        dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, plasma_current_var_number, 0);
        Real existing_current_density = element_loads_system.current_local_solution->el(elem_current_density_dof);
        Real adjusted_current_density = existing_current_density * (plasma_current_Ip / total_plasma_current);
        plasma_current_vector.set(elem_current_density_dof, adjusted_current_density);
    }
}

bool APECMEquilibriumAnalysis::check_poloidal_flux_convergence()
{
    auto eq_sys = this->get_parent_object();

    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");

    const auto &mesh = eq_sys->get_mesh();
    auto &comm = mesh.comm();
    dof_id_type vacuum_chamber_id = mesh.get_elem_integer_index("vacuum_chamber");

    if (first_absolute_iteration)
    {
        for (auto &elem : mesh.active_local_element_ptr_range())
        {
            if (elem->get_extra_integer(vacuum_chamber_id))
            {
                dof_id_type psi_dof = elem->dof_number(Globals::SystemType::ELEMENT_SOLUTION, psi_var_number, 0);
                vacuum_chamber_psi_dofs.push_back(psi_dof);
            }
        }
        prev_it_vacuum_chamber_psi_vector.resize(vacuum_chamber_psi_dofs.size(), 0.0);
        return false;
    }

    bool converged = true;

    element_solution_system.current_local_solution->get(vacuum_chamber_psi_dofs, curr_it_vacuum_chamber_psi_vector);

    for (std::size_t i = 0; i < curr_it_vacuum_chamber_psi_vector.size(); ++i)
    {
        Real rel_diff = std::abs(curr_it_vacuum_chamber_psi_vector[i] - prev_it_vacuum_chamber_psi_vector[i]) / std::max(std::abs(curr_it_vacuum_chamber_psi_vector[i]), std::abs(prev_it_vacuum_chamber_psi_vector[i]));
        if (rel_diff > tol)
        {
            converged = false;
            break;
        }
    }

    // Check if converged in all processors
    comm.min(converged);

    std::swap(prev_it_vacuum_chamber_psi_vector, curr_it_vacuum_chamber_psi_vector);
    return converged;
}

bool APECMEquilibriumAnalysis::picard_loop(Substep &substep)
{
    if (!first_absolute_iteration)
        find_separatrix();

    if (first_absolute_iteration)
        constant_plasma_current_density();
    else
        apply_jplasma_gradshafranov();

    bool converged = check_poloidal_flux_convergence();

    if (converged)
        return true;

    auto correction_currents = compute_correction_currents();
    apply_control_currents(correction_currents);

    solve_equation_systems(substep);
    return false;
}

void APECMEquilibriumAnalysis::constant_plasma_current_density()
{
    auto eq_sys = this->get_parent_object();

    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int plasma_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &plasma_current_vector = element_loads_system.get_vector("ramped_loads");

    const auto &mesh = eq_sys->get_mesh();
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("vacuum_chamber");

    auto real_constants_map = this->context<RealConstantsObject>();
    dof_id_type real_constant_mesh_id = mesh.get_elem_integer_index("ansys_real_constant");

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!elem->get_extra_integer(vacuum_chamber_component))
            continue;

        dof_id_type real_constant_elem_id = elem->get_extra_integer(real_constant_mesh_id);
        Real plasma_area = real_constants_map->get_constant(real_constant_elem_id)[0];
        Real plasma_current_density = plasma_current_Ip / plasma_area;
        dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, plasma_current_var_number, 0);
        plasma_current_vector.set(elem_current_density_dof, plasma_current_density);
    }
}

#endif // APECM_EQUILIBRIUM_ANALYSIS_H
