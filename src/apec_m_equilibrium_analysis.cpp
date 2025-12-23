#include "apec_m_equilibrium_analysis.h"

#include "apec_m_settings.h"
#include "session.h"

#include <algorithm>
#include <libmesh/libmesh.h>

APECMEquilibriumAnalysis::APECMEquilibriumAnalysis(APEC_M *program,
                                                   std::shared_ptr<AnalysisBase::ContextView> context_view,
                                                   AnalysisBase::Type type)
    : TransientAnalysis(type, program, context_view)
{
}

bool APECMEquilibriumAnalysis::run_picard_iterations(Substep &substep)
{
    auto settings = this->context<ApecMSettingsObject>();
    const auto iterations_limit = std::max<unsigned int>(1u, settings->max_picard_iterations());
    bool converged = false;

    for (unsigned int iteration = 0; iteration < iterations_limit; ++iteration)
    {
        if (picard_loop(substep))
        {
            libMesh::out << "Equilibrium analysis converged in " << iteration + 1 << " Picard iterations." << std::endl;
            converged = true;
            break;
        }
    }

    if (!converged)
    {
        libMesh::out << "Equilibrium analysis did not converge within " << iterations_limit
                     << " Picard iterations." << std::endl;
        return false;
    }

    return converged;
}

void APECMEquilibriumAnalysis::reset_plasma_components()
{
    libMesh::out << "Resetting plasma region and boundary components..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("PLASMA");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");
    dof_id_type vacuum_boundary_component = mesh.get_node_integer_index("VACUUM_CHAMBER_BOUNDARY");
    dof_id_type new_vacuum_boundary_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");

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
    libMesh::out << "Finding B field at axis point..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto session = this->context<Session>();
    auto &main_sys = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int bx_var_number = nodal_solution_system.variable_number("B_x");
    unsigned int by_var_number = nodal_solution_system.variable_number("B_y");
    const auto &mesh = eq_sys->get_mesh();
    auto builder = ElementBuilder::make_element_builder(eq_sys);
    auto point_locator = mesh.sub_point_locator();
    auto settings = this->context<ApecMSettingsObject>();
    const auto &axis = settings->target_axis();
    // Here you will be able to restrict to the subdomain of stranded conductors for speedup
    // You can also find before the element where the center is in pre processing
    const Elem *target_elem = point_locator->operator()(axis);
    auto &comm = mesh.comm();
    // target_axis(0) = -std::numeric_limits<Real>::max();
    // target_axis(1) = -std::numeric_limits<Real>::max();
    // target_axis(2) = -std::numeric_limits<Real>::max();
    // if (target_elem)
    // {
    builder->build_solution_vectors_from_emat_files(target_elem->id(), session->get_path("emat_folderpath"));
    // }
    b_target_axis.resize(2);
    b_target_axis(0) = nodal_solution_system.point_value(bx_var_number, axis, target_elem, nodal_solution_system.solution.get());
    b_target_axis(1) = nodal_solution_system.point_value(by_var_number, axis, target_elem, nodal_solution_system.solution.get());
}

void APECMEquilibriumAnalysis::remove_plasma_out_of_X_point()
{
    libMesh::out << "Removing plasma region outside X-point..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    // unsigned int psi_var_number = element_solution_system.variable_number("PSI");
    unsigned int bx_var_number = element_solution_system.variable_number("Be_x");
    unsigned int by_var_number = element_solution_system.variable_number("Be_y");
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("PLASMA");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");
    // dof_id_type vacuum_boundary_component = mesh.get_node_integer_index("VACUUM_CHAMBER_BOUNDARY");
    dof_id_type new_vacuum_boundary_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");

    auto settings_ptr = this->context<ApecMSettingsObject>();

    int current_sign = (settings_ptr->plasma_currents() >= 0) ? 1 : -1;

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
    libMesh::out << "Finding magnetic axis..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();

    dof_id_type vacuum_chamber_region_component = mesh.get_node_integer_index("VACUUM_CHAMBER_NODES");

    std::pair<dof_id_type, Real> magnetic_axis_data;

    auto &comm = mesh.comm();

    auto settings_ptr = this->context<ApecMSettingsObject>();
    min_distance_to_axis = settings_ptr->min_distance_to_axis();

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
    libMesh::out << "Identifying plasma region..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &mesh = eq_sys->get_mesh();

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("PLASMA");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");

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
    libMesh::out << "Checking limiter configuration..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    unsigned int bx_var_number = nodal_solution_system.variable_number("B_x");
    unsigned int by_var_number = nodal_solution_system.variable_number("B_y");
    const auto &mesh = eq_sys->get_mesh();
    Real psi_limiter = -std::numeric_limits<Real>::max();

    dof_id_type vacuum_boundary_region_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");

    auto settings_ptr = this->context<ApecMSettingsObject>();
    int current_sign = (settings_ptr->plasma_currents() >= 0) ? 1 : -1;

    for (const auto &node : mesh.local_node_ptr_range())
    {
        bool is_vacuum_chamber = node->get_extra_integer(vacuum_boundary_region_component);
        if (!is_vacuum_chamber)
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
    libMesh::out << "Checking X-point configuration..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();
    Real psi_X_point = std::numeric_limits<Real>::max();

    dof_id_type vacuum_chamber_region_component = mesh.get_node_integer_index("VACUUM_CHAMBER_NODES");
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
    libMesh::out << "Finding separatrix..." << std::endl;
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
    libMesh::out << "Computing correction currents..." << std::endl;
    auto settings = this->context<ApecMSettingsObject>();
    const Real bradx = settings->bradx();
    const Real brady = settings->brady();
    const Real bverx = settings->bverx();
    const Real bvery = settings->bvery();
    const Real cpgain = settings->cpgain();

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
    libMesh::out << "Applying control currents..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int vertical_coil_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &coils_current_vector = element_loads_system.get_vector("ramped_loads");

    const auto &mesh = eq_sys->get_mesh();
    dof_id_type vertical_control_coils_component = mesh.get_elem_integer_index("VERTICAL_CONTROL_COILS");
    dof_id_type radial_control_coils_component = mesh.get_elem_integer_index("RADIAL_CONTROL_COILS");
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

        char correction_current_sign = real_constants_map->get_constant(real_constant_elem_id)[22]; // CPSIGN
        dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, vertical_coil_current_var_number, 0);
        Real existing_current_density = element_loads_system.current_local_solution->el(elem_current_density_dof);
        Real new_current_density = existing_current_density + correction_current_sign * correction_current_density;
        coils_current_vector.set(elem_current_density_dof, new_current_density);
    }

    first_absolute_iteration = false;
}

void APECMEquilibriumAnalysis::apply_jplasma_gradshafranov()
{
    libMesh::out << "Applying J. Plasma Grad-Shafranov current profile..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto settings = this->context<ApecMSettingsObject>();
    const Real beta_profile = settings->beta_profile();
    const Real alpha = settings->alpha();
    const Real R0 = settings->major_radius();

    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");

    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int plasma_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &plasma_current_vector = element_loads_system.get_vector("ramped_loads");

    const auto &mesh = eq_sys->get_mesh();
    auto &comm = mesh.comm();
    dof_id_type plasma_component = mesh.get_elem_integer_index("PLASMA");

    auto settings_ptr = this->context<ApecMSettingsObject>();
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
                                              (beta_profile * elem_centroid_rad / R0 + (1 - beta_profile) * R0 / elem_centroid_rad);
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
        Real adjusted_current_density = existing_current_density * (settings_ptr->plasma_currents() / total_plasma_current);
        plasma_current_vector.set(elem_current_density_dof, adjusted_current_density);
    }
}

bool APECMEquilibriumAnalysis::check_poloidal_flux_convergence()
{
    libMesh::out << "Checking poloidal flux convergence..." << std::endl;
    auto eq_sys = this->active_equation_system();

    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");

    const auto &mesh = eq_sys->get_mesh();
    auto &comm = mesh.comm();
    dof_id_type vacuum_chamber_id = mesh.get_elem_integer_index("VACUUM_CHAMBER");

    auto settings = this->context<ApecMSettingsObject>();

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
    const Real tolerance = settings->tolerance();

    element_solution_system.current_local_solution->get(vacuum_chamber_psi_dofs, curr_it_vacuum_chamber_psi_vector);

    for (std::size_t i = 0; i < curr_it_vacuum_chamber_psi_vector.size(); ++i)
    {
        Real rel_diff = std::abs(curr_it_vacuum_chamber_psi_vector[i] - prev_it_vacuum_chamber_psi_vector[i]) / std::max(std::abs(curr_it_vacuum_chamber_psi_vector[i]), std::abs(prev_it_vacuum_chamber_psi_vector[i]));
        if (rel_diff > tolerance)
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
    {
        find_separatrix();
        apply_jplasma_gradshafranov();
    }
    bool converged = check_poloidal_flux_convergence();

    set_b_target_axis();

    if (converged)
        return true;

    auto correction_currents = compute_correction_currents();
    apply_control_currents(correction_currents);

    solve_equation_systems(substep);
    return false;
}

// void APECMEquilibriumAnalysis::constant_plasma_current_density()
// {
//     auto eq_sys = this->active_equation_system();

//     auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
//     unsigned int plasma_current_var_number = element_loads_system.variable_number("Js_z");
//     NumericVector<Number> &plasma_current_vector = element_loads_system.get_vector("ramped_loads");

//     const auto &mesh = eq_sys->get_mesh();
//     dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");

//     auto real_constants_map = this->context<RealConstantsObject>();
//     dof_id_type real_constant_mesh_id = mesh.get_elem_integer_index("ansys_real_constant");

//     auto settings_ptr = this->context<ApecMSettingsObject>();

//     for (auto &elem : mesh.active_local_element_ptr_range())
//     {
//         if (!elem->get_extra_integer(vacuum_chamber_component))
//             continue;

//         dof_id_type real_constant_elem_id = elem->get_extra_integer(real_constant_mesh_id);
//         Real plasma_area = real_constants_map->get_constant(real_constant_elem_id)[0];
//         Real plasma_current_density = settings_ptr->plasma_currents() / plasma_area;
//         dof_id_type elem_current_density_dof = elem->dof_number(Globals::SystemType::ELEMENT_LOADS, plasma_current_var_number, 0);
//         plasma_current_vector.set(elem_current_density_dof, plasma_current_density);
//     }
// }
