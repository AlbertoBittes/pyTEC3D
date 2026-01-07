#include "apec_m_equilibrium_analysis.h"

#include "apec_m_settings.h"
#include "session.h"

#include <algorithm>
#include <libmesh/libmesh.h>

bool APECMEquilibriumAnalysis::run_picard_iterations(Substep &substep)
{
    auto settings = this->context<ApecMSettingsObject>();
    const auto iterations_limit = std::max<unsigned int>(1u, settings->max_picard_iterations());
    bool converged = false;

    substep.reset_number_picard_iterations(0);
    for (unsigned int iteration = 0; iteration < iterations_limit; ++iteration)
    {
        if (picard_loop(substep))
        {
            substep.increment_number_picard_iterations();
            libMesh::out << "APEC algorithm converged in " << iteration + 1 << " Picard iterations." << std::endl;
            converged = true;
            break;
        }
        substep.increment_number_picard_iterations();
    }

    if (!converged)
    {
        libMesh::out << "APEC algorithm did not converge within " << iterations_limit
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

    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int plasma_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &plasma_current_vector = *element_loads_system.solution;

    dof_id_type plasma_region_component = mesh.get_elem_integer_index("PLASMA");
    dof_id_type vacuum_chamber_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");
    dof_id_type vacuum_boundary_component = mesh.get_node_integer_index("VACUUM_CHAMBER_BOUNDARY");
    dof_id_type new_vacuum_boundary_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");

    for (auto &elem : mesh.active_local_element_ptr_range())
        if (elem->get_extra_integer(vacuum_chamber_component) == 1)
        {
            elem->set_extra_integer(plasma_region_component, 1);

            dof_id_type elem_current_density_dof = elem->dof_number(ApecMSystemBase::ELEMENT_LOADS, plasma_current_var_number, 0);

            plasma_current_vector.set(elem_current_density_dof, 0.0);

            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(vacuum_boundary_component) == 1)
                    node.set_extra_integer(new_vacuum_boundary_component, 1);
        }

    plasma_current_vector.close();
    element_loads_system.update();
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
    // auto builder = ElementBuilder::make_element_builder(eq_sys);
    // auto point_locator = mesh.sub_point_locator();
    // point_locator->enable_out_of_mesh_mode();
    auto settings = this->context<ApecMSettingsObject>();
    dof_id_type target_axis_element_id = settings->target_axis_element_id();
    const auto &axis = settings->target_axis();
    // Here you will be able to restrict to the subdomain of stranded conductors for speedup
    // You can also find before the element where the center is in pre processing
    auto &comm = mesh.comm();
    // const Elem *target_elem = point_locator->operator()(axis);

    b_target_axis(0) = -std::numeric_limits<Real>::max();
    b_target_axis(1) = -std::numeric_limits<Real>::max();
    // b_target_axis(2) = -std::numeric_limits<Real>::max();
    // if (target_elem)
    // {
    // builder->build_solution_vectors_from_emat_files(target_elem->id(), session->get_path("emat_folderpath"));
    // }
    b_target_axis(0) = nodal_solution_system.point_value(bx_var_number, axis, target_axis_element_id, nodal_solution_system.solution.get());
    b_target_axis(1) = nodal_solution_system.point_value(by_var_number, axis, target_axis_element_id, nodal_solution_system.solution.get());
    // }
    comm.max(b_target_axis(0));
    comm.max(b_target_axis(1));
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
        if (elem->get_extra_integer(vacuum_chamber_component) != 1)
            continue;

        dof_id_type bx_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, bx_var_number, 0);
        dof_id_type by_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, by_var_number, 0);
        Real bx_value = element_solution_system.current_local_solution->el(bx_dof);
        Real by_value = element_solution_system.current_local_solution->el(by_dof);

        Real x_centroid = elem->centroid()(0);
        Real y_centroid = elem->centroid()(1);

        Real b_rot = current_sign * (by_value * (x_centroid - magnetic_axis(0)) - bx_value * (y_centroid - magnetic_axis(1)));

        if (b_rot >= 0)
        {
            elem->set_extra_integer(plasma_region_component, 0);
            for (auto &node : elem->node_ref_range())
                if (node.get_extra_integer(new_vacuum_boundary_component) == 1)
                    node.set_extra_integer(new_vacuum_boundary_component, 0);
        }
        else
        {
            continue;
        }
    }
}

void APECMEquilibriumAnalysis::set_psi_magnetic_axis()
{
    libMesh::out << "Finding magnetic axis..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &support_solution_system = eq_sys->get_system<ExplicitSystem>("SolutionSupportSystem");
    unsigned int psi_var_number = support_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();

    dof_id_type vacuum_chamber_region_component = mesh.get_node_integer_index("VACUUM_CHAMBER_NODES");

    std::pair<dof_id_type, Real> magnetic_axis_data;

    auto &comm = mesh.comm();

    auto settings_ptr = this->context<ApecMSettingsObject>();
    min_distance_to_axis = settings_ptr->min_distance_to_axis();

    for (const auto &node : mesh.local_node_ptr_range())
    {
        if (node->get_extra_integer(vacuum_chamber_region_component) != 1)
            continue;

        // if (!check_distance_magnetic_axis(node))
        //     continue;
        // Process the element as needed
        dof_id_type psi_dof = node->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, psi_var_number, 0);
        Real psi_value = support_solution_system.current_local_solution->el(psi_dof);
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

    // auto &nodal_solution_system = eq_sys->get_system<ExplicitSystem>("NodalSolutionSystem");
    // unsigned int psi_var_number = nodal_solution_system.variable_number("PSI");

    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(plasma_region_component) != 1)
            continue;

        dof_id_type psi_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, psi_var_number, 0);
        Real psi_value = element_solution_system.current_local_solution->el(psi_dof);
        if (psi_value < psi_min)
            elem->set_extra_integer(plasma_region_component, 0);

        // for (auto &node : elem->node_ref_range())
        // {
        //     dof_id_type psi_dof = node.dof_number(ApecMSystemBase::NODAL_SOLUTION, psi_var_number, 0);
        //     Real psi_value = nodal_solution_system.current_local_solution->el(psi_dof);
        //     if (psi_value < psi_min)
        //     {
        //         elem->set_extra_integer(plasma_region_component, 0);
        //         break;
        //     }
        // }
    }
}

Real APECMEquilibriumAnalysis::check_limiter()
{
    libMesh::out << "Checking limiter configuration..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &support_solution_system = eq_sys->get_system<ExplicitSystem>("SolutionSupportSystem");
    unsigned int psi_var_number = support_solution_system.variable_number("PSI");

    const auto &mesh = eq_sys->get_mesh();
    Real psi_limiter = -std::numeric_limits<Real>::max();

    dof_id_type vacuum_boundary_region_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");

    auto settings_ptr = this->context<ApecMSettingsObject>();
    int current_sign = (settings_ptr->plasma_currents() >= 0) ? 1 : -1;

    unsigned int limiter_node = 0;
    for (const auto &node : mesh.local_node_ptr_range())
    {
        if (node->get_extra_integer(vacuum_boundary_region_component) != 1)
            continue;
        // Process the element as needed
        dof_id_type psi_dof = node->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, psi_var_number, 0);
        Real psi_value = support_solution_system.current_local_solution->el(psi_dof);
        if (psi_value > psi_limiter)
        {
            limiter_node = node->id();
            psi_limiter = psi_value;
        }
    }
    mesh.comm().max(psi_limiter);
    return psi_limiter;
}

Real APECMEquilibriumAnalysis::check_X_point()
{
    libMesh::out << "Checking X-point configuration..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int psi_var_number = element_solution_system.variable_number("PSIe");
    unsigned int bx_var_number = element_solution_system.variable_number("Be_x");
    unsigned int by_var_number = element_solution_system.variable_number("Be_y");
    const auto &mesh = eq_sys->get_mesh();
    Real min_b_sum = std::numeric_limits<Real>::max();
    Real psi_X_point;

    dof_id_type vacuum_chamber_region_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");

    unsigned int X_point_elem = 0;
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(vacuum_chamber_region_component) != 1)
            continue;

        Real dist = (elem->centroid() - magnetic_axis).norm();
        if (dist < min_distance_to_axis)
            continue;
        // Process the element as needed
        dof_id_type bx_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, bx_var_number, 0);
        dof_id_type by_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, by_var_number, 0);
        Real bx_value = element_solution_system.current_local_solution->el(bx_dof);
        Real by_value = element_solution_system.current_local_solution->el(by_dof);
        Real b_sum = std::sqrt(bx_value * bx_value + by_value * by_value);
        if (b_sum < min_b_sum)
        {
            min_b_sum = b_sum;
            dof_id_type psi_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, psi_var_number, 0);
            psi_X_point = element_solution_system.current_local_solution->el(psi_dof);
            X_point_elem = elem->id();
        }
    }
    mesh.comm().min(psi_X_point);
    return psi_X_point;
}

void APECMEquilibriumAnalysis::find_separatrix()
{
    libMesh::out << "Finding separatrix..." << std::endl;
    this->reset_plasma_components();
    this->set_psi_magnetic_axis();
    this->remove_plasma_out_of_X_point();

    Real psi_limiter = this->check_limiter();
    Real psi_X_point = this->check_X_point();
    if (psi_X_point < psi_limiter)
        psi_min = psi_limiter;
    else
        psi_min = psi_X_point;

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

    // libMesh::out << "Radial correction current: " << xirad << " A" << std::endl;
    // libMesh::out << "Vertical correction current: " << xiver << " A" << std::endl;

    return std::make_pair(xirad, xiver);
}

void APECMEquilibriumAnalysis::apply_control_currents(std::pair<Real, Real> &correction_currents)
{
    libMesh::out << "Applying control currents..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &element_loads_system = eq_sys->get_system<ExplicitSystem>("ElementLoadsSystem");
    unsigned int vertical_coil_current_var_number = element_loads_system.variable_number("Js_z");
    NumericVector<Number> &coils_current_vector = *element_loads_system.solution;

    const auto &mesh = eq_sys->get_mesh();
    dof_id_type vertical_control_coils_component = mesh.get_elem_integer_index("VERTICAL_CONTROL_COILS");
    dof_id_type radial_control_coils_component = mesh.get_elem_integer_index("RADIAL_CONTROL_COILS");
    dof_id_type real_constant_mesh_id = mesh.get_elem_integer_index("ansys_real_constant");

    auto real_constants_map = this->context<RealConstantsObject>();

    Real correction_current;

    for (auto &elem : eq_sys->get_mesh().active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(vertical_control_coils_component) == 1)
            correction_current = correction_currents.second;
        else if (elem->get_extra_integer(radial_control_coils_component) == 1)
            correction_current = correction_currents.first;
        else
            continue;

        dof_id_type real_constant_elem_id = elem->get_extra_integer(real_constant_mesh_id);
        Real coil_area = real_constants_map->get_constant(real_constant_elem_id)[1];
        Real correction_current_density = correction_current / coil_area;

        Real correction_current_sign = real_constants_map->get_constant(real_constant_elem_id)[19]; // CPSIGN
        dof_id_type elem_current_density_dof = elem->dof_number(ApecMSystemBase::ELEMENT_LOADS, vertical_coil_current_var_number, 0);

        Real base_current = real_constants_map->get_constant(real_constant_elem_id)[20]; // CURR
        Real new_current_density = base_current / coil_area + correction_current_sign * correction_current_density;
        coils_current_vector.set(elem_current_density_dof, new_current_density);
    }
    coils_current_vector.close();
    element_loads_system.update();
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
    NumericVector<Number> &plasma_current_vector = *element_loads_system.solution;

    const auto &mesh = eq_sys->get_mesh();
    auto &comm = mesh.comm();
    dof_id_type plasma_component = mesh.get_elem_integer_index("PLASMA");

    auto settings_ptr = this->context<ApecMSettingsObject>();
    Real total_plasma_current = 0.0;

    unsigned int plasma_elem_count = 0;
    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(plasma_component) != 1)
            continue;

        dof_id_type elem_current_density_dof = elem->dof_number(ApecMSystemBase::ELEMENT_LOADS, plasma_current_var_number, 0);
        dof_id_type psi_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, psi_var_number, 0);

        Real psi_elem = element_solution_system.current_local_solution->el(psi_dof);
        Real elem_centroid_rad = elem->centroid()(0);

        Real psi_segn = (psi_elem - psi_magnetic_axis) / (psi_min - psi_magnetic_axis);
        Real grad_shafranov_current_density = (1 - std::pow(psi_segn, alpha)) *
                                              (beta_profile * elem_centroid_rad / R0 + (1 - beta_profile) * R0 / elem_centroid_rad);
        total_plasma_current += grad_shafranov_current_density * elem->volume();
        plasma_current_vector.set(elem_current_density_dof, grad_shafranov_current_density);
        plasma_elem_count++;
    }

    // Take total current across all processors
    comm.sum(total_plasma_current);

    Real lambda = settings_ptr->plasma_currents() / total_plasma_current;
    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(plasma_component) != 1)
            continue;

        dof_id_type elem_current_density_dof = elem->dof_number(ApecMSystemBase::ELEMENT_LOADS, plasma_current_var_number, 0);
        Real existing_current_density = plasma_current_vector.el(elem_current_density_dof);
        Real adjusted_current_density = existing_current_density * lambda;
        plasma_current_vector.set(elem_current_density_dof, adjusted_current_density);
    }

    plasma_current_vector.close();
    element_loads_system.update();
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
            if (elem->get_extra_integer(vacuum_chamber_id) == 1)
            {
                dof_id_type psi_dof = elem->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, psi_var_number, 0);
                vacuum_chamber_psi_dofs.push_back(psi_dof);
            }
        }
        prev_it_vacuum_chamber_psi_vector.resize(vacuum_chamber_psi_dofs.size(), 0.0);
        return false;
    }

    bool converged = true;
    const Real tolerance = settings->tolerance();

    element_solution_system.solution->get(vacuum_chamber_psi_dofs, curr_it_vacuum_chamber_psi_vector);
    // Removed incorrect vectorized relative difference calculation.
    Real tot_norm = 0.0;
    for (std::size_t i = 0; i < curr_it_vacuum_chamber_psi_vector.size(); ++i)
    {
        Real rel_diff = std::abs(curr_it_vacuum_chamber_psi_vector[i] - prev_it_vacuum_chamber_psi_vector[i]) / std::min(std::abs(curr_it_vacuum_chamber_psi_vector[i]), std::abs(prev_it_vacuum_chamber_psi_vector[i]));
        tot_norm += rel_diff;
    }

    if (tot_norm > tolerance * curr_it_vacuum_chamber_psi_vector.size())
    {
        converged = false;
    }
    // Check if converged in all processors
    comm.min(converged);

    std::swap(prev_it_vacuum_chamber_psi_vector, curr_it_vacuum_chamber_psi_vector);

    return converged;
}

bool APECMEquilibriumAnalysis::picard_loop(Substep &substep)
{
    auto settings = this->context<ApecMSettingsObject>();
    bool print_full_output = settings->print_full_output();

    if (!first_absolute_iteration)
    {
        find_separatrix();
        apply_jplasma_gradshafranov();
        set_b_target_axis();
        auto correction_currents = compute_correction_currents();
        apply_control_currents(correction_currents);
    }

    solve_equation_systems(substep, true);

    bool converged = check_poloidal_flux_convergence();

    if (converged)
        return true;

    if (print_full_output)
        print_output(substep, true);

    first_absolute_iteration = false;
    return false;
}
