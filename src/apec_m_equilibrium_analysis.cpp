#include "apec_m_equilibrium_analysis.h"

#include "apec_m_settings.h"
#include "session.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <libmesh/libmesh.h>

void APECMEquilibriumAnalysis::append_iteration_diagnostics(Substep &substep,
                                                            bool converged)
{
    auto eq_sys = this->active_equation_system();
    if (!eq_sys)
        return;

    const auto &mesh = eq_sys->get_mesh();
    if (mesh.comm().rank() != 0)
        return;

    auto session = this->context<Session>();
    auto settings = this->context<ApecMSettingsObject>();
    if (!session || !settings)
        return;

    const std::filesystem::path log_path = session->get_path("time_steps_folderpath") / "apec_iteration_diagnostics.csv";
    std::filesystem::create_directories(log_path.parent_path());

    const bool write_header = !std::filesystem::exists(log_path) ||
                              std::filesystem::file_size(log_path) == 0;

    std::ofstream out(log_path, std::ios::app);
    if (!out)
        return;

    if (write_header)
    {
        out << "substep_time,picard_iteration,"
            << "converged,flux_metric,flux_threshold,"
            << "psi_magnetic_axis,psi_min,psi_limiter,psi_x_point,separatrix_configuration,"
            << "magnetic_axis_x,magnetic_axis_y,magnetic_axis_z,"
            << "limiter_x,limiter_y,limiter_z,x_point_x,x_point_y,x_point_z,"
            << "b_target_axis_x,b_target_axis_y,"
            << "radial_correction_current,vertical_correction_current,"
            << "plasma_current_target,beta_profile,"
            << "total_plasma_current_before_scaling,plasma_current_scaling_lambda,plasma_element_count"
            << '\n';
    }

    out << std::setprecision(17);
    const unsigned int picard_iteration = substep.get_number_picard_iterations() + 1;

    const auto point_or_sentinel = [](const Point &p, bool valid) -> Point
    {
        if (valid)
            return p;
        return Point(-std::numeric_limits<Real>::max(),
                     -std::numeric_limits<Real>::max(),
                     -std::numeric_limits<Real>::max());
    };

    const Point limiter_point = point_or_sentinel(last_limiter_point, has_last_limiter_point);
    const Point x_point = point_or_sentinel(last_x_point, has_last_x_point);

    out << substep.get_time() << ','
        << picard_iteration << ','
        << (converged ? 1 : 0) << ','
        << last_poloidal_flux_convergence_metric << ','
        << last_poloidal_flux_convergence_threshold << ','
        << psi_magnetic_axis << ','
        << psi_min << ','
        << last_psi_limiter << ','
        << last_psi_x_point << ','
        << last_separatrix_configuration << ','
        << magnetic_axis(0) << ','
        << magnetic_axis(1) << ','
        << magnetic_axis(2) << ','
        << limiter_point(0) << ','
        << limiter_point(1) << ','
        << limiter_point(2) << ','
        << x_point(0) << ','
        << x_point(1) << ','
        << x_point(2) << ','
        << b_target_axis(0) << ','
        << b_target_axis(1) << ','
        << last_radial_correction_current << ','
        << last_vertical_correction_current << ','
        << settings->plasma_currents() << ','
        << settings->beta_profile() << ','
        << last_total_plasma_current_before_scaling << ','
        << last_plasma_current_scaling_lambda << ','
        << last_plasma_element_count
        << '\n';
}

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
    const auto &axis = settings->target_axis();

    // Here you will be able to restrict to the subdomain of stranded conductors for speedup
    // You can also find before the element where the center is in pre processing
    auto &comm = mesh.comm();
    // const Elem *target_elem = point_locator->operator()(axis);

    b_target_axis(0) = 0.0;
    b_target_axis(1) = 0.0;
    // b_target_axis(2) = -std::numeric_limits<Real>::max();
    // if (target_elem)
    // {
    // builder->build_solution_vectors_from_emat_files(target_elem->id(), session->get_path("emat_folderpath"));
    // }
    b_target_axis(0) = nodal_solution_system.point_value(bx_var_number, axis, true);
    b_target_axis(1) = nodal_solution_system.point_value(by_var_number, axis, true);
    // }
    comm.sum(b_target_axis(0));
    comm.sum(b_target_axis(1));
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

        Real x_centroid = elem->vertex_average()(0);
        Real y_centroid = elem->vertex_average()(1);

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

    auto &support_solution_system = eq_sys->get_system<ExplicitSystem>("SolutionSupportSystem");
    unsigned int psi_var_number = support_solution_system.variable_number("PSI");

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(plasma_region_component) != 1)
            continue;

        for (auto &node : elem->node_ref_range())
        {
            dof_id_type psi_dof = node.dof_number(ApecMSystemBase::SUPPORT_SOLUTION, psi_var_number, 0);
            Real psi_value = support_solution_system.current_local_solution->el(psi_dof);
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
    auto &support_solution_system = eq_sys->get_system<ExplicitSystem>("SolutionSupportSystem");
    unsigned int psi_var_number = support_solution_system.variable_number("PSI");

    const auto &mesh = eq_sys->get_mesh();
    Real psi_limiter = -std::numeric_limits<Real>::max();

    dof_id_type vacuum_boundary_region_component = mesh.get_node_integer_index("NEW_VACUUM_CHAMBER_BOUNDARY");
    std::pair<dof_id_type, Real> limiter_data{libMesh::invalid_uint, -std::numeric_limits<Real>::max()};
    for (const auto &node : mesh.local_node_ptr_range())
    {
        if (node->get_extra_integer(vacuum_boundary_region_component) != 1)
            continue;
        dof_id_type psi_dof = node->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, psi_var_number, 0);
        Real psi_value = support_solution_system.current_local_solution->el(psi_dof);
        if (psi_value > psi_limiter)
        {
            psi_limiter = psi_value;
            limiter_data = std::make_pair(node->id(), psi_value);
        }
    }
    std::set<std::pair<dof_id_type, Real>> limiter_data_all;
    limiter_data_all.insert(limiter_data);
    mesh.comm().set_union(limiter_data_all);

    dof_id_type limiter_node = libMesh::invalid_uint;
    psi_limiter = -std::numeric_limits<Real>::max();
    for (const auto &data : limiter_data_all)
    {
        if (data.second > psi_limiter)
        {
            psi_limiter = data.second;
            limiter_node = data.first;
        }
    }

    Real lx = -std::numeric_limits<Real>::max();
    Real ly = -std::numeric_limits<Real>::max();
    Real lz = -std::numeric_limits<Real>::max();
    if (const Node *limiter_node_ptr = mesh.query_node_ptr(limiter_node))
    {
        lx = (*limiter_node_ptr)(0);
        ly = (*limiter_node_ptr)(1);
        lz = (*limiter_node_ptr)(2);
    }
    mesh.comm().max(lx);
    mesh.comm().max(ly);
    mesh.comm().max(lz);

    last_limiter_point = Point(lx, ly, lz);
    has_last_limiter_point = (limiter_node != libMesh::invalid_uint);

    last_psi_limiter = psi_limiter;

    return psi_limiter;
}

Real APECMEquilibriumAnalysis::check_X_point()
{
    libMesh::out << "Checking X-point configuration..." << std::endl;
    auto eq_sys = this->active_equation_system();
    auto &element_solution_system = eq_sys->get_system<ExplicitSystem>("ElementSolutionSystem");
    unsigned int bx_var_number = element_solution_system.variable_number("Be_x");
    unsigned int by_var_number = element_solution_system.variable_number("Be_y");
    auto &support_solution_system = eq_sys->get_system<ExplicitSystem>("SolutionSupportSystem");
    unsigned int psi_var_number = support_solution_system.variable_number("PSI");
    const auto &mesh = eq_sys->get_mesh();
    const auto proc_id = mesh.comm().rank();

    // 1) Find local candidate element with minimum |B| from element solution.
    // Real min_b_sum = std::numeric_limits<Real>::max();
    // dof_id_type x_point_elem_id = libMesh::invalid_uint;
    Real psi_X_point = std::numeric_limits<Real>::max();
    std::pair<dof_id_type, Real> x_point_data{libMesh::invalid_uint, std::numeric_limits<Real>::max()};

    Real min_b_sum = std::numeric_limits<Real>::max();

    dof_id_type vacuum_chamber_region_component = mesh.get_elem_integer_index("VACUUM_CHAMBER");

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->get_extra_integer(vacuum_chamber_region_component) != 1)
            continue;

        Real dist = (elem->vertex_average() - magnetic_axis).norm();
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
            x_point_data = std::make_pair(elem->id(), b_sum);
        }
    }

    // Global element selection by minimum field across processors.
    std::set<std::pair<dof_id_type, Real>> x_point_elem_data_all;
    x_point_elem_data_all.insert(x_point_data);
    mesh.comm().set_union(x_point_elem_data_all);

    min_b_sum = std::numeric_limits<Real>::max();
    dof_id_type x_point_elem_id = libMesh::invalid_uint;
    for (const auto &data : x_point_elem_data_all)
    {
        if (data.second < min_b_sum)
        {
            min_b_sum = data.second;
            x_point_elem_id = data.first;
        }
    }

    // 2-4) On the owning processor: accumulate element Be contributions on nodes and scan nodes.
    Real local_min_node_b_norm = std::numeric_limits<Real>::max();
    std::pair<dof_id_type, Real> x_point_data_nodes{libMesh::invalid_uint, std::numeric_limits<Real>::max()};
    const Elem *x_elem_ptr = mesh.query_elem_ptr(x_point_elem_id);
    if (x_elem_ptr && x_elem_ptr->processor_id() == proc_id)
    {
        std::vector<const Elem *> candidate_elements;
        std::unordered_set<dof_id_type> candidate_element_ids;
        candidate_elements.push_back(x_elem_ptr);
        candidate_element_ids.insert(x_elem_ptr->id());
        for (unsigned int k = 0; k < x_elem_ptr->n_neighbors(); ++k)
        {
            const Elem *neighbor = x_elem_ptr->neighbor_ptr(k);
            if (!neighbor)
                continue;
            if (!candidate_element_ids.insert(neighbor->id()).second)
                continue;
            candidate_elements.push_back(neighbor);
        }

        std::unordered_map<dof_id_type, std::pair<Real, Real>> node_b_sums;
        std::unordered_map<dof_id_type, unsigned int> node_counts;

        for (const Elem *elem_ptr : candidate_elements)
        {
            dof_id_type bx_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, bx_var_number, 0);
            dof_id_type by_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, by_var_number, 0);
            const Real bx_value = element_solution_system.current_local_solution->el(bx_dof);
            const Real by_value = element_solution_system.current_local_solution->el(by_dof);

            for (const auto &node : elem_ptr->node_ref_range())
            {
                auto &sum_pair = node_b_sums[node.id()];
                sum_pair.first += bx_value;
                sum_pair.second += by_value;
                node_counts[node.id()] += 1;
            }
        }

        for (const auto &node_sum_entry : node_b_sums)
        {
            const dof_id_type node_id = node_sum_entry.first;
            const auto count_it = node_counts.find(node_id);
            if (count_it == node_counts.end() || count_it->second == 0)
                continue;

            const Real inv_count = 1.0 / static_cast<Real>(count_it->second);
            const Real avg_bx = node_sum_entry.second.first * inv_count;
            const Real avg_by = node_sum_entry.second.second * inv_count;
            const Real b_norm = std::sqrt(avg_bx * avg_bx + avg_by * avg_by);
            if (b_norm < local_min_node_b_norm)
            {
                const Node *node_ptr = mesh.query_node_ptr(node_id);
                if (!node_ptr)
                    continue;

                local_min_node_b_norm = b_norm;
                dof_id_type psi_dof = node_ptr->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, psi_var_number, 0);
                x_point_data_nodes = std::make_pair(node_id, support_solution_system.current_local_solution->el(psi_dof));
            }
        }
    }

    std::set<std::pair<dof_id_type, Real>> x_point_node_data_all;
    x_point_node_data_all.insert(x_point_data_nodes);
    mesh.comm().set_union(x_point_node_data_all);

    dof_id_type x_point_node_id = libMesh::invalid_uint;
    for (const auto &data : x_point_node_data_all)
    {
        if (data.first != libMesh::invalid_uint)
        {
            psi_X_point = data.second;
            x_point_node_id = data.first;
        }
    }

    Real xpx = -std::numeric_limits<Real>::max();
    Real xpy = -std::numeric_limits<Real>::max();
    Real xpz = -std::numeric_limits<Real>::max();

    const Node *x_node_ptr = mesh.query_node_ptr(x_point_node_id);
    if (x_node_ptr)
    {
        xpx = x_node_ptr->operator()(0);
        xpy = x_node_ptr->operator()(1);
        xpz = x_node_ptr->operator()(2);
    }

    mesh.comm().max(xpx);
    mesh.comm().max(xpy);
    mesh.comm().max(xpz);

    last_x_point = Point(xpx, xpy, xpz);
    has_last_x_point = (x_point_node_id != libMesh::invalid_uint);
    last_psi_x_point = psi_X_point;

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
    {
        psi_min = psi_limiter;
        last_separatrix_configuration = "limiter";
    }
    else
    {
        psi_min = psi_X_point;
        last_separatrix_configuration = "x_point";
    }

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

    last_radial_correction_current = xirad;
    last_vertical_correction_current = xiver;
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
    const Real alpha_mp = settings->alpha_mp_();
    const Real alpha_np = settings->alpha_np_();
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
        Real elem_centroid_rad = elem->vertex_average()(0);

        Real psi_segn = (psi_elem - psi_magnetic_axis) / (psi_min - psi_magnetic_axis);
        Real grad_shafranov_current_density = std::pow((1 - std::pow(psi_segn, alpha_mp)), alpha_np) *
                                              (beta_profile * elem_centroid_rad / R0 + (1 - beta_profile) * R0 / elem_centroid_rad);
        // Real grad_shafranov_current_density = std::pow((1 - std::pow(psi_segn, alpha_mp)), alpha_np) *
        //                                           (beta_profile * elem_centroid_rad / R0) +
        //                                       std::pow((1 - std::pow(psi_segn, alpha_mp)), alpha_np + 1) * ((1 - beta_profile) * R0 / elem_centroid_rad);
        if (std::isnan(grad_shafranov_current_density))
            grad_shafranov_current_density = 0.0;
        total_plasma_current += grad_shafranov_current_density * elem->volume();
        plasma_current_vector.set(elem_current_density_dof, grad_shafranov_current_density);
        plasma_elem_count++;
    }

    // Take total current across all processors
    comm.sum(total_plasma_current);
    comm.sum(plasma_elem_count);

    last_total_plasma_current_before_scaling = total_plasma_current;
    last_plasma_element_count = plasma_elem_count;

    Real lambda = 0.0;
    if (std::abs(total_plasma_current) > std::numeric_limits<Real>::epsilon())
        lambda = settings_ptr->plasma_currents() / total_plasma_current;
    last_plasma_current_scaling_lambda = lambda;
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

    // if(0 == substep.get_number_picard_iterations())
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
        element_solution_system.solution->get(vacuum_chamber_psi_dofs, prev_it_vacuum_chamber_psi_vector);
        last_poloidal_flux_convergence_metric = std::numeric_limits<Real>::quiet_NaN();
        last_poloidal_flux_convergence_threshold = settings->tolerance();
        return false;
    }

    bool converged = true;
    const Real tolerance = settings->tolerance();
    last_poloidal_flux_convergence_threshold = tolerance;

    element_solution_system.solution->get(vacuum_chamber_psi_dofs, curr_it_vacuum_chamber_psi_vector);
    // Removed incorrect vectorized relative difference calculation.
    Real max_norm = 0.0;

    for (std::size_t i = 0; i < curr_it_vacuum_chamber_psi_vector.size(); ++i)
    {
        const Real denom = std::max(std::numeric_limits<Real>::epsilon(),
                                    std::min(std::abs(curr_it_vacuum_chamber_psi_vector[i]),
                                             std::abs(prev_it_vacuum_chamber_psi_vector[i])));
        // Real rel_diff = std::abs(curr_it_vacuum_chamber_psi_vector[i] - prev_it_vacuum_chamber_psi_vector[i]) / denom;
        Real rel_diff = std::abs(curr_it_vacuum_chamber_psi_vector[i] - prev_it_vacuum_chamber_psi_vector[i]);
        max_norm = std::max(max_norm, rel_diff);
    }
    comm.max(max_norm);
    last_poloidal_flux_convergence_metric = max_norm;

    if (max_norm > tolerance)
    {
        converged = false;
    }
    // if (std::abs(max_norm - prev_norm) < 0.1 * tolerance)
    // {
    //     converged = true;
    // }
    // // Check if converged in all processors
    // comm.min(converged);

    std::swap(prev_it_vacuum_chamber_psi_vector, curr_it_vacuum_chamber_psi_vector);
    prev_norm = max_norm;

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

    append_iteration_diagnostics(substep, converged);

    if (converged)
        return true;

    if (print_full_output)
        print_output(substep, true);

    first_absolute_iteration = false;
    return false;
}
