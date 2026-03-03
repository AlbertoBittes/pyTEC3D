#include "apec_m_evolution_analysis.h"
#include "apec_m.h"

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>

#include <libmesh/libmesh.h>

namespace
{
    std::vector<std::string> get_j_b_labels(unsigned int mesh_dim)
    {
        if (mesh_dim == 2)
        {
            return {"_z", "_x", "_y"};
        }
        else
        {
            return {"_y", "_x", "_z"};
        }
    }
} // namespace

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : APECMEquilibriumAnalysis(program, context_view),
      cylindrical_coordinate_system(CylindricalCoordinateSystem())
{
}

void APECMEvolutionAnalysis::set_poloidal_currents()
{
    // libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // // Rotate the J and B to cylindrical coordinates
    // auto settings = this->context<ApecMSettingsObject>();

    // APEC_M *program = this->get_parent_object();

    // const std::string &evolution_mesh = settings->evolution_mesh();
    // const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    // auto main_eq_sys = program->get_equation_system(evolution_mesh);
    // auto static_eq_sys = program->get_equation_system(equilibrium_mesh);

    // ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    // ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    // ApecMExplicitSystem &main_sys_element_loads = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    // ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    // auto toroidal_proj = toroidal_projection.get();

    // auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    // auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());
}

void APECMEvolutionAnalysis::set_toroidal_slices()
{
    auto settings = this->context<ApecMSettingsObject>();
    if (!settings)
    {
        toroidal_slices.clear();
        return;
    }

    const unsigned int n_slices = settings->has_column("n_slices")
                                      ? static_cast<unsigned int>(settings->get_int("n_slices"))
                                      : 0;
    toroidal_slices.assign(n_slices, {});
    if (n_slices == 0)
        return;

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh = settings->evolution_mesh();
    auto evolution_eq_sys = program->get_equation_system(evolution_mesh);
    if (!evolution_eq_sys)
        return;

    const MeshBase &mesh = evolution_eq_sys->get_mesh();

    std::vector<dof_id_type> slice_integer_indices(n_slices, libMesh::invalid_uint);
    for (unsigned int slice = 0; slice < n_slices; ++slice)
    {
        const std::string integer_name = "slice_" + std::to_string(slice);
        if (mesh.has_node_integer(integer_name))
            slice_integer_indices[slice] = mesh.get_node_integer_index(integer_name);
    }

    const dof_id_type slice0_integer_index = slice_integer_indices[0];
    if (slice0_integer_index == libMesh::invalid_uint)
        return;

    std::vector<dof_id_type> slice0_nodes;
    std::unordered_map<dof_id_type, std::size_t> slice0_to_row;
    for (const auto &node_ptr : mesh.node_ptr_range())
    {
        const long long slice0_node_id_signed = node_ptr->get_extra_datum<int>(slice0_integer_index);
        if (slice0_node_id_signed < 0)
            continue;

        const dof_id_type slice0_node_id = static_cast<dof_id_type>(slice0_node_id_signed);
        if (!mesh.query_node_ptr(slice0_node_id))
            continue;

        if (slice0_to_row.find(slice0_node_id) == slice0_to_row.end())
        {
            const std::size_t row = slice0_nodes.size();
            slice0_nodes.push_back(slice0_node_id);
            slice0_to_row.emplace(slice0_node_id, row);
        }
    }

    for (unsigned int slice = 0; slice < n_slices; ++slice)
        toroidal_slices[slice].assign(slice0_nodes.size(), libMesh::invalid_uint);

    for (std::size_t row = 0; row < slice0_nodes.size(); ++row)
        toroidal_slices[0][row] = slice0_nodes[row];

    for (const auto &node_ptr : mesh.node_ptr_range())
    {
        for (unsigned int slice = 1; slice < n_slices; ++slice)
        {
            const dof_id_type integer_index = slice_integer_indices[slice];
            if (integer_index == libMesh::invalid_uint)
                continue;

            const long long slice0_node_id_signed = node_ptr->get_extra_datum<int>(integer_index);
            if (slice0_node_id_signed < 0)
                continue;

            const dof_id_type slice0_node_id = static_cast<dof_id_type>(slice0_node_id_signed);
            const auto row_it = slice0_to_row.find(slice0_node_id);
            if (row_it == slice0_to_row.end())
                continue;

            toroidal_slices[slice][row_it->second] = node_ptr->id();
        }
    }
}

void APECMEvolutionAnalysis::toroidal_average()
{
    libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh = settings->evolution_mesh();
    const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    auto main_eq_sys = program->get_equation_system(evolution_mesh);
    auto static_eq_sys = program->get_equation_system(equilibrium_mesh);

    ApecMLinearSystem &main_sys = main_eq_sys->get_system<ApecMLinearSystem>(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    ApecMLinearSystem &static_sys = static_eq_sys->get_system<ApecMLinearSystem>(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());

    {
        VectorFieldCoordinateTransformation element_sol_rotation(main_sys_element_sol,
                                                                 cylindrical_coordinate_system,
                                                                 {"Be", "Jt"});
        VectorFieldCoordinateTransformation static_element_sol_rotation(static_sys_element_sol,
                                                                        cylindrical_coordinate_system,
                                                                        {"Be", "Jt"});

        // Project back to the static mesh
        // Project J, from plasma, coils and vessel
        // Project PSI
        const std::pair<std::string, std::string> toroidal_avg_pairs[] = {
            {"PSI", "PSI"},
            {"Bn" + j_b_labels_evol[1], "Bn" + j_b_labels_stat[1]},
            {"Bn" + j_b_labels_evol[2], "Bn" + j_b_labels_stat[2]},
        };
    }
}

void APECMEvolutionAnalysis::toroidal_projection_current()
{
    libMesh::out << "Projection from static mesh to main mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh = settings->evolution_mesh();
    const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    auto main_eq_sys = program->get_equation_system(evolution_mesh);
    auto static_eq_sys = program->get_equation_system(equilibrium_mesh);

    ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &main_sys_element_loads = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());

    {
        VectorFieldCoordinateTransformation static_element_loads_rotation(static_sys_element_loads,
                                                                          cylindrical_coordinate_system,
                                                                          {"Jt"});
        VectorFieldCoordinateTransformation element_loads_rotation(main_sys_element_loads,
                                                                   cylindrical_coordinate_system,
                                                                   {"Js"});

        // Project back to the static mesh
        // Project J, from plasma, coils and vessel
        const std::pair<std::string, std::string> toroidal_proj_pairs[] = {
            {"Jt" + j_b_labels_stat[0], "Js" + j_b_labels_evol[0]},
        };

        for (const auto &toroidal_proj_pair : toroidal_proj_pairs)
            toroidal_proj->project(
                *static_eq_sys,
                *main_eq_sys,
                toroidal_proj_pair.first,
                toroidal_proj_pair.second,
                cylindrical_coordinate_system);
    }
}

bool APECMEvolutionAnalysis::picard_loop(Substep &substep)
{
    APEC_M *program = this->get_parent_object();

    auto settings = this->context<ApecMSettingsObject>();
    bool print_full_output = settings->print_full_output();

    const std::string &evolution_mesh = settings->evolution_mesh();
    const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    if (first_absolute_iteration)
        set_toroidal_slices();

    if (substep.get_parent_object()->get_previous_time_step() == nullptr)
    {
        program->set_active_mesh(equilibrium_mesh);
        bool conv = APECMEquilibriumAnalysis::picard_loop(substep);
        program->set_active_mesh(evolution_mesh);
        toroidal_projection_current();
        solve_equation_systems(substep, false);
        return conv;
    }

    program->set_active_mesh(evolution_mesh);

    toroidal_projection_current();

    set_poloidal_currents();

    libMesh::out << "Solving evolution step..." << std::endl;
    solve_equation_systems(substep, false);

    if (print_full_output)
        print_output(substep, true);

    program->set_active_mesh(equilibrium_mesh);

    toroidal_average(); // With this you set the averaged psi, eddy currents, and B fields on the static mesh

    find_separatrix();

    apply_jplasma_gradshafranov();

    libMesh::out << "Solving equilibrium step..." << std::endl;
    solve_equation_systems(substep, false);

    bool converged = check_poloidal_flux_convergence();

    append_iteration_diagnostics(substep, converged);

    if (converged)
        return true;

    if (print_full_output)
        print_output(substep, true);

    return false;
}
