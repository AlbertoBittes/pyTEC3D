#include "apec_m_evolution_analysis.h"
#include "apec_m.h"

#include <filesystem>
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
      toroidal_projection(std::make_unique<ToroidalProjection>()),
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

void APECMEvolutionAnalysis::set_toroidal_projection_map()
{
    auto settings = this->context<ApecMSettingsObject>();

    const auto &toroidal_map_nodes = settings->toroidal_projection_map_nodes();
    const auto &toroidal_map_elems = settings->toroidal_projection_map_elems();

    ToroidalProjection::Map nodal_map, elem_map;

    if (!toroidal_map_nodes.empty())
    {
        nodal_map = toroidal_map_nodes;
    }
    if (!toroidal_map_elems.empty())
    {
        elem_map = toroidal_map_elems;
    }

    toroidal_projection->set_maps(nodal_map, elem_map);
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

    ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);

    auto toroidal_proj = toroidal_projection.get();

    auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());

    {
        ApecMFieldRotation element_sol_rotation(main_sys_element_sol,
                                                cylindrical_coordinate_system,
                                                {"Be", "Jt"});

        // Project back to the static mesh
        // Project J, from plasma, coils and vessel
        // Project PSI
        const std::pair<std::string, std::string> toroidal_avg_pairs[] = {
            {"PSI", "PSI"},
            {"Jt" + j_b_labels_evol[0], "Js" + j_b_labels_stat[0]},
            {"PSIe", "PSIe"},
            {"Be" + j_b_labels_evol[1], "Be" + j_b_labels_stat[1]},
            {"Be" + j_b_labels_evol[2], "Be" + j_b_labels_stat[2]},
        };

        for (const auto &toroidal_avg_pair : toroidal_avg_pairs)
            toroidal_proj->toroidal_average(
                *main_eq_sys,
                *static_eq_sys,
                toroidal_avg_pair.first,
                toroidal_avg_pair.second,
                cylindrical_coordinate_system);
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

    auto toroidal_proj = toroidal_projection.get();

    auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());

    {
        ApecMFieldRotation element_loads_rotation(main_sys_element_loads,
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
    auto settings = this->context<ApecMSettingsObject>();
    bool print_full_output = settings->print_full_output();

    auto toroidal_proj = toroidal_projection.get();

    if (first_absolute_iteration)
        set_toroidal_projection_map();

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh = settings->evolution_mesh();
    const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    program->set_active_mesh(evolution_mesh);

    if (!first_absolute_iteration)
    {
        toroidal_projection_current();
        set_poloidal_currents();
    }

    libMesh::out << "Solving evolution step..." << std::endl;
    solve_equation_systems(substep, false);

    if (print_full_output)
        print_output(substep, false);

    program->set_active_mesh(equilibrium_mesh);

    toroidal_average(); // With this you set the averaged psi, eddy currents, and B fields on the static mesh

    find_separatrix();

    apply_jplasma_gradshafranov();

    libMesh::out << "Solving equilibrium step..." << std::endl;
    solve_equation_systems(substep, false);

    bool converged = check_poloidal_flux_convergence();

    if (converged)
        return true;

    if (print_full_output)
        print_output(substep, false);

    first_absolute_iteration = false;
    return false;
}
