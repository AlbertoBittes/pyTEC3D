#include "apec_m_evolution_analysis.h"
#include "apec_m.h"

#include <filesystem>

#include <libmesh/libmesh.h>

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : APECMEquilibriumAnalysis(program, context_view),
      toroidal_projection(std::make_unique<ToroidalProjection>()),
      cylindrical_coordinate_system(CylindricalCoordinateSystem())
{
}

bool APECMEvolutionAnalysis::picard_loop(Substep &substep)
{
    auto settings = this->context<ApecMSettingsObject>();
    bool print_full_output = settings->print_full_output();

    APEC_M *program = this->get_parent_object();

    const std::string &main_mesh = settings->main_mesh();
    const std::string &static_mesh = settings->main_mesh();

    auto main_eq_sys = program->get_equation_system(main_mesh);
    auto static_eq_sys = program->get_equation_system(static_mesh);

    ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &main_sys_nodal_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::NODAL_SOLUTION);
    ApecMExplicitSystem &static_sys_nodal_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::NODAL_SOLUTION);
    ApecMExplicitSystem &main_sys_element_loads = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    program->set_active_mesh(main_mesh);

    auto toroidal_proj = toroidal_projection.get();

    unsigned int main_mesh_dim = main_eq_sys->get_mesh().mesh_dimension();

    const std::string jtor_lab = (main_mesh_dim == 2) ? "_z" : "_y";
    const std::string brad_lab = "_x";
    const std::string bver_lab = (main_mesh_dim == 2) ? "_y" : "_z";

    if (!first_absolute_iteration)
    {
        libMesh::out << "Projecting from static mesh to main mesh..." << std::endl;
        toroidal_proj->project(
            *static_eq_sys,
            *main_eq_sys,
            "Jt_z",
            *static_sys_element_sol.current_local_solution.get(),
            "Js" + jtor_lab,
            *main_sys_element_loads.rhs,
            cylindrical_coordinate_system);
        libMesh::out << "Projection completed." << std::endl;
        libMesh::out << "Rotating to global coordinates..." << std::endl;
        main_sys_element_loads.local_to_global_coordinate_system_solution(
            cylindrical_coordinate_system,
            "Js");
    }

    libMesh::out << "Solving evolution step..." << std::endl;
    solve_equation_systems(substep, true);

    if (print_full_output)
        print_output(substep, true);

    program->set_active_mesh(static_mesh);

    libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    main_sys_nodal_sol.global_to_local_coordinate_system_solution(
        cylindrical_coordinate_system,
        "B");
    main_sys_element_sol.global_to_local_coordinate_system_solution(
        cylindrical_coordinate_system,
        "Jt");

    // Project back to the static mesh
    // Project J, from plasma, coils and vessel
    toroidal_proj->toroidal_average(
        *main_eq_sys,
        *static_eq_sys,
        "Jt" + jtor_lab,
        *main_sys_element_loads.rhs,
        "Js_z",
        *static_sys_element_loads.current_local_solution.get(),
        cylindrical_coordinate_system);

    // Project PSI
    toroidal_proj->toroidal_average(
        *main_eq_sys,
        *static_eq_sys,
        "PSI",
        *main_sys_nodal_sol.current_local_solution.get(),
        "PSI",
        *static_sys_nodal_sol.current_local_solution.get(),
        cylindrical_coordinate_system);

    toroidal_proj->toroidal_average(
        *main_eq_sys,
        *static_eq_sys,
        "PSIe",
        *main_sys_element_sol.current_local_solution.get(),
        "PSIe",
        *static_sys_element_sol.current_local_solution.get(),
        cylindrical_coordinate_system);

    // Project B
    toroidal_proj->toroidal_average(
        *main_eq_sys,
        *static_eq_sys,
        "B" + brad_lab,
        *main_sys_nodal_sol.rhs,
        "B_y",
        *static_sys_nodal_sol.current_local_solution.get(),
        cylindrical_coordinate_system);

    toroidal_proj->toroidal_average(
        *main_eq_sys,
        *static_eq_sys,
        "B" + bver_lab,
        *main_sys_nodal_sol.rhs,
        "B_y",
        *static_sys_nodal_sol.current_local_solution.get(),
        cylindrical_coordinate_system);

    find_separatrix();

    apply_jplasma_gradshafranov();

    bool converged = check_poloidal_flux_convergence();

    if (converged)
        return true;

    libMesh::out << "Solving equilibrium step..." << std::endl;
    solve_equation_systems(substep, true);
    if (print_full_output)
        print_output(substep, true);

    first_absolute_iteration = false;
    return false;
}
