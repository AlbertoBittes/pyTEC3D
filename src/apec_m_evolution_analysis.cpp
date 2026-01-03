#include "apec_m_evolution_analysis.h"
#include "apec_m.h"

#include <filesystem>

#include <libmesh/libmesh.h>

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : APECMEquilibriumAnalysis(program, context_view)
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

    program->set_active_mesh(main_mesh);

    auto toroidal_proj = toroidal_projection.get();

    if (!first_absolute_iteration)
    {
        reset_plasma_components();
        toroidal_proj->project(
            *static_eq_sys,
            *main_eq_sys,
            "Jt_z",
            *static_eq_sys->get_system(ApecMSystemBase::ELEMENT_SOLUTION).current_local_solution.get(),
            "Js_z",
            *main_eq_sys->get_system(ApecMSystemBase::ELEMENT_LOADS).current_local_solution.get());
    }

    solve_equation_systems(substep, true);

    if (print_full_output)
        print_output(substep, true);

    program->set_active_mesh(static_mesh);

    toroidal_proj->project_inverse(
        *main_eq_sys,
        *static_eq_sys,
        "Jt_z",
        *main_eq_sys->get_system(ApecMSystemBase::ELEMENT_SOLUTION).current_local_solution.get(),
        "Js_z",
        *static_eq_sys->get_system(ApecMSystemBase::ELEMENT_LOADS).current_local_solution.get());

    toroidal_proj->project_inverse(
        *main_eq_sys,
        *static_eq_sys,
        "PSI",
        *main_eq_sys->get_system(ApecMSystemBase::ELEMENT_SOLUTION).current_local_solution.get(),
        "PSI",
        *static_eq_sys->get_system(ApecMSystemBase::ELEMENT_LOADS).current_local_solution.get());

    toroidal_proj->project_inverse(
        *main_eq_sys,
        *static_eq_sys,
        "Jt_z",
        *main_eq_sys->get_system(ApecMSystemBase::ELEMENT_SOLUTION).current_local_solution.get(),
        "Js_z",
        *static_eq_sys->get_system(ApecMSystemBase::ELEMENT_LOADS).current_local_solution.get());

    find_separatrix();
    apply_jplasma_gradshafranov();

    bool converged = check_poloidal_flux_convergence();

    if (converged)
        return true;

    solve_equation_systems(substep, true);
    if (print_full_output)
        print_output(substep, true);

    first_absolute_iteration = false;
    return false;
}
