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
    const std::string &main_mesh = settings->main_mesh();

    APEC_M *program = this->get_parent_object();
    program->set_active_mesh(main_mesh);

    apply_jplasma_gradshafranov();

    solve_equation_systems(substep, true);

    if (print_full_output)
        print_output(substep, true);

    // Static part
    program->set_active_mesh(main_mesh);

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