#include "apec_m_evolution_analysis.h"

#include <filesystem>

#include <libmesh/libmesh.h>

APECMEvolutionAnalysis::APECMEvolutionAnalysis(ApecMEquationSystems *eq_systems,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : APECMEquilibriumAnalysis(eq_systems, std::move(context_view))
{
    set_name("APECMEvolutionAnalysis");
}

bool APECMEvolutionAnalysis::solve_transient(const std::string &base_name)
{
    return this->solve(base_name);
}

bool APECMEvolutionAnalysis::solve(const std::string &filename)
{
    std::filesystem::path path;
    if (filename.empty())
    {
        auto base = resolve_default_base_name();
        if (base.empty())
        {
            libMesh::out << "Unable to determine default load step base name for evolution analysis." << std::endl;
            return false;
        }

        const auto next_index = load_steps_completed() + 1;
        path = base;
        path.replace_filename(base.stem().string() + "_" + std::to_string(next_index));
        path.replace_extension(".bin");
    }
    else
        path = resolve_load_step_path(filename);

    if (path.empty())
    {
        libMesh::out << "Unable to resolve load step path for evolution analysis." << std::endl;
        return false;
    }

    if (!std::filesystem::exists(path))
    {
        libMesh::out << "Evolution load step '" << path << "' not found." << std::endl;
        return false;
    }

    if (!read_time_step(path))
        return false;

    auto time_step = this->context<TimeStep>();
    if (!time_step)
    {
        libMesh::out << "Evolution analysis requires load step data." << std::endl;
        return false;
    }

    set_plasma_state_from_profile(0);
    if (!solve_equilibrium_step(*time_step))
        return false;
    mark_load_step_completed();

    if (plasma_current_profile().size() <= 1)
        return true;

    const Real base_dt = time_step->get_delta_t();
    const Real base_theta = time_step->get_theta();
    const int substeps = static_cast<int>(time_step->get_n_substeps() ? time_step->get_n_substeps() : 1);
    const int write_results = time_step->get_write_results();
    const int ramped = time_step->get_ramped_load();

    Real current_time = time_step->get_time();
    int step_number = time_step->get_step();

    for (std::size_t idx = 1; idx < plasma_current_profile().size(); ++idx)
    {
        auto new_step = std::make_shared<TimeStep>();
        new_step->set_step(++step_number);
        new_step->set_time_beginning(current_time);
        new_step->set_time(current_time + base_dt);
        new_step->set_delta_t(base_dt);
        new_step->set_theta(base_theta);
        new_step->set_write_results(write_results);
        new_step->set_ramped_load(ramped);
        new_step->set_transient_effect(1);
        new_step->set_substeps(substeps);

        update_time_step_context(new_step);
        current_time = new_step->get_time();
        set_plasma_state_from_profile(idx);

        if (!solve_evolution_step(*new_step, idx))
            return false;
        mark_load_step_completed();
    }

    return true;
}

bool APECMEvolutionAnalysis::solve_evolution_step(TimeStep &time_step, std::size_t profile_index)
{
    prepare_time_step(time_step, false);
    libMesh::out << "APECM evolution step " << profile_index + 1 << " using I_p = "
                 << plasma_current_Ip << " (beta_p = " << beta << ")." << std::endl;
    return run_picard_iterations(time_step);
}
