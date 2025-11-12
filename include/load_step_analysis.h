#ifndef LOAD_STEP_ANALYSIS_H
#define LOAD_STEP_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "apec_m_analysis.h"

class ApecMLinearSystem;
class ApecMExplicitSystem;

class LoadStepAnalysis : public AnalysisBase
{
public:
    explicit LoadStepAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    std::filesystem::path load_steps_folder();
    std::vector<std::filesystem::path> list_load_step_files();
    std::size_t load_step_file_count();

    std::shared_ptr<TimeStep> read_time_step(const std::filesystem::path &filepath);

    bool solve(const std::string &filename) override;
    bool solve_transient(const std::string &base_name) override;

private:
    void configure_services();
    bool execute_solve();

    std::filesystem::path resolve_load_step_path(const std::string &filename);
    std::filesystem::path resolve_default_base_name();

    unsigned int cumulative_substep_index = 0;
};

#endif // LOAD_STEP_ANALYSIS_H
