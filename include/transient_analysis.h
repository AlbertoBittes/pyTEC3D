#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "analysis_base.h"

class ApecMLinearSystem;
class ApecMExplicitSystem;

class TransientAnalysis : public AnalysisBase
{
public:
    explicit TransientAnalysis(ApecMEquationSystems *eq_systems = nullptr,
                               std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

    std::filesystem::path load_steps_folder() const;
    std::vector<std::filesystem::path> list_load_step_files() const;
    std::size_t load_step_file_count() const;

    std::shared_ptr<TimeStep> read_time_step(const std::filesystem::path &filepath);

    bool solve(const std::string &filename) override;
    bool solve_transient(const std::string &base_name) override;
    virtual bool run_picard_iterations(Substep &substep);
    bool finish() override;
    void solve_equation_systems(Substep &substep);
    void prepare_time_substep_solve(Substep &substep);
    virtual bool execute_solve() override;
    void print_output(Substep &substep);
    void configure_services() const;

private:
protected:
    std::filesystem::path resolve_load_step_path(const std::string &filename) const;
    std::filesystem::path resolve_default_base_name() const;

private:
    unsigned int cumulative_substep_index = 0;

protected:
    TransientAnalysis(AnalysisBase::Type type,
                      ApecMEquationSystems *eq_systems,
                      std::string analysis_name,
                      std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
};

#endif // TRANSIENT_ANALYSIS_H
