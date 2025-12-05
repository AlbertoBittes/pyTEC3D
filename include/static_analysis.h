#ifndef STATIC_ANALYSIS_H
#define STATIC_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <utility>

#include "analysis_base.h"

class StaticAnalysis : public AnalysisBase
{
public:
    explicit StaticAnalysis(APEC_M *program = nullptr,
                            std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr,
                            Type type = Type::STATIC,
                            std::string analysis_name = "StaticAnalysis");

    virtual bool solve() override;
    virtual bool solve_transient() override;
    virtual void prepare_time_substep_solve(Substep &substep) override;
    virtual void solve_equation_systems(Substep &substep) override;
    virtual bool run_picard_iterations(Substep &substep) override;
    virtual void print_output(Substep &substep) override;
    void configure_services() const;

protected:
    StaticAnalysis(Type type,
                   APEC_M *program,
                   std::shared_ptr<AnalysisBase::ContextView> context_view,
                   std::string analysis_name);
    bool execute_solve() override;

    std::filesystem::path resolve_grid_step_path(const std::string &filename) const;
    std::filesystem::path resolve_time_step_path(const std::string &filename) const;
    std::filesystem::path resolve_default_grid_base_name() const;
    std::filesystem::path resolve_default_time_base_name() const;
};

class DOFInterpolationAnalysis : public StaticAnalysis
{
public:
    explicit DOFInterpolationAnalysis(APEC_M *program,
                                      std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

    bool solve() override;
};

inline DOFInterpolationAnalysis::DOFInterpolationAnalysis(APEC_M *program,
                                                          std::shared_ptr<AnalysisBase::ContextView> context_view)
    : StaticAnalysis(Type::DOF_INTERPOLATION, program, std::move(context_view), "DOFInterpolationAnalysis")
{
}

#endif // STATIC_ANALYSIS_H
