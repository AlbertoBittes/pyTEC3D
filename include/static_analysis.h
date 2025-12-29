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
    using AnalysisBase::AnalysisBase;

    virtual bool solve() override;
    virtual bool solve_transient() override;
    virtual void prepare_time_substep_solve(Substep &substep) override;
    virtual void solve_equation_systems(Substep &substep, bool nodal_solution = true) override;
    virtual bool run_picard_iterations(Substep &substep) override;
    virtual void print_output(Substep &substep) override;
    void configure_services() const;

protected:
    bool execute_solve() override;
};

class DOFInterpolationAnalysis : public StaticAnalysis
{
public:
    using StaticAnalysis::StaticAnalysis;

    bool solve() override;
};

#endif // STATIC_ANALYSIS_H
