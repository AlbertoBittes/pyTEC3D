#ifndef ANALYSIS_BASE_H
#define ANALYSIS_BASE_H

#include <cstddef>
#include <functional>
#include <memory>
#include <string>

#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "fem_data_repository.h"
#include "fem_service.h"
#include "element_types_factory.h"
#include "time_step.h"

class ApecMEquationSystems;
class APEC_M;

class AnalysisBase : public FEDataObjectAbstract,
                     public FEParentObjectBase<APEC_M>,
                     public FEMService<ElementTypesFactory, TimeStep, Materials, Constants>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory, TimeStep, Materials, Constants>;

    virtual ~AnalysisBase() = default;

    // bool initialize_context(std::shared_ptr<ContextView> context_view);
    bool has_valid_context() const { return context_initialized; }

    virtual bool solve() = 0;
    virtual void prepare_time_substep_solve(Substep &substep) = 0;
    virtual void solve_equation_systems(Substep &substep, bool nodal_solution) = 0;
    virtual bool run_picard_iterations(Substep &substep) = 0;
    virtual void print_output(Substep &substep, bool picard_iteration_print = false);

    // public:
    //     using AnalysisBase::AnalysisBase;

    //     virtual bool solve() override;
    //     virtual bool solve_transient() override;
    //     virtual void prepare_time_substep_solve(Substep &substep) override;
    //     virtual void solve_equation_systems(Substep &substep, bool nodal_solution = true) override;
    //     virtual bool run_picard_iterations(Substep &substep) override;
    //     virtual void print_output(Substep &substep, bool picard_iteration_print = false) override;
    //     void configure_services() const;

    // protected:
    //     bool execute_solve() override;

    virtual bool finish();

    // void update_time_data_context(const std::shared_ptr<TimeStep> &time_step,
    //                               const std::shared_ptr<ApecMSettingsObject> &settings);

    static std::shared_ptr<AnalysisBase> build(APEC_M *program,
                                               std::shared_ptr<ContextView> context_view);

    AnalysisBase(APEC_M *program, std::shared_ptr<ContextView> context_view);

protected:
    virtual bool execute_solve() = 0;

    bool context_initialized = false;
};

#endif // ANALYSIS_BASE_H
