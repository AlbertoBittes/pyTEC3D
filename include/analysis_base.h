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
    virtual void prepare_time_step_solve(TimeStep &time_step);
    virtual void solve_equation_systems(TimeStep &time_step, bool nodal_solution);
    virtual bool run_picard_iterations(TimeStep &time_step);
    virtual void print_output(TimeStep &time_step, bool picard_iteration_print = false);

    virtual bool finish();

    static std::shared_ptr<AnalysisBase> build(APEC_M *program,
                                               std::shared_ptr<ContextView> context_view);

    AnalysisBase(APEC_M *program, std::shared_ptr<ContextView> context_view);

protected:
    virtual bool execute_solve() = 0;

    bool context_initialized = false;
};

#endif // ANALYSIS_BASE_H
