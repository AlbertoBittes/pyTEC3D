#ifndef ANALYSIS_BASE_H
#define ANALYSIS_BASE_H

#include <cstddef>
#include <memory>
#include <string>

#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "fem_data_repository.h"
#include "fem_service.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "session.h"
#include "time_step.h"
#include "apec_m_settings.h"

class ApecMEquationSystems;
class APEC_M;

class AnalysisBase : public FEDataObjectAbstract,
                     public FEParentObjectBase<APEC_M>,
                     public FEMService<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>;

    enum class Type
    {
        STATIC,
        TRANSIENT,
        APECM_EQUILIBRIUM,
        APECM_EVOLUTION,
        DOF_INTERPOLATION
    };

    AnalysisBase(Type type = Type::STATIC,
                 APEC_M *program = nullptr,
                 std::shared_ptr<ContextView> context_view = nullptr);
    virtual ~AnalysisBase() = default;

    Type type() const { return analysis_type; }

    // bool initialize_context(std::shared_ptr<ContextView> context_view);
    bool has_valid_context() const { return context_initialized; }
    // static std::shared_ptr<ContextView> create_context_view(const std::shared_ptr<FEMDataRepository> &model_repo);
    ApecMEquationSystems *active_equation_system() const;
    std::shared_ptr<ApecMEquationSystems> active_equation_system_ptr() const;

    virtual bool solve() = 0;
    virtual bool solve_transient() = 0;
    virtual void prepare_time_substep_solve(Substep &substep) = 0;
    virtual void solve_equation_systems(Substep &substep) = 0;
    virtual bool run_picard_iterations(Substep &substep) = 0;
    virtual void print_output(Substep &substep) = 0;

    virtual bool finish();

    // void update_time_data_context(const std::shared_ptr<TimeStep> &time_step,
    //                               const std::shared_ptr<ApecMSettingsObject> &settings);

    static std::shared_ptr<AnalysisBase> build(APEC_M *program,
                                               std::shared_ptr<ContextView> context_view,
                                               Type requested_type);

protected:
    virtual bool execute_solve() = 0;

    bool context_initialized = false;
    Type analysis_type;
};

#endif // ANALYSIS_BASE_H
