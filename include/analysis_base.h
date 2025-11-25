#ifndef ANALYSIS_BASE_H
#define ANALYSIS_BASE_H

#include <cstddef>
#include <memory>
#include <string>

#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "fem_data_repository.h"
#include "fem_service_abstract.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "session.h"
#include "time_step.h"
#include "apec_m_settings.h"

class ApecMEquationSystems;

class AnalysisBase : public FEDataObjectAbstract,
                     public FEParentObjectBase<ApecMEquationSystems>,
                     public FEMServiceWithOldContext<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>
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
                 ApecMEquationSystems *eq_systems = nullptr,
                 std::string analysis_name = "",
                 std::shared_ptr<ContextView> context_view = nullptr);
    virtual ~AnalysisBase() = default;

    const std::string &get_name() const { return name; }
    void set_name(const std::string &analysis_name) { name = analysis_name; }

    Type type() const { return analysis_type; }

    bool initialize_context(std::shared_ptr<ContextView> context_view);
    bool has_valid_context() const { return context_initialized; }
    static std::shared_ptr<ContextView> create_context_view(const std::shared_ptr<FEMDataRepository> &model_repo);

    virtual bool solve(const std::string &identifier) = 0;
    virtual bool solve_transient(const std::string &base_name);
    virtual bool finish();

    std::size_t load_steps_completed() const { return load_steps_completed_count; }

    void update_time_step_context(const std::shared_ptr<TimeStep> &time_step);

    static Type determine_analysis_type(const std::shared_ptr<FEMDataRepository> &model_repo);
    static std::shared_ptr<AnalysisBase> build(ApecMEquationSystems *eq_systems,
                                               const std::shared_ptr<FEMDataRepository> &model_repo,
                                               std::shared_ptr<ContextView> context_view);

protected:
    void mark_load_step_completed();
    void reset_load_step_counter();
    std::size_t load_steps_completed_count = 0;
    bool context_initialized = false;

private:
    virtual bool execute_solve() = 0;
    Type analysis_type;
    std::string name;
};

#endif // ANALYSIS_BASE_H
