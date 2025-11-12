#ifndef ANALYSIS_BASE_H
#define ANALYSIS_BASE_H

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

class ApecMEquationSystems;

class AnalysisBase : public FEDataObjectAbstract,
                     public FEParentObjectBase<ApecMEquationSystems>,
                     public FEMServiceWithOldContext<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session>;

    enum class Type
    {
        Static,
        APECM,
        LoadStep,
        DOFInterpolation
    };

    AnalysisBase(Type type = Type::Static,
                 ApecMEquationSystems *eq_systems = nullptr,
                 std::string analysis_name = "");
    virtual ~AnalysisBase() = default;

    const std::string &get_name() const { return name; }
    void set_name(const std::string &analysis_name) { name = analysis_name; }

    Type type() const { return analysis_type; }

    bool initialize_model_context(const std::shared_ptr<FEMDataRepository> &model_repo);

    virtual bool solve(const std::string &identifier) = 0;
    virtual bool solve_transient(const std::string &base_name);

    void update_time_step_context(const std::shared_ptr<TimeStep> &time_step);

    static Type determine_analysis_type(const std::shared_ptr<FEMDataRepository> &model_repo);
    static std::shared_ptr<AnalysisBase> build(ApecMEquationSystems *eq_systems,
                                               const std::shared_ptr<FEMDataRepository> &model_repo);

private:
    Type analysis_type;
    std::string name;
};

#endif // ANALYSIS_BASE_H
