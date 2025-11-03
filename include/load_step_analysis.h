#ifndef LOAD_STEP_ANALYSIS_H
#define LOAD_STEP_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "apec_m_analysis.h"
#include "fe_parent_object.h"
#include "fem_data_repository.h"
#include "fem_service_abstract.h"
#include "load_step.h"
#include "session.h"

class ApecMEquationSystems;
class ApecMLinearSystem;
class ApecMExplicitSystem;

class LoadStepAnalysis : public ApecMAnalysis,
                         public FEDataObjectAbstract,
                         public FEParentObjectBase<ApecMEquationSystems>,
                         public FEMServiceWithOldContext<Session, LoadStep>
{
public:
    using Base = FEMServiceWithOldContext<Session, LoadStep>;
    using ContextView = typename Base::ContextView;

    LoadStepAnalysis(ApecMEquationSystems *eq_systems = nullptr,
                     std::shared_ptr<FEMDataRepository> model_data_repo = nullptr);

    void set_model_data(const std::shared_ptr<FEMDataRepository> &model_data_repo);

    std::filesystem::path load_steps_folder() const;
    std::vector<std::filesystem::path> list_load_step_files() const;
    std::size_t load_step_file_count() const;

    std::shared_ptr<LoadStep> read_load_step(const std::filesystem::path &filepath);

    bool solve(const std::string &filename);
    bool solve_transient(const std::string &base_name);

private:
    std::shared_ptr<ApecMEquationSystems> equation_systems_ptr() const;
    void prepare_new_context(const std::shared_ptr<LoadStep> &load_step);
    void configure_services(const std::shared_ptr<LoadStep> &load_step);
    bool execute_solve(const std::shared_ptr<LoadStep> &load_step);

    std::filesystem::path resolve_load_step_path(const std::string &filename) const;
    std::filesystem::path resolve_default_base_name() const;

    std::shared_ptr<FEMDataRepository> model_data;
    unsigned int cumulative_substep_index = 0;
};

#endif // LOAD_STEP_ANALYSIS_H
