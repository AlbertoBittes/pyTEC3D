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
#include "session.h"

class ApecMEquationSystems;
class ApecMLinearSystem;
class ApecMExplicitSystem;

class LoadStepAnalysis : public ApecMAnalysis,
                         public FEDataObjectAbstract,
                         public FEParentObjectBase<ApecMEquationSystems>,
                         public FEMServiceWithOldContext<Session, TimeStep>
{
public:
    using Base = FEMServiceWithOldContext<Session, TimeStep>;
    using ContextView = typename Base::ContextView;

    LoadStepAnalysis(ApecMEquationSystems *eq_systems = nullptr,
                     std::shared_ptr<FEMDataRepository> model_data_repo = nullptr);

    void set_model_repository(const std::shared_ptr<FEMDataRepository> &repo);
    void set_session(const std::shared_ptr<Session> &session_ptr);

    std::filesystem::path load_steps_folder() const;
    std::vector<std::filesystem::path> list_load_step_files() const;
    std::size_t load_step_file_count() const;

    std::shared_ptr<TimeStep> read_time_step(const std::filesystem::path &filepath);

    bool solve(const std::string &filename);
    bool solve_transient(const std::string &base_name);

private:
    std::shared_ptr<ApecMEquationSystems> equation_systems_ptr() const;
    void set_active_load_step_repository(const std::shared_ptr<FEMDataRepository> &repo);
    std::shared_ptr<FEMDataRepository> model_repository() const;
    std::shared_ptr<Session> session_ptr() const;
    void configure_services();
    bool execute_solve();

    std::filesystem::path resolve_load_step_path(const std::string &filename) const;
    std::filesystem::path resolve_default_base_name() const;

    std::weak_ptr<FEMDataRepository> model_data_handle;
    std::weak_ptr<Session> session_handle;
    std::shared_ptr<FEMDataRepository> active_load_step_repository;
    unsigned int cumulative_substep_index = 0;
};

#endif // LOAD_STEP_ANALYSIS_H
