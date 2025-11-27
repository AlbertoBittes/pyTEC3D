#ifndef APEC_M_H
#define APEC_M_H

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <optional>

#include <libmesh/libmesh.h>

#include "apec_m_eq_system_object.h"
#include "apec_m_mesh.h"
#include "fem_data_repository.h"
#include "transient_analysis.h"
#include "session.h"

class CommandLineInterface
{
public:
    virtual ~CommandLineInterface() = default;
    virtual void run_interactive(class APEC_M &program) = 0;
    virtual void run_batch(class APEC_M &program, const std::string &batch_file) = 0;
};

class APEC_M
{
public:
    explicit APEC_M(libMesh::Parallel::Communicator &comm);

    std::shared_ptr<FEMDataRepository> model_data() const { return model_data_repository; }
    std::shared_ptr<FEMDataRepository> model_data() { return model_data_repository; }

    std::shared_ptr<ApecMMesh> get_mesh(const std::string &name) const;
    std::shared_ptr<ApecMMesh> get_active_mesh() const;
    std::shared_ptr<ApecMEquationSystems> get_equation_system(const std::string &name) const;
    std::shared_ptr<ApecMEquationSystems> get_active_equation_system() const;
    std::shared_ptr<AnalysisBase> get_analysis(const std::string &name) const;
    std::shared_ptr<AnalysisBase> get_analysis() const;

    bool has_mesh(const std::string &name) const;
    bool set_active_mesh(const std::string &name);
    bool load_model_data(const std::string &filepath);
    bool load_data(const std::vector<std::string> &filepaths);
    bool set_mode(Session::Mode mode, std::optional<AnalysisBase::Type> analysis_type = std::nullopt);
    Session::Mode get_mode() const
    {
        auto session = model_data_repository->get_object<Session>();
        return session ? session->get_mode() : Session::Mode::PreProcessing;
    }

    bool distribute_boundary_ids();
    bool generate_infinite_elements();
    bool generate_nodal_couplings();
    bool generate_periodic_boundaries();

    bool solve_load_step();
    bool solve_transient_load_steps();
    bool finish_analysis();

    void print_model_data() const;

    void read_mesh(const std::string &name, const std::string &filepath);
    void clear_mesh(const std::string &name);

    void set_cli(std::unique_ptr<CommandLineInterface> cli_instance);
    void run_cli();
    void run_batch_file(const std::string &filepath);

    libMesh::Parallel::Communicator &communicator() { return comm; }

private:
    FEMDataRepository &ensure_mesh_repository(const std::string &name);
    std::shared_ptr<FEMDataRepository> find_mesh_repository(const std::string &name) const;
    bool generate_nodal_couplings_for_mesh(const std::string &name,
                                           const std::shared_ptr<ApecMMesh> &mesh);

    libMesh::Parallel::Communicator &comm;
    std::shared_ptr<FEMDataRepository> model_data_repository;
    std::vector<std::pair<std::string, std::shared_ptr<FEMDataRepository>>> mesh_repositories;
    std::unordered_map<std::string, std::size_t> mesh_repository_lookup;
    std::unique_ptr<CommandLineInterface> cli;
    std::shared_ptr<AnalysisBase> shared_analysis;

    bool enter_analysis_mode(std::optional<AnalysisBase::Type> analysis_type);
    void enter_preprocessing_mode();
};

#endif // APEC_M_H
