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

#include "apec_m_equation_systems.h"
#include "apec_m_mesh.h"
#include "fem_data_repository.h"

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
    explicit APEC_M(libMesh::Parallel::Communicator &comm) : comm(comm),
                                                             model_data_repository(std::make_shared<FEMDataRepository>()) {}

    std::shared_ptr<ApecMMesh> get_mesh(const std::string &name) const;
    std::shared_ptr<ApecMMesh> get_active_mesh() const;
    std::shared_ptr<ApecMEquationSystems> get_equation_system(const std::string &name) const;
    std::shared_ptr<ApecMEquationSystems> get_active_equation_system() const;
    std::vector<std::shared_ptr<ApecMEquationSystems>> get_all_equation_systems() const;
    std::vector<std::pair<std::string, std::shared_ptr<ApecMEquationSystems>>> get_all_equation_systems_with_names() const;

    bool has_mesh(const std::string &name) const;
    bool set_active_mesh(const std::string &name);
    bool load_data(const std::vector<std::string> &filepaths);

    bool solve_load_step();
    bool finish_analysis();

    void read_mesh(const std::string &name, const std::string &filepath);
    void clear_mesh(const std::string &name);

    void set_cli(std::unique_ptr<CommandLineInterface> cli_instance);
    void run_cli();

    libMesh::Parallel::Communicator &communicator() { return comm; }

    libMesh::Parallel::Communicator &comm;
    std::shared_ptr<FEMDataRepository> model_data_repository;
    std::vector<std::pair<std::string, std::shared_ptr<FEMDataRepository>>> mesh_repositories;
    std::unordered_map<std::string, std::size_t> mesh_repository_lookup;
    std::unique_ptr<CommandLineInterface> cli;
    std::shared_ptr<AnalysisBase> shared_analysis;
    std::string active_mesh_name = "";
    bool analysis_mode = false;

    FEMDataRepository &ensure_mesh_repository(const std::string &name);
    std::shared_ptr<FEMDataRepository> find_mesh_repository(const std::string &name) const;

    bool enter_analysis_mode();
    void enter_preprocessing_mode();
};

#endif // APEC_M_H
