#ifndef APEC_M_H
#define APEC_M_H

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

#include <libmesh/libmesh.h>

#include "apec_m_eq_system_object.h"
#include "apec_m_mesh.h"
#include "fem_data_repository.h"
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

    Session &session();
    const Session &session() const;

    std::shared_ptr<ApecMMesh> get_mesh(const std::string &name) const;
    std::shared_ptr<ApecMEquationSystems> get_equation_system(const std::string &name) const;

    bool has_mesh(const std::string &name) const;
    bool set_active_mesh(const std::string &name);
    std::optional<ApecMEquationSystems::AnalysisActivationState> activate_analysis(const std::string &mesh_name, const std::string &analysis_name);
    bool load_model_data(const std::string &filepath);

    void read_mesh(const std::string &name, const std::string &filepath);
    void clear_mesh(const std::string &name);

    void set_cli(std::unique_ptr<CommandLineInterface> cli_instance);
    void run_cli();
    void run_batch_file(const std::string &filepath);

    libMesh::Parallel::Communicator &communicator() { return comm; }

private:
    struct MeshContext
    {
        std::shared_ptr<ApecMMesh> mesh;
        std::shared_ptr<ApecMEquationSystems> equation_system;
    };

    MeshContext &ensure_context(const std::string &name);

    libMesh::Parallel::Communicator &comm;
    std::shared_ptr<FEMDataRepository> model_data_repository;
    std::shared_ptr<Session> session_object;
    std::unordered_map<std::string, MeshContext> mesh_registry;
    std::unique_ptr<CommandLineInterface> cli;
};

#endif // APEC_M_H
