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

#include "analysis_base.h"
#include "apec_m_equation_systems.h"
#include "apec_m_mesh.h"
#include "fem_data_repository.h"
#include "cli.h"

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
                                                             cli(std::make_unique<BasicCommandLineInterface>()),
                                                             model_data_repository(std::make_shared<FEMDataRepository>()) {}

    bool load_data(const std::vector<std::string> &filepaths);

    bool solve_load_step();

    void read_mesh(const std::string &name, const std::string &filepath);

    libMesh::Parallel::Communicator &communicator() { return comm; }

    libMesh::Parallel::Communicator &comm;
    std::shared_ptr<FEMDataRepository> model_data_repository;
    std::unique_ptr<CommandLineInterface> cli;
    std::shared_ptr<AnalysisBase> shared_analysis;
    bool analysis_mode = false;

    bool enter_analysis_mode();
};

#endif // APEC_M_H
