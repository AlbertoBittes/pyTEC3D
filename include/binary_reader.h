#ifndef BINARY_READER_H
#define BINARY_READER_H

#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <vector>

#include "binary_labels.h"
#include "fem_data_repository.h"
#include "session.h"

#include <libmesh/libmesh.h>
#include <libmesh/petsc_linear_solver.h>

class ApecMEquationSystems;
class LoadStep;

class BinaryReader
{
public:
    explicit BinaryReader(Session &session);

    std::unique_ptr<std::ios> open_file_stream(const std::filesystem::path &folder,
                                               const std::string &filename,
                                               FileMode mode);

    std::vector<std::unique_ptr<std::ios>> open_file_streams(const std::filesystem::path &folder,
                                                             const std::vector<std::string> &filenames,
                                                             FileMode mode);

    void close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams);

    void read_binary_fem_data(const std::filesystem::path &filepath,
                              FEMDataRepository &repository,
                              const std::shared_ptr<ApecMEquationSystems> &equation_systems = nullptr,
                              libMesh::PetscLinearSolver<Number> *solver = nullptr,
                              const std::shared_ptr<LoadStep> &previous_load_step = nullptr);

private:
    Session &session;

    void dispatch_read(BinaryLabelId label,
                       std::ifstream &file,
                       FEMDataRepository &repository,
                       const std::shared_ptr<ApecMEquationSystems> &equation_systems,
                       libMesh::PetscLinearSolver<Number> *solver,
                       const std::shared_ptr<LoadStep> &previous_load_step);
};

#endif // BINARY_READER_H
