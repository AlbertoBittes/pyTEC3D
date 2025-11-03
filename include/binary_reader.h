#ifndef BINARY_READER_H
#define BINARY_READER_H

#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <vector>

#include "binary_labels.h"
#include "fem_data_repository.h"

#include <libmesh/libmesh.h>

class ApecMEquationSystems;
class LoadStep;

class BinaryReader
{
public:
    explicit BinaryReader() = default;

    std::unique_ptr<std::ios> open_file_stream(const std::filesystem::path &folder,
                                               const std::string &filename,
                                               FileMode mode);

    std::vector<std::unique_ptr<std::ios>> open_file_streams(const std::filesystem::path &folder,
                                                             const std::vector<std::string> &filenames,
                                                             FileMode mode);

    void close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams);

    void read_binary_fem_data(const std::filesystem::path &filepath,
                              const std::shared_ptr<FEMDataRepository> &repository = nullptr,
                              const std::shared_ptr<ApecMEquationSystems> &equation_systems = nullptr);

private:
    void dispatch_read(BinaryLabelId label,
                       std::ifstream &file,
                       const std::shared_ptr<FEMDataRepository> &repository,
                       const std::shared_ptr<ApecMEquationSystems> &equation_systems);
};

#endif // BINARY_READER_H
