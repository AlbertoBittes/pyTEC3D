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

class ElementBuilder;
class ApecMEquationSystems;

enum class FileMode
{
    Read,
    Write
};
class BinaryReader
{
public:
    explicit BinaryReader() = default;

    std::unique_ptr<std::ios> open_file_stream(const std::filesystem::path &filepath,
                                               FileMode mode);

    std::vector<std::unique_ptr<std::ios>> open_file_streams(const std::filesystem::path &folder,
                                                             const std::vector<std::string> &filenames,
                                                             FileMode mode);

    void close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams);

    void read_binary_fem_data(const std::string &filepath,
                              const std::shared_ptr<FEMDataRepository> &repository = nullptr,
                              const std::shared_ptr<ApecMEquationSystems> &equation_systems = nullptr);

    static void write_element_matrices(ElementBuilder &builder, std::ofstream &outFile, unsigned int matrices_set);
    static void read_element_matrices(ElementBuilder &builder, std::ifstream &inFile, unsigned int matrices_set);
    static void write_element_vectors(ElementBuilder &builder, std::ofstream &outFile, unsigned int vectors_set);
    static void read_element_vectors(ElementBuilder &builder, std::ifstream &inFile, unsigned int vectors_set);

private:
    void read_binary_fem_data_from_stream(std::ifstream &file,
                                          const std::shared_ptr<FEMDataRepository> &repository,
                                          const std::shared_ptr<ApecMEquationSystems> &equation_systems);
    void dispatch_read(BinaryLabelId label,
                       std::ifstream &file,
                       const std::shared_ptr<FEMDataRepository> &repository,
                       const std::shared_ptr<ApecMEquationSystems> &equation_systems);
};

#endif // BINARY_READER_H
