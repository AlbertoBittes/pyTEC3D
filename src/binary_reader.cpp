/**
 * @file binary_reader.cpp
 * @brief BinaryReader: rank-local .emat file I/O for element matrix and vector caches.
 */
#include "binary_reader.h"

#include <fstream>
#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "element_builder.h"
#include "element_types_factory.h"
#include "solver_configuration.h"
#include "apec_m_logging.h"

BinaryReader::BinaryReader() = default;

BinaryReader::InputStreamPtr BinaryReader::open_input_stream(const std::filesystem::path &filepath)
{
    apec::logging::trace("Opening binary file for read: {}", filepath.string());
    auto stream = std::make_unique<std::ifstream>(filepath, std::ios::binary);
    if (!stream->is_open())
        apec::logging::error_and_abort("BinaryReader: failed to open file for reading '" + filepath.string() + "'.");
    return stream;
}

BinaryReader::OutputStreamPtr BinaryReader::open_output_stream(const std::filesystem::path &filepath)
{
    const auto folder = filepath.parent_path();
    if (!folder.empty() && !std::filesystem::exists(folder))
        std::filesystem::create_directories(folder);

    apec::logging::trace("Opening binary file for write: {}", filepath.string());
    auto stream = std::make_unique<std::ofstream>(filepath, std::ios::binary | std::ios::trunc);
    if (!stream->is_open())
        apec::logging::error_and_abort("BinaryReader: failed to open file for writing '" + filepath.string() + "'.");
    return stream;
}

std::vector<BinaryReader::InputStreamPtr> BinaryReader::open_input_streams(
    const std::filesystem::path &folder,
    const std::vector<std::string> &filenames)
{
    std::vector<InputStreamPtr> streams;
    streams.reserve(filenames.size());
    for (const auto &name : filenames)
        streams.push_back(open_input_stream(folder / name));
    return streams;
}

std::vector<BinaryReader::OutputStreamPtr> BinaryReader::open_output_streams(
    const std::filesystem::path &folder,
    const std::vector<std::string> &filenames)
{
    std::vector<OutputStreamPtr> streams;
    streams.reserve(filenames.size());
    for (const auto &name : filenames)
        streams.push_back(open_output_stream(folder / name));
    return streams;
}

void BinaryReader::write_element_matrices(ElementBuilder &builder, std::ofstream &outFile, unsigned int matrices_set)
{
    auto &elem_matrices_vector = builder.matrices[matrices_set];

    const unsigned int elem_id = builder.elem_data.elem ? builder.elem_data.elem->id() : 0;
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));

    const unsigned int num_containers = static_cast<unsigned int>(elem_matrices_vector.size());
    outFile.write(reinterpret_cast<const char *>(&num_containers), sizeof(num_containers));

    for (const auto &elem_matrix : elem_matrices_vector)
    {
        const unsigned int rows = elem_matrix.m();
        const unsigned int cols = elem_matrix.n();
        outFile.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
        outFile.write(reinterpret_cast<const char *>(&cols), sizeof(cols));

        const auto &values = elem_matrix.get_values();
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }
}

void BinaryReader::read_element_matrices(ElementBuilder &builder, std::ifstream &inFile, unsigned int matrices_set)
{
    auto &elem_matrices_vector = builder.matrices[matrices_set];

    unsigned int elem_id = 0;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));

    unsigned int num_containers = 0;
    inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers));
    elem_matrices_vector.resize(num_containers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
        unsigned int rows = 0, cols = 0;
        inFile.read(reinterpret_cast<char *>(&rows), sizeof(rows));
        inFile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

        elem_matrix.resize(rows, cols);
        auto &values = elem_matrix.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }
}

void BinaryReader::write_element_vectors(ElementBuilder &builder, std::ofstream &outFile, unsigned int vectors_set)
{
    auto &elem_vectors_vector = builder.vectors[vectors_set];

    const unsigned int elem_id = builder.elem_data.elem ? builder.elem_data.elem->id() : 0;
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));

    const unsigned int num_containers = static_cast<unsigned int>(elem_vectors_vector.size());
    outFile.write(reinterpret_cast<const char *>(&num_containers), sizeof(num_containers));

    for (const auto &elem_vector : elem_vectors_vector)
    {
        const unsigned int values_size = elem_vector.size();
        outFile.write(reinterpret_cast<const char *>(&values_size), sizeof(values_size));

        const auto &values = elem_vector.get_values();
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }
}

void BinaryReader::read_element_vectors(ElementBuilder &builder, std::ifstream &inFile, unsigned int vectors_set)
{
    auto &elem_vectors_vector = builder.vectors[vectors_set];

    unsigned int elem_id = 0;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));

    unsigned int num_containers = 0;
    inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers));
    elem_vectors_vector.resize(num_containers);

    for (auto &elem_vector : elem_vectors_vector)
    {
        unsigned int values_size = 0;
        inFile.read(reinterpret_cast<char *>(&values_size), sizeof(values_size));

        elem_vector.resize(values_size);
        auto &values = elem_vector.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }
}
