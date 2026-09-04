/**
 * @file binary_reader.cpp
 * @brief BinaryReader: rank-local .emat file I/O for element matrix and vector caches.
 */
#include "apec_m_solver/binary_reader.h"

#include <fstream>
#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>
#include "apec_m_solver/apec_m_equation_systems.h"
#include "apec_m_solver/apec_m_system.h"
#include "apec_m_solver/apec_m_logging.h"

BinaryReader::BinaryReader() = default;

BinaryReader::InputStreamPtr BinaryReader::open_input_stream(const std::filesystem::path &filepath)
{
    apec::logging::trace("Opening binary file for read: {}", filepath.string());
    auto stream = std::make_unique<std::ifstream>(filepath, std::ios::binary);
    return stream;
}

BinaryReader::OutputStreamPtr BinaryReader::open_output_stream(const std::filesystem::path &filepath)
{
    const auto folder = filepath.parent_path();
    if (!folder.empty() && !std::filesystem::exists(folder))
        std::filesystem::create_directories(folder);

    apec::logging::trace("Opening binary file for write: {}", filepath.string());
    auto stream = std::make_unique<std::ofstream>(filepath, std::ios::binary | std::ios::trunc);
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

void BinaryReader::write_element_matrices(
    std::ofstream &outFile,
    unsigned int elem_id,
    const std::vector<libMesh::DenseMatrix<libMesh::Number>> &elem_matrices_vector)
{
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

bool BinaryReader::read_element_matrices(
    std::ifstream &inFile,
    unsigned int &elem_id,
    std::vector<libMesh::DenseMatrix<libMesh::Number>> &elem_matrices_vector)
{
    if (!inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id)))
        return false;

    unsigned int num_containers = 0;
    if (!inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers)))
        return false;

    elem_matrices_vector.resize(num_containers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
        unsigned int rows = 0;
        unsigned int cols = 0;
        inFile.read(reinterpret_cast<char *>(&rows), sizeof(rows));
        inFile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

        elem_matrix.resize(rows, cols);
        auto &values = elem_matrix.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }

    return static_cast<bool>(inFile);
}

void BinaryReader::write_element_vectors(
    std::ofstream &outFile,
    unsigned int elem_id,
    const std::vector<libMesh::DenseVector<libMesh::Number>> &elem_vectors_vector)
{
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

bool BinaryReader::read_element_vectors(
    std::ifstream &inFile,
    unsigned int &elem_id,
    std::vector<libMesh::DenseVector<libMesh::Number>> &elem_vectors_vector)
{
    if (!inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id)))
        return false;

    unsigned int num_containers = 0;
    if (!inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers)))
        return false;
    elem_vectors_vector.resize(num_containers);

    for (auto &elem_vector : elem_vectors_vector)
    {
        unsigned int values_size = 0;
        inFile.read(reinterpret_cast<char *>(&values_size), sizeof(values_size));

        elem_vector.resize(values_size);
        auto &values = elem_vector.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(libMesh::Number));
    }

    return static_cast<bool>(inFile);
}
