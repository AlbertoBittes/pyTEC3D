/**
 * @file binary_reader.h
 * @brief Binary I/O helpers for per-element matrix and vector caches.
 *
 * Element matrices (stiffness Ke, capacitance Ce) and right-hand-side vectors
 * are expensive to recompute every time step.  BinaryReader serialises them to
 * rank-local .emat files after the first assembly pass and deserialises them on
 * subsequent time steps when no recomputation is needed.
 */
#ifndef BINARY_READER_H
#define BINARY_READER_H

#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <stdexcept>
#include <vector>
#include <unordered_set>

#include <libmesh/libmesh.h>

class ElementBuilder;
class ApecMEquationSystems;

/** @brief Unused file-mode tag reserved for future streaming extensions. */
enum class FileMode
{
    Read,
    Write
};

/**
 * @brief Serialises and deserialises per-element matrix and vector caches.
 *
 * Opens rank-local binary files under the current working directory's
 * @c emat/ sub-folder and streams libMesh DenseMatrix / DenseVector data
 * to or from them, allowing the assembler to skip recomputation when the
 * element matrices have not changed.
 */
class BinaryReader
{
public:
    using InputStreamPtr = std::unique_ptr<std::ifstream>;
    using OutputStreamPtr = std::unique_ptr<std::ofstream>;

    BinaryReader();

    /** @brief Opens a single binary file for reading. */
    InputStreamPtr open_input_stream(const std::filesystem::path &filepath);

    /** @brief Opens a single binary file for writing, creating parent directories as needed. */
    OutputStreamPtr open_output_stream(const std::filesystem::path &filepath);

    /** @brief Opens multiple binary files for reading from @p folder. */
    std::vector<InputStreamPtr> open_input_streams(const std::filesystem::path &folder,
                                                   const std::vector<std::string> &filenames);

    /** @brief Opens multiple binary files for writing inside @p folder. */
    std::vector<OutputStreamPtr> open_output_streams(const std::filesystem::path &folder,
                                                     const std::vector<std::string> &filenames);

    /** @brief Closes and destroys a vector of open file streams. */
    template <typename StreamT>
    static void close_file_streams(std::vector<std::unique_ptr<StreamT>> &streams)
    {
        for (auto &stream : streams)
            stream.reset();
        streams.clear();
    }

    /** @brief Serialises one element's matrices for @p matrices_set to a binary stream. */
    static void write_element_matrices(ElementBuilder &builder, std::ofstream &outFile, unsigned int matrices_set);

    /** @brief Deserialises one element's matrices for @p matrices_set from a binary stream. */
    static void read_element_matrices(ElementBuilder &builder, std::ifstream &inFile, unsigned int matrices_set);

    /** @brief Serialises one element's vectors for @p vectors_set to a binary stream. */
    static void write_element_vectors(ElementBuilder &builder, std::ofstream &outFile, unsigned int vectors_set);

    /** @brief Deserialises one element's vectors for @p vectors_set from a binary stream. */
    static void read_element_vectors(ElementBuilder &builder, std::ifstream &inFile, unsigned int vectors_set);
};

#endif // BINARY_READER_H
