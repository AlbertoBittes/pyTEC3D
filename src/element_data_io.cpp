#include "../include/element_data_io.h"

void buildVectors(std::map<std::string, std::tuple<unsigned int, unsigned int>> &sizes_map, std::vector<DenseVector<Number>> &vectors)
{
    vectors.resize(sizes_map.size());
    for (const auto &[key, value] : sizes_map)
    {
        unsigned int index = std::distance(sizes_map.begin(), sizes_map.find(key));
        vectors[index] = DenseVector<Number>(std::get<1>(value));
    }
}

void buildMatrices(std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &sizes_map, std::vector<DenseMatrix<Number>> &matrices)
{
    matrices.resize(sizes_map.size());
    for (const auto &[key, value] : sizes_map)
    {
        unsigned int index = std::distance(sizes_map.begin(), sizes_map.find(key));
        matrices[index] = DenseMatrix<Number>(std::get<1>(value), std::get<2>(value));
    }
}

void writeElementMatrices(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector)
{
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(elem_matrices_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_matrix : elem_matrices_vector)
    {
        const auto &values = elem_matrix.get_values();

        // Write the size of the values vector as an unsigned int
        // unsigned int valuesSize = static_cast<unsigned int>(values.size());
        unsigned int rows = elem_matrix.m();
        unsigned int cols = elem_matrix.n();

        // Write the number of rows and columns as unsigned ints
        outFile.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
        outFile.write(reinterpret_cast<const char *>(&cols), sizeof(cols));

        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void readElementMatrices(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector)
{
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_matrices_vector.resize(numContainers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
        // Read the size of the values vector as an unsigned int
        unsigned int rows, cols;
        inFile.read(reinterpret_cast<char *>(&rows), sizeof(rows));
        inFile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

        // Resize the container's values vector to match the size
        elem_matrix.resize(rows, cols);
        auto &values = elem_matrix.get_values();

        // Read the actual values into the container's values vector
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}

void readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Read the matrices and vectors
    readElementMatrices(inFile, elem_matrices_vector);
    readElementVectors(inFile, elem_vectors_vector);
}

void writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Write the matrices and vectors
    writeElementMatrices(outFile, elem_matrices_vector);
    writeElementVectors(outFile, elem_vectors_vector);
}

void writeElementVectors(std::ofstream &outFile, const std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(elem_vectors_vector.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &elem_vector : elem_vectors_vector)
    {
        const auto &values = elem_vector.get_values();

        // Write the size of the values vector as an unsigned int
        unsigned int valuesSize = elem_vector.size();
        outFile.write(reinterpret_cast<const char *>(&valuesSize), sizeof(valuesSize));

        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void readElementVectors(std::ifstream &inFile, std::vector<DenseVector<Number>> &elem_vectors_vector)
{
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    elem_vectors_vector.resize(numContainers);

    for (auto &elem_vector : elem_vectors_vector)
    {
        // Read the size of the values vector as an unsigned int
        unsigned int valuesSize;
        inFile.read(reinterpret_cast<char *>(&valuesSize), sizeof(valuesSize));

        // Resize the container's values vector to match the size
        elem_vector.resize(valuesSize);
        auto &values = elem_vector.get_values();

        // Read the actual values into the container's values vector
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}