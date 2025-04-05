// Prepare the file element_data_io.h
#ifndef ELEMENT_DATA_IO_H
#define ELEMENT_DATA_IO_H
//
#include <fstream>
#include <vector>
#include <stdexcept>
#include <utility> // for std::pair
#include <cstdint>
#include "../include/element_types.h" // Include the header file for ElementTypesCollector
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

using namespace libMesh;

void buildVectors(std::map<std::string, std::tuple<unsigned int, unsigned int>> &sizes_map, std::vector<DenseVector<Number>> &vectors);

void buildMatrices(std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &sizes_map, std::vector<DenseMatrix<Number>> &matrices);

void writeElementMatrices(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector);

void readElementMatrices(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector);

void writeElementVectors(std::ofstream &outFile, const std::vector<DenseVector<Number>> &elem_vectors_vector);

void readElementVectors(std::ifstream &inFile, std::vector<DenseVector<Number>> &elem_vectors_vector);

void readElementData(std::ifstream &inFile, std::vector<DenseMatrix<Number>> &elem_matrices_vector, std::vector<DenseVector<Number>> &elem_vectors_vector);

void writeElementData(std::ofstream &outFile, const std::vector<DenseMatrix<Number>> &elem_matrices_vector, const std::vector<DenseVector<Number>> &elem_vectors_vector);

// int main() {
//     try {
//         std::string filename = "data.bin";
//         unsigned int numElements;

//         std::cout << "Enter the number of elements: ";
//         std::cin >> numElements;

//         FileHandler fileHandler(filename, numElements);
//         std::ifstream& inFile = fileHandler.getStream();

//         ElementHeader header;
//         std::vector<DenseMatrix<Number>> matrices;
//         std::vector<DenseVector<Number>> vectors;

//         for (unsigned int i = 0; i < numElements; ++i) {
//             // Assume the header is populated manually by the user
//             std::cout << "Enter the number of matrices for element " << i + 1 << ": ";
//             std::cin >> header.numMatrices;
//             header.matrixDimensions.resize(header.numMatrices);
//             header.matrixSizes.resize(header.numMatrices);
//             for (unsigned int j = 0; j < header.numMatrices; ++j) {
//                 std::cout << "Enter rows and columns for matrix " << j + 1 << ": ";
//                 std::cin >> header.matrixDimensions[j].first >> header.matrixDimensions[j].second;
//                 header.matrixSizes[j] = header.matrixDimensions[j].first * header.matrixDimensions[j].second;
//             }

//             std::cout << "Enter the number of vectors for element " << i + 1 << ": ";
//             std::cin >> header.numVectors;
//             header.vectorSizes.resize(header.numVectors);
//             for (unsigned int j = 0; j < header.numVectors; ++j) {
//                 std::cout << "Enter size for vector " << j + 1 << ": ";
//                 std::cin >> header.vectorSizes[j];
//             }

//             buildMatrices(header, matrices);
//             readElementData<DenseMatrix<Number>>(inFile, header.numMatrices, header.matrixSizes, matrices);

//             // Perform operations on matrices
//             // ...

//             // Example: Print the first value of the first matrix
//             if (!matrices.empty() && !matrices[0].get_values().empty()) {
//                 std::cout << "First value of the first matrix: " << matrices[0].get_values()[0] << std::endl;
//             }

//             buildVectors(header, vectors);
//             readElementData<DenseVector<Number>>(inFile, header.numVectors, header.vectorSizes, vectors);

//             // Perform operations on vectors
//             // ...

//             // Example: Print the first value of the first vector
//             if (!vectors.empty() && !vectors[0].get_values().empty()) {
//                 std::cout << "First value of the first vector: " << vectors[0].get_values()[0] << std::endl;
//             }
//         }

//         // Writing data back to a file
//         std::ofstream outFile("output_data.bin", std::ios::binary);
//         if (!outFile) {
//             throw std::runtime_error("Cannot open file for writing");
//         }

//         outFile.write(reinterpret_cast<const char*>(&numElements), sizeof(numElements));

//         for (unsigned int i = 0; i < numElements; ++i) {
//             outFile.write(reinterpret_cast<const char*>(&header.numMatrices), sizeof(header.numMatrices));
//             for (unsigned int j = 0; j < header.numMatrices; ++j) {
//                 outFile.write(reinterpret_cast<const char*>(&header.matrixDimensions[j].first), sizeof(header.matrixDimensions[j].first));
//                 outFile.write(reinterpret_cast<const char*>(&header.matrixDimensions[j].second), sizeof(header.matrixDimensions[j].second));
//             }
//             writeElementData<DenseMatrix<Number>>(outFile, header.numMatrices, header.matrixSizes, matrices);
//             writeElementData<DenseVector<Number>>(outFile, header.numVectors, header.vectorSizes, vectors);
//         }

//         outFile.close();
//     } catch (const std::exception& e) {
//         std::cerr << "Error: " << e.what() << std::endl;
//         return 1;
//     }

//     return 0;
// }
#endif // ELEMENT_DATA_IO_H