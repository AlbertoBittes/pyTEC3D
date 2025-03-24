#include "../include/element_data_io.h"

void buildVectors(const ElementHeader &header, std::vector<DenseVector<Number>> &vectors)
{
    vectors.resize(header.numVectors);
    for (unsigned int i = 0; i < header.numVectors; ++i)
    {
        vectors[i] = DenseVector<Number>(header.vectorSizes[i]);
    }
}

void buildMatrices(const ElementHeader &header, std::vector<DenseMatrix<Number>> &matrices)
{
    matrices.resize(header.numMatrices);
    for (unsigned int i = 0; i < header.numMatrices; ++i)
    {
        matrices[i] = DenseMatrix<Number>(header.matrixDimensions[i].first, header.matrixDimensions[i].second);
    }
}

template <typename Container>
void writeElementData(std::ofstream &outFile, unsigned int numItems, const std::vector<unsigned int> &sizes, const std::vector<Container> &containers)
{
    outFile.write(reinterpret_cast<const char *>(&numItems), sizeof(numItems));

    for (const auto &size : sizes)
    {
        outFile.write(reinterpret_cast<const char *>(&size), sizeof(size));
    }

    for (const auto &container : containers)
    {
        const auto &values = container.get_values();
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

template <typename Container>
void readElementData(std::ifstream &inFile, unsigned int numItems, std::vector<unsigned int> &sizes, std::vector<Container> &containers)
{
    if (containers.size() != numItems)
    {
        throw std::runtime_error("Size of containers array does not match the number of items in the header");
    }

    sizes.resize(numItems);

    for (auto &size : sizes)
    {
        inFile.read(reinterpret_cast<char *>(&size), sizeof(size));
    }

    for (unsigned int i = 0; i < numItems; ++i)
    {
        inFile.read(reinterpret_cast<char *>(containers[i].get_values().data()), sizes[i] * sizeof(Number));
    }
}