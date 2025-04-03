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

template <typename Container>
void writeElementData(std::ofstream &outFile, const std::vector<Container> &containers)
{
    // Write the number of containers as an unsigned int
    unsigned int numContainers = static_cast<unsigned int>(containers.size());
    outFile.write(reinterpret_cast<const char *>(&numContainers), sizeof(numContainers));

    for (const auto &container : containers)
    {
        const auto &values = container.get_values();

        // Write the size of the values vector as an unsigned int
        unsigned int valuesSize = static_cast<unsigned int>(values.size());
        outFile.write(reinterpret_cast<const char *>(&valuesSize), sizeof(valuesSize));

        // Write the actual values
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

template <typename Container>
void readElementData(std::ifstream &inFile, std::vector<Container> &containers)
{
    // Read the number of containers as an unsigned int
    unsigned int numContainers;
    inFile.read(reinterpret_cast<char *>(&numContainers), sizeof(numContainers));

    // Resize the containers vector to match the number of containers
    containers.resize(numContainers);

    for (auto &container : containers)
    {
        // Read the size of the values vector as an unsigned int
        unsigned int valuesSize;
        inFile.read(reinterpret_cast<char *>(&valuesSize), sizeof(valuesSize));

        // Resize the container's values vector to match the size
        auto &values = container.get_values();
        values.resize(valuesSize);

        // Read the actual values into the container's values vector
        inFile.read(reinterpret_cast<char *>(values.data()), valuesSize * sizeof(Number));
    }
}