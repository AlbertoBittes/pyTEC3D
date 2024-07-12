#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include "libmesh/libmesh.h"
#include "../include/loads_bc_object.h"
#include "libmesh/linear_implicit_system.h"

using namespace libMesh;

template <typename Sys>
CurrentLoadsMap LoadsBCObject<Sys>::buildCurrentLoadsMapBFElist(const std::string &filename)
{
    std::ifstream file(filename);
    CurrentLoadsMap map;
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return map; // Return empty map on failure
    }

    while (getline(file, line))
    {
        if (line.find("ELEMENT=") != std::string::npos)
        {
            std::istringstream iss(line);
            std::string temp;
            dof_id_type elementNumber; // Changed type to dof_id_type

            // Extract element number
            iss >> temp >> elementNumber; // temp will hold "ELEMENT=" string

            // Assuming the next line contains the data
            if (getline(file, line))
            {
                std::vector<double> values;
                std::istringstream dataStream(line);
                double jsX, jsY, jsZ;

                // Extract the first three values
                dataStream >> jsX >> jsY >> jsZ;
                values.push_back(jsX);
                values.push_back(jsY);
                values.push_back(jsZ);

                // Use elementNumber as key to insert values into the map
                map[elementNumber] = values;
            }
        }
    }

    return map;
}

template <typename Sys>
NodeDataMatrix LoadsBCObject<Sys>::readInputDlist(const std::string &filename)
{
    std::ifstream file(filename);
    NodeDataMatrix matrix;
    std::string line;
    std::regex linePattern(R"(^\s*\d+)"); // Regex to match lines starting with whitespace(s) followed by a number

    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return matrix;
    }

    while (getline(file, line))
    {
        if (std::regex_search(line, linePattern))
        {
            std::istringstream iss(line);
            double node, real;
            std::string label;

            if (iss >> node >> label >> real)
            {
                // Cast node to int and retain the label as a string
                matrix.push_back(std::make_tuple(static_cast<dof_id_type>(node), label, real));
            }
        }
    }

    return matrix;
}

// Function to process the vector of tuples and build the table
template <typename Sys>
NodeBCMap LoadsBCObject<Sys>::buildNodeBCMapDlist(const std::string &dlist_filepath)
{
    NodeBCMap table;
    NodeDataMatrix bc_matrix = this->readInputDlist(dlist_filepath);
    for (const auto &entry : bc_matrix)
    {
        dof_id_type node = std::get<0>(entry);
        std::string label = std::get<1>(entry);
        double value = std::get<2>(entry);

        // Initialize the node in the map if it doesn't exist
        if (table.find(node) == table.end())
        {
            table[node] = NodeData();
        }

        // Update the corresponding value based on the label
        if (label == "AX")
        {
            table[node].AX = value;
        }
        else if (label == "AY")
        {
            table[node].AY = value;
        }
        else if (label == "AZ")
        {
            table[node].AZ = value;
        }
        else if (label == "VOLT")
        {
            table[node].VOLT = value;
        }
    }

    return table;
}

void processInput(LinearImplicitSystem &system, std::string bc_filename, std::string loads_filename)
{

    // Print the matrix for verification
    // for (const auto &pair : loads_map)
    // {
    //     std::cout << pair.first << ": ";    // Print the key
    //     for (const auto &val : pair.second) // Iterate over the vector (or container) of values
    //     {
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Assuming readInput is defined and reads data into a NodeDataMatrix

    // Example: Print the table
    // for (const auto &entry : bc_map)
    // {
    //     std::cout << "Node " << entry.first << ": "
    //               << "AX=" << entry.second.AX << ", "
    //               << "AY=" << entry.second.AY << ", "
    //               << "AZ=" << entry.second.AZ << ", "
    //               << "VOLT=" << entry.second.VOLT << std::endl;
    // }
}