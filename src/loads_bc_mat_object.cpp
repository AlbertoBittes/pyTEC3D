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
#include "../include/loads_bc_mat_object.h"
#include "libmesh/linear_implicit_system.h"

using namespace libMesh;

CurrentLoadsMap LoadsBCMatObject::buildCurrentLoadsMapBFElist(const std::string &filename)
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
                std::vector<Real> values;
                std::istringstream dataStream(line);
                Real jsX, jsY, jsZ;

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

NodeDataMatrix LoadsBCMatObject::readInputDlist(const std::string &filename)
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
            Real node, real;
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
NodeBCMap LoadsBCMatObject::buildNodeBCMapDlist(const std::string &dlist_filepath)
{
    NodeBCMap table;
    NodeDataMatrix bc_matrix = this->readInputDlist(dlist_filepath);
    for (const auto &entry : bc_matrix)
    {
        dof_id_type node = std::get<0>(entry);
        std::string label = std::get<1>(entry);
        Real value = std::get<2>(entry);

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

LoadsBCMatObject::LoadsBCMatObject()
{
    // Constructor implementation
}

void LoadsBCMatObject::loadMPLIST(const std::string &filePath)
{
    std::ifstream file(filePath);
    std::string line;
    unsigned short int materialNumber = 0;
    Real RSVX = 0.0, MURX = 0.0;

    if (file.is_open())
    {
        while (getline(file, line))
        {
            std::istringstream iss(line);
            std::string key;
            iss >> key;

            if (key == "MATERIAL")
            {
                iss >> key; // Skip "NUMBER"
                iss >> materialNumber;
            }
            else if (key == "TEMP")
            {
                iss >> key; // Read the next key to determine if it's RSVX or MURX
                if (key == "RSVX")
                {
                    iss >> RSVX;
                }
                else if (key == "MURX")
                {
                    iss >> MURX;
                    // Use the setter to add the material properties to the map
                    setMaterialProperty(materialNumber, RSVX, MURX);
                }
            }
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << filePath << std::endl;
    }
}

void LoadsBCMatObject::loadMaterialsFromCSV(const std::string &filePath)
{
    std::ifstream file(filePath);
    std::string line;
    // Skip the header line if present
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;
        unsigned short int material_id;
        Real rsvx, murx;

        // Read material_id
        std::getline(iss, token, ',');
        material_id = static_cast<unsigned short int>(std::stoi(token));

        // Read rsvx
        std::getline(iss, token, ',');
        rsvx = static_cast<Real>(std::stod(token));

        // Read murx
        std::getline(iss, token, ',');
        murx = static_cast<Real>(std::stod(token));

        // Use the existing method to set the material property
        setMaterialProperty(material_id, rsvx, murx);
    }

    file.close();
}