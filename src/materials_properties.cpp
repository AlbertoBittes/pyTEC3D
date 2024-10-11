#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <limits> // for std::numeric_limits
#include "materials_properties.h"

std::unordered_map<int, MaterialsProperties> readMaterialsProperties(const std::string &filename)
{
    std::unordered_map<int, MaterialsProperties> materials;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return materials;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        MaterialsProperties material;
        int material_number;
        char comma; // to consume the commas

        std::string value;

        ss >> material_number >> comma;

        // Read murx
        std::getline(ss, value, ',');
        material.murx = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read mury
        std::getline(ss, value, ',');
        material.mury = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read murz
        std::getline(ss, value, ',');
        material.murz = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read rsvx
        std::getline(ss, value, ',');
        material.rsvx = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read rsvy
        std::getline(ss, value, ',');
        material.rsvy = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read rsvz
        std::getline(ss, value, ',');
        material.rsvz = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read perx
        std::getline(ss, value, ',');
        material.perx = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read pery
        std::getline(ss, value, ',');
        material.pery = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        // Read perz
        std::getline(ss, value, ',');
        material.perz = value.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(value);

        materials[material_number] = material;
    }

    file.close();
    return materials;
}