#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "../include/element_types.h"

ElementTypesMapCSV readElementTypes(const std::string &filename)
{
    ElementTypesMapCSV elementTypesMap;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string ansys_element_type_str, element_type;
        std::vector<int> keyopts(9, 0); // Initialize with 9 zeros
        unsigned int ansys_element_type;

        if (!(std::getline(iss, ansys_element_type_str, ',') && std::getline(iss, element_type, ',')))
        {
            throw std::runtime_error("Error parsing line: " + line);
        }

        ansys_element_type = std::stoul(ansys_element_type_str);

        for (int i = 0; i < 9; ++i)
        {
            std::string keyopt_str;
            if (std::getline(iss, keyopt_str, ',') && !keyopt_str.empty())
            {
                keyopts[i] = std::stoi(keyopt_str);
            }
        }

        ElementTypeProperties properties = {element_type, keyopts};
        (*elementTypesMap)[ansys_element_type] = properties;
    }

    file.close();
    return elementTypesMap;
}