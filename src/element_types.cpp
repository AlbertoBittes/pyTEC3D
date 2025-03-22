#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "../include/element_types.h"

std::shared_ptr<ElementTypesMap> readElementTypes(const std::string &filename)
{
    auto elementTypes = std::make_shared<ElementTypesMap>();
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return elementTypes;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        ElementTypeProperties elementType;
        int element_number;
        char comma; // to consume the commas

        ss >> element_number >> comma;

        // Read ansys_element_type
        std::getline(ss, elementType.ansys_element_type, ',');

        // Read keyopts
        std::string value;
        while (std::getline(ss, value, ','))
        {
            if (!value.empty())
            {
                elementType.keyopts.push_back(std::stod(value));
            }
            else
            {
                elementType.keyopts.push_back(std::numeric_limits<double>::quiet_NaN()); // Default value for empty keyopts
            }
        }

        (*elementTypes)[element_number] = elementType;
    }

    file.close();
    return elementTypes;
}