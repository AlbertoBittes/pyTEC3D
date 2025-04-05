#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "../include/element_types.h"

ElementTypesMapCSV &ElementTypesCollector::readElementTypes(const std::string &filename)
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
        elementTypesMap[ansys_element_type] = properties;
    }

    file.close();
    return elementTypesMap;
}

void ElementTypesCollector::collectElementTypes(ElementTypesMapCSV &elementTypesMapCSV)
{
    for (const auto &pair : elementTypesMapCSV)
    {
        const auto &properties = pair.second;

        // Create the appropriate ElementType object
        if (properties.element_type == "SOLID97")
        {
            elementTypes[pair.first] = std::make_shared<SOLID97Type>(SOLID97Type(properties));
        }
        else if (properties.element_type == "SOLID236")
        {
            elementTypes[pair.first] = std::make_shared<SOLID236Type>(SOLID236Type(properties));
        }
        else if (properties.element_type == "PLANE13")
        {
            elementTypes[pair.first] = std::make_shared<PLANE13Type>(PLANE13Type(properties));
        }
        else if (properties.element_type == "PLANE233")
        {
            elementTypes[pair.first] = std::make_shared<PLANE233Type>(PLANE233Type(properties));
        }
        // Add more conditions for other ElementType subclasses as needed
    }
}

VariableToElementIdsMapType &ElementTypesCollector::associateVariablesWithElementIds(std::string variables_set) const
{
    VariableToElementIdsMapType variableToElementIds;

    for (const auto &pair : elementTypes)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        for (const auto &variable : elementType.getActiveVariables(variables_set))
        {
            variableToElementIds[variable].insert(elementId);
        }
    }

    return variableToElementIds;
}