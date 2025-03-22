#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <string>
#include <memory>
#include <vector>

struct ElementTypeProperties
{
    std::string ansys_element_type;
    std::vector<int> keyopts;
};

typedef std::unordered_map<int, ElementTypeProperties> ElementTypesMap;

std::shared_ptr<ElementTypesMap> readElementTypes(const std::string &filename);

#endif // ELEMENT_TYPES_H