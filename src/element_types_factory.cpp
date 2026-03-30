#include "element_types_factory.h"

#include <algorithm>
#include <cmath>
#include <regex>
#include <string>
#include <vector>

#include <pybind11/stl.h>

namespace py = pybind11;

namespace
{
    unsigned int parse_option_column_index(const std::string &name)
    {
        std::smatch match;
        if (std::regex_search(name, match, std::regex(R"(\d+)")))
            return static_cast<unsigned int>(std::stoul(match.str()));
        return 0;
    }
}

void ElementTypesFactory::load_from_python(const py::object &element_types_attribute)
{
    if (!element_types_map)
        element_types_map = std::make_shared<ElementTypesMap>();
    element_types_map->clear();

    if (element_types_attribute.is_none())
        return;

    py::object dataframe = element_types_attribute;

    std::vector<std::string> OPTION_columns;
    for (const py::handle &name_obj : dataframe.attr("columns"))
    {
        const std::string name = py::cast<std::string>(name_obj);
        if (name.rfind("OPTION", 0) == 0)
            OPTION_columns.push_back(name);
    }

    std::sort(OPTION_columns.begin(), OPTION_columns.end(), [](const std::string &a, const std::string &b)
              { return parse_option_column_index(a) < parse_option_column_index(b); });

    const py::object dataframe_iloc = dataframe.attr("iloc");
    const std::size_t count = py::len(dataframe);

    for (std::size_t row = 0; row < count; ++row)
    {
        const py::object row_data = dataframe_iloc.attr("__getitem__")(py::int_(row));
        const int ansys_type = py::cast<int>(row_data.attr("__getitem__")("TYPE_ID"));
        const unsigned int element_number =
            static_cast<unsigned int>(py::cast<int>(row_data.attr("__getitem__")("ELEMENT_TYPE")));

        std::vector<unsigned int> OPTIONs;
        OPTIONs.reserve(1 + OPTION_columns.size());
        OPTIONs.push_back(element_number);

        for (const auto &column_name : OPTION_columns)
        {
            OPTIONs.push_back(
                static_cast<unsigned int>(py::cast<int>(row_data.attr("__getitem__")(column_name))));
        }

        (*element_types_map)[ansys_type] = create_element_type(
            static_cast<ElementTypeEnum>(element_number),
            OPTIONs);
    }
}

// --------------------------------------------------------
// ElementType lookup
// --------------------------------------------------------
ElementType *ElementTypesFactory::get_element_type(subdomain_id_type element_id)
{
    auto it = element_types_map->find(element_id);
    if (it != element_types_map->end())
        return it->second.get();
    libmesh_error_msg("Element type not found for subdomain ID: " + std::to_string(element_id));
}

// --------------------------------------------------------
// Add element type
// --------------------------------------------------------
void ElementTypesFactory::add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type)
{
    (*element_types_map)[element_id] = std::move(element_type);
}

// --------------------------------------------------------
// Create element type instance
// --------------------------------------------------------
std::unique_ptr<ElementType> ElementTypesFactory::create_element_type(ElementTypeEnum type, const std::vector<unsigned int> &OPTIONs)
{
    switch (type)
    {
    case ElementTypeEnum::AV3DNODAL:
        return std::make_unique<AV3DNODALType>(OPTIONs);
    case ElementTypeEnum::AV3DEDGE:
        return std::make_unique<AV3DEDGEType>(OPTIONs);
    case ElementTypeEnum::AV2DNODAL:
        return std::make_unique<AV2DNODALType>(OPTIONs);
    default:
        libmesh_error_msg("Unknown element type");
    }
}
