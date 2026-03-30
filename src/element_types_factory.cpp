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
    unsigned int parse_element_type_number(const std::string &label)
    {
        std::smatch match;
        if (std::regex_search(label, match, std::regex(R"(\d+)")))
            return static_cast<unsigned int>(std::stoul(match.str()));

        libmesh_error_msg("ElementTypesFactory: unable to parse element type from label '" << label << "'.");
    }

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
    const auto ansys_ids = py::cast<std::vector<int>>(
        dataframe.attr("__getitem__")("TYPE_ID").attr("tolist")());
    const auto element_labels = py::cast<std::vector<std::string>>(
        dataframe.attr("__getitem__")("element_type").attr("tolist")());

    std::vector<std::string> OPTION_columns;
    for (const py::handle &name_obj : dataframe.attr("columns"))
    {
        const std::string name = py::cast<std::string>(name_obj);
        if (name.rfind("OPTION", 0) == 0)
            OPTION_columns.push_back(name);
    }

    std::sort(OPTION_columns.begin(), OPTION_columns.end(), [](const std::string &a, const std::string &b)
              { return parse_option_column_index(a) < parse_option_column_index(b); });

    std::vector<std::vector<int>> option_values;
    option_values.reserve(OPTION_columns.size());
    for (const auto &column_name : OPTION_columns)
    {
        option_values.push_back(
            py::cast<std::vector<int>>(dataframe.attr("__getitem__")(column_name).attr("tolist")()));
    }

    const std::size_t count = ansys_ids.size();

    for (std::size_t row = 0; row < count; ++row)
    {
        const int ansys_type = ansys_ids[row];
        const unsigned int element_number = parse_element_type_number(element_labels[row]);

        std::vector<unsigned int> OPTIONs;
        OPTIONs.reserve(1 + OPTION_columns.size());
        OPTIONs.push_back(element_number);

        for (const auto &option_column_values : option_values)
            OPTIONs.push_back(static_cast<unsigned int>(option_column_values[row]));

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
    case ElementTypeEnum::SOLID97:
        return std::make_unique<SOLID97Type>(OPTIONs);
    case ElementTypeEnum::SOLID236:
        return std::make_unique<SOLID236Type>(OPTIONs);
    case ElementTypeEnum::PLANE233:
        return std::make_unique<PLANE233Type>(OPTIONs);
    default:
        libmesh_error_msg("Unknown element type");
    }
}
