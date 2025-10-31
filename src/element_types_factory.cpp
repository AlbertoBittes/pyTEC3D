#include "element_types_factory.h"
#include "globals.h"

#include <cmath>
#include <regex>

namespace
{
unsigned int parse_element_type_number(const std::string &label)
{
    std::smatch match;
    if (std::regex_search(label, match, std::regex(R"(\d+)")))
        return static_cast<unsigned int>(std::stoul(match.str()));

    libmesh_error_msg("ElementTypesFactory: unable to parse element type from label '" << label << "'.");
}

unsigned int coerce_to_unsigned(const std::variant<int, Real, std::string> &value)
{
    if (std::holds_alternative<int>(value))
        return static_cast<unsigned int>(std::get<int>(value));
    if (std::holds_alternative<Real>(value))
        return static_cast<unsigned int>(std::llround(std::get<Real>(value)));
    if (std::holds_alternative<std::string>(value))
        return static_cast<unsigned int>(std::stoul(std::get<std::string>(value)));
    return 0u;
}
}

ElementTypesFactory::ElementTypesFactory()
    : FEDataObjectMap<DataframeMap>(std::make_shared<DataframeMap>()),
      element_types_map(std::make_shared<ElementTypesMap>())
{
}

ElementTypesFactory::ElementTypesFactory(std::shared_ptr<DataframeMap> table)
    : FEDataObjectMap<DataframeMap>(table ? table : std::make_shared<DataframeMap>()),
      element_types_map(std::make_shared<ElementTypesMap>())
{
    populate_element_types_map();
}

void ElementTypesFactory::populate_element_types_map()
{
    const DataframeMap *table = this->get_data_map();
    if (!table)
        libmesh_error_msg("ElementTypesFactory: data table not set.");

    for (const auto &row : *table)
    {
        if (row.size() < 2)
            libmesh_error_msg("ElementTypesFactory: expected at least two columns per row.");

        if (!std::holds_alternative<int>(row[0]))
            libmesh_error_msg("ElementTypesFactory: ansys_etype column must be integer.");
        if (!std::holds_alternative<std::string>(row[1]))
            libmesh_error_msg("ElementTypesFactory: element_type column must be string.");

        int ansys_type = std::get<int>(row[0]);
        unsigned int element_number = parse_element_type_number(std::get<std::string>(row[1]));

        std::vector<unsigned int> keyopts;
        keyopts.reserve(row.size());
        keyopts.push_back(element_number);

        for (std::size_t col = 2; col < row.size(); ++col)
            keyopts.push_back(coerce_to_unsigned(row[col]));

        (*element_types_map)[ansys_type] = create_element_type(
            static_cast<ElementTypeEnum>(element_number),
            keyopts);
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
// Variable → Element ID associations
// --------------------------------------------------------
std::unique_ptr<VariableToElementIdsMapType>
ElementTypesFactory::associate_variables_with_element_IDs(Globals::SystemType system_type) const
{
    auto variableToElementIds = std::make_unique<VariableToElementIdsMapType>();
    for (const auto &pair : *element_types_map)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        auto active_vars = elementType->get_active_variables(system_type);
        for (const auto &var : active_vars)
            (*variableToElementIds)[var].insert(elementId);
    }
    return variableToElementIds;
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
std::unique_ptr<ElementType> ElementTypesFactory::create_element_type(ElementTypeEnum type, const std::vector<unsigned int> &keyopts)
{
    switch (type)
    {
    case ElementTypeEnum::SOLID97:
        return std::make_unique<SOLID97Type>(keyopts);
    case ElementTypeEnum::SOLID236:
        return std::make_unique<SOLID236Type>(keyopts);
    case ElementTypeEnum::PLANE233:
        return std::make_unique<PLANE233Type>(keyopts);
    case ElementTypeEnum::INFIN110:
        return std::make_unique<INFIN110Type>(keyopts);
    case ElementTypeEnum::INFIN111:
        return std::make_unique<INFIN111Type>(keyopts);
    default:
        libmesh_error_msg("Unknown element type");
    }
}
