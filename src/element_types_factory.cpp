#include "element_types_factory.h"

#include "variables_map.h"

#include <algorithm>
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

}

void ElementTypesFactory::read_data(std::ifstream &file)
{
    BinaryDataFEDataObject::read_data(file);
    populate_element_types_map();
}

void ElementTypesFactory::populate_element_types_map()
{
    if (element_types_map == nullptr)
        element_types_map = std::make_shared<ElementTypesMap>();
    element_types_map->clear();

    if (empty())
        return;

    if (!has_column("ansys_etype"))
        libmesh_error_msg("ElementTypesFactory: missing 'ansys_etype' column.");
    if (!has_column("element_type"))
        libmesh_error_msg("ElementTypesFactory: missing 'element_type' column.");

    const auto &ansys_ids = get_int_column("ansys_etype");
    const auto &element_labels = get_string_column("element_type");
    const std::size_t count = entry_count();

    if (ansys_ids.size() < count || element_labels.size() < count)
        libmesh_error_msg("ElementTypesFactory: inconsistent element type data size.");

    std::vector<std::string> keyopt_columns;
    keyopt_columns.reserve(this->get_data_map()->columns.size());
    for (const auto &pair : this->get_data_map()->columns)
    {
        if (pair.first.rfind("KEYOPT", 0) == 0)
            keyopt_columns.push_back(pair.first);
    }

    std::sort(keyopt_columns.begin(), keyopt_columns.end(), [](const std::string &a, const std::string &b)
              {
        auto extract_index = [](const std::string &name) {
            auto pos = name.find_first_of("0123456789");
            return pos == std::string::npos ? 0 : std::stoul(name.substr(pos));
        };
        return extract_index(a) < extract_index(b); });

    for (std::size_t row = 0; row < count; ++row)
    {
        int ansys_type = ansys_ids[row];
        unsigned int element_number = parse_element_type_number(element_labels[row]);

        std::vector<unsigned int> keyopts;
        keyopts.reserve(1 + keyopt_columns.size());
        keyopts.push_back(element_number);

        for (const auto &column_name : keyopt_columns)
            keyopts.push_back(static_cast<unsigned int>(get_int(column_name, row)));

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
ElementTypesFactory::associate_variables_with_element_IDs(ApecMSystemBase::SystemType system_type) const
{
    auto variableToElementIds = std::make_unique<VariableToElementIdsMapType>();
    for (const auto &pair : *element_types_map)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;
        const auto vars = VariablesMap::get_variables_for_element(*elementType);
        for (const auto &var : vars)
        {
            if (VariablesMap::get_info_by_key(var).system_type == system_type)
                (*variableToElementIds)[var].insert(elementId);
        }
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

void ElementTypesFactory::print_data_info() const
{
    if (!element_types_map || element_types_map->empty())
    {
        libMesh::out << "  Element types: none loaded." << std::endl;
        return;
    }

    libMesh::out << "  Element type definitions: " << element_types_map->size() << std::endl;
    for (const auto &entry : *element_types_map)
    {
        libMesh::out << "    ANSYS etype " << entry.first;
        const auto &element_type = entry.second;
        if (element_type)
        {
            const auto &keyopts = element_type->get_keyopts();
            if (!keyopts.empty())
                libMesh::out << " (Element " << keyopts.front() << ")";
            if (keyopts.size() > 1)
            {
                libMesh::out << " KEYOPTs=[";
                for (std::size_t idx = 1; idx < keyopts.size(); ++idx)
                {
                    if (idx > 1)
                        libMesh::out << ", ";
                    libMesh::out << keyopts[idx];
                }
                libMesh::out << "]";
            }
        }
        libMesh::out << std::endl;
    }
}
