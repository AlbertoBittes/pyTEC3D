#include "element_types_factory.h"
#include "globals.h"
#include <regex>
#include <iostream>

ElementTypesFactory::ElementTypesFactory()
    : element_types_map(std::make_shared<ElementTypesMap>())
{
}

// --------------------------------------------------------
// Stage 1: Read field_data arrays
// --------------------------------------------------------
void ElementTypesFactory::populate_intermediate_map_from_vtk(vtkFieldData *field_data)
{
    if (!field_data)
        libmesh_error_msg("No vtkFieldData provided.");

    std::vector<std::string> array_names = {
        "ansys_etype",
        "element_type",
        "KEYOPT1", "KEYOPT2", "KEYOPT3", "KEYOPT4", "KEYOPT5"};

    this->read_from_vtk_field_data(field_data, array_names);
}

// --------------------------------------------------------
// Stage 2: Build final element_types_map
// --------------------------------------------------------
void ElementTypesFactory::populate_final_map()
{
    auto *intermediate_map = this->get_data_map();
    if (!intermediate_map)
        libmesh_error_msg("Intermediate map is empty. Call populate_intermediate_map_from_vtk first.");

    element_types_map->clear();

    for (const auto &[row_id, fields] : *intermediate_map)
    {
        if (fields.size() < 2)
        {
            libmesh_warning("Row " << row_id << " has too few fields; skipping.");
            continue;
        }

        // --- Extract ansys_type ---
        int ansys_type = std::get<int>(fields[0]);

        // --- Extract element_type string ---
        std::string elem_type_str = std::get<std::string>(fields[1]);

        // --- Parse numeric element type ---
        std::smatch match;
        int elem_type_number = 0;
        if (std::regex_search(elem_type_str, match, std::regex(R"(\d+)")))
            elem_type_number = std::stoi(match.str());
        else
            libmesh_warning("Failed to extract element type number from: " << elem_type_str);

        // --- Build KEYOPT vector ---
        std::vector<unsigned int> keyopts;
        keyopts.push_back(static_cast<unsigned int>(elem_type_number));

        for (size_t i = 2; i < fields.size(); ++i)
        {
            if (std::holds_alternative<int>(fields[i]))
                keyopts.push_back(static_cast<unsigned int>(std::get<int>(fields[i])));
            else
                keyopts.push_back(0);
        }

        // --- Construct and store ElementType ---
        element_types_map->emplace(ansys_type,
                                   ElementTypesFactory::create_element_type(
                                       static_cast<ElementTypeEnum>(elem_type_number),
                                       keyopts));
    }
}

// --------------------------------------------------------
// Broadcast intermediate data map
// --------------------------------------------------------
// void ElementTypesFactory::broadcast(const libMesh::Parallel::Communicator &comm)
// {
//     comm.broadcast<std::unordered_map<unsigned int,
//                                       std::vector<std::variant<
//                                           int, Real, std::string,
//                                           std::vector<std::vector<int>>,
//                                           std::vector<std::vector<Real>>,
//                                           std::vector<std::vector<std::string>>>>>>(*get_data_map());
// }

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

template class vtkFieldDataFEDataObject<
    std::unordered_map<unsigned int,
                       std::vector<std::variant<
                           int,
                           Real,
                           std::string,
                           std::vector<std::vector<int>>,
                           std::vector<std::vector<Real>>,
                           std::vector<std::vector<std::string>>>>>>;