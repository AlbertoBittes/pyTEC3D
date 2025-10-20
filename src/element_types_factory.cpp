#include "element_types_factory.h"
#include "globals.h"
#include <regex>
#include <iostream>

ElementTypesFactory::ElementTypesFactory(std::shared_ptr<vtkFieldData> field_data)
    : FEDataObjectMap<vtkFieldData>(field_data),
      element_types_map(std::make_shared<ElementTypesMap>())
{
    populate_element_types_map();
}
void ElementTypesFactory::populate_element_types_map()
{
    vtkFieldData *field_data = this->get_data_map();
    if (!field_data)
        libmesh_error_msg("No vtkFieldData available.");

    // Get ansys_etype and element_type arrays
    vtkIntArray *ansys_etype_arr = vtkIntArray::SafeDownCast(field_data->GetArray("ansys_etype"));
    vtkStringArray *element_type_arr = vtkStringArray::SafeDownCast(field_data->GetAbstractArray("element_type"));

    // Collect all KEYOPT arrays and sort them numerically
    std::vector<std::pair<int, vtkIntArray *>> keyopt_arrays;
    int num_field_arrays = field_data->GetNumberOfArrays();
    for (int i = 0; i < num_field_arrays; ++i)
    {
        std::string arr_name = field_data->GetArrayName(i);
        std::smatch match;
        if (std::regex_match(arr_name, match, std::regex(R"(KEYOPT(\d+))")))
        {
            int keyopt_num = std::stoi(match.str(1));
            vtkIntArray *keyopt_arr = vtkIntArray::SafeDownCast(field_data->GetArray(arr_name.c_str()));
            if (keyopt_arr)
                keyopt_arrays.emplace_back(keyopt_num, keyopt_arr);
        }
    }
    std::sort(keyopt_arrays.begin(), keyopt_arrays.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    if (!ansys_etype_arr || !element_type_arr)
        libmesh_error_msg("Required arrays not found in vtkFieldData.");

    vtkIdType num_rows = ansys_etype_arr->GetNumberOfTuples();
    if (element_type_arr->GetNumberOfTuples() != num_rows)
        libmesh_error_msg("Mismatch in number of tuples between ansys_etype and element_type arrays.");

    for (const auto &pair : keyopt_arrays)
    {
        if (pair.second->GetNumberOfTuples() != num_rows)
            libmesh_error_msg("Mismatch in number of tuples in KEYOPT arrays.");
    }

    element_types_map->clear();

    for (vtkIdType i = 0; i < num_rows; ++i)
    {
        int ansys_type = ansys_etype_arr->GetValue(i);
        std::string elem_type_str = element_type_arr->GetValue(i);

        std::smatch match;
        int elem_type_number = 0;
        if (std::regex_search(elem_type_str, match, std::regex(R"(\d+)")))
            elem_type_number = std::stoi(match.str());
        else
            libmesh_warning("Failed to extract element type number from: " << elem_type_str);

        std::vector<unsigned int> keyopts;
        keyopts.push_back(static_cast<unsigned int>(elem_type_number));
        for (const auto &pair : keyopt_arrays)
            keyopts.push_back(static_cast<unsigned int>(pair.second->GetValue(i)));

        element_types_map->emplace(ansys_type,
                                   ElementTypesFactory::create_element_type(
                                       static_cast<ElementTypeEnum>(elem_type_number),
                                       keyopts));
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