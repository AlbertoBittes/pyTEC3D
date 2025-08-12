#include "element_types_factory.h"
#include "system_type.h"
#include "globals.h"
#include <regex>

void ElementTypesFactory::add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type)
{
    (*element_types_map)[element_id] = std::move(element_type);
}

std::unique_ptr<VariableToElementIdsMapType>
ElementTypesFactory::associate_variables_with_element_IDs(Globals::VariableCategory variable_category) const
{
    auto variableToElementIds = std::make_unique<VariableToElementIdsMapType>();
    for (const auto &pair : *element_types_map)
    {
        subdomain_id_type elementId = pair.first;
        const auto &elementType = pair.second;

        auto active_vars = elementType->get_active_variables(variable_category);
        for (const auto &var : active_vars)
        {
            (*variableToElementIds)[var].insert(elementId);
        }
    }
    return variableToElementIds;
}

std::unique_ptr<ElementType> ElementTypesFactory::create_element_type(ElementTypeEnum type, const std::vector<unsigned int> &keyopts)
{
    // Convert the string type to ElementTypeEnum
    switch (type)
    {
    case ElementTypeEnum::SOLID97:
        return std::make_unique<SOLID97Type>(keyopts);
    case ElementTypeEnum::SOLID236:
        return std::make_unique<SOLID236Type>(keyopts);
    case ElementTypeEnum::PLANE13:
        return std::make_unique<PLANE13Type>(keyopts);
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

void ElementTypesFactory::populate_data_map_from_vtk(vtkFieldData *field_data)
{
    if (!field_data)
        libmesh_error_msg("No vtkFieldData provided to populate_element_types_from_vtk.");

    // Retrieve the arrays from vtkFieldData
    vtkIntArray *ansys_element_type_num = vtkIntArray::SafeDownCast(field_data->GetArray("ansys_etype"));
    vtkStringArray *element_type = vtkStringArray::SafeDownCast(field_data->GetAbstractArray("element_type"));

    // Collect all KEYOPT arrays dynamically
    std::vector<vtkIntArray *> keyopt_arrays;
    std::vector<std::string> keyopt_names;
    int num_field_arrays = field_data->GetNumberOfArrays();
    for (int i = 0; i < num_field_arrays; ++i)
    {
        std::string arr_name = field_data->GetArrayName(i);
        std::smatch match;
        if (std::regex_match(arr_name, match, std::regex(R"(KEYOPT\d+)")))
        {
            vtkIntArray *keyopt_arr = vtkIntArray::SafeDownCast(field_data->GetArray(arr_name.c_str()));
            if (keyopt_arr)
            {
                keyopt_arrays.push_back(keyopt_arr);
                keyopt_names.push_back(arr_name);
            }
        }
    }

    // Sort keyopt_names and keyopt_arrays by KEYOPT number for consistency
    std::vector<std::pair<std::string, vtkIntArray *>> keyopt_pairs;
    for (size_t i = 0; i < keyopt_names.size(); ++i)
        keyopt_pairs.emplace_back(keyopt_names[i], keyopt_arrays[i]);
    std::sort(keyopt_pairs.begin(), keyopt_pairs.end(),
              [](const auto &a, const auto &b)
              {
                  int anum = std::stoi(a.first.substr(6));
                  int bnum = std::stoi(b.first.substr(6));
                  return anum < bnum;
              });
    keyopt_arrays.clear();
    for (const auto &p : keyopt_pairs)
        keyopt_arrays.push_back(p.second);

    // Validate that required arrays are present
    if (!ansys_element_type_num || !element_type || keyopt_arrays.empty())
    {
        libmesh_error_msg("Missing required vtkFieldData arrays: ansys_etype, element_type, or KEYOPT arrays.");
        return;
    }

    // Ensure all arrays have the same number of tuples
    vtkIdType num_tuples = ansys_element_type_num->GetNumberOfTuples();
    if (element_type->GetNumberOfTuples() != num_tuples)
    {
        libmesh_error_msg("Mismatch in the number of tuples between ansys_etype and element_type arrays.");
        return;
    }
    for (auto arr : keyopt_arrays)
    {
        if (arr->GetNumberOfTuples() != num_tuples)
        {
            libmesh_error_msg("Mismatch in the number of tuples between KEYOPT arrays.");
            return;
        }
    }

    // --- Populate data_map ---
    this->data_map->clear();
    for (vtkIdType i = 0; i < num_tuples; ++i)
    {
        std::vector<std::variant<int, Real, std::string>> row;
        // ansys_etype (int)
        row.push_back(static_cast<int>(ansys_element_type_num->GetValue(i)));
        // element_type (string)
        row.push_back(element_type->GetValue(i));
        // KEYOPT values (int)
        for (auto arr : keyopt_arrays)
            row.push_back(static_cast<int>(arr->GetValue(i)));
        this->data_map->push_back(std::move(row));
    }

    // --- Populate dataframe_map from data_map ---
    this->populate_dataframe_map();
}

void ElementTypesFactory::populate_element_types_map_from_dataframe()
{
    element_types_map->clear();

    auto *df_map = this->get_dataframe_map();
    if (!df_map)
        libmesh_error_msg("No dataframe_map available.");

    for (const auto &pair : *df_map)
    {
        int ansys_type = pair.first;
        const std::vector<int> &keyopts_int = pair.second;

        if (keyopts_int.empty())
            libmesh_error_msg("Empty keyopts vector in dataframe_map.");

        // The first entry is the element type number
        unsigned int elem_type_number = static_cast<unsigned int>(keyopts_int[0]);

        // The rest (if any) are KEYOPT values
        std::vector<unsigned int> keyopts;
        for (size_t j = 1; j < keyopts_int.size(); ++j)
            keyopts.push_back(static_cast<unsigned int>(keyopts_int[j]));

        add_element_type(ansys_type, create_element_type(static_cast<ElementTypeEnum>(elem_type_number), keyopts));
    }
}