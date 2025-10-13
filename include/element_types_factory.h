#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"
#include "fe_data_object.h"
// ElementTypeFactory holds all ElementTypes

class ElementTypesFactory
    : public vtkFieldDataFEDataObject<std::unordered_map<unsigned int, std::tuple<int, std::string, int, int, int, int, int>>>
{
public:
    ElementTypesFactory()
        : element_types_map(std::make_shared<ElementTypesMap>())
    {
    }

    // -------------------------------
    // Stage 1: Populate intermediate data_map from vtkFieldData
    // -------------------------------
    void populate_intermediate_map_from_vtk(vtkFieldData *field_data) override
    {
        if (!field_data)
            libmesh_error_msg("No vtkFieldData provided.");

        // Specify arrays to read
        std::vector<std::string> array_names = {
            "ansys_etype",
            "element_type",
            "KEYOPT1", "KEYOPT2", "KEYOPT3", "KEYOPT4", "KEYOPT5"};

        // Base class method handles the reading and filling data_map
        this->read_from_vtk_field_data(field_data, array_names);
    }

    // -------------------------------
    // Stage 2: Populate element_types_map from intermediate map
    // -------------------------------
    void populate_final_map() override
    {
        auto *intermediate_map = this->get_data_map();
        if (!intermediate_map)
            libmesh_error_msg("Intermediate map is empty. Call populate_intermediate_map_from_vtk first.");

        element_types_map->clear();

        for (const auto &[row_id, row_tuple] : *intermediate_map)
        {
            int ansys_type = std::get<0>(row_tuple);
            std::string elem_type_str = std::get<1>(row_tuple);

            // Extract element type number from string
            std::smatch match;
            int elem_type_number = 0;
            if (std::regex_search(elem_type_str, match, std::regex(R"(\d+)")))
                elem_type_number = std::stoi(match.str());
            else
                libmesh_warning("Failed to extract element type number from: " << elem_type_str);

            // Build KEYOPT vector
            std::vector<unsigned int> keyopts;
            keyopts.push_back(static_cast<unsigned int>(elem_type_number));

            // Extract KEYOPTs from tuple
            constexpr size_t tuple_size = std::tuple_size<decltype(row_tuple)>::value;
            for (size_t i = 2; i < tuple_size; ++i)
            {
                int kopt = 0;
                switch (i)
                {
                case 2:
                    kopt = std::get<2>(row_tuple);
                    break;
                case 3:
                    kopt = std::get<3>(row_tuple);
                    break;
                case 4:
                    kopt = std::get<4>(row_tuple);
                    break;
                case 5:
                    kopt = std::get<5>(row_tuple);
                    break;
                case 6:
                    kopt = std::get<6>(row_tuple);
                    break;
                default:
                    break;
                }
                keyopts.push_back(static_cast<unsigned int>(kopt));
            }

            // Construct and store ElementType
            element_types_map->emplace(ansys_type,
                                       ElementTypesFactory::create_element_type(static_cast<ElementTypeEnum>(elem_type_number), keyopts));
        }
    }

    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

    inline ElementTypesMap *get_element_types_map() const
    {
        return element_types_map.get();
    }
    inline ElementType *get_element_type(subdomain_id_type element_id)
    {
        auto it = element_types_map->find(element_id);
        if (it != element_types_map->end())
        {
            return it->second.get();
        }
        else
        {
            libmesh_error_msg("Element type not found for subdomain ID: " + std::to_string(element_id));
        }
    }
    std::unique_ptr<VariableToElementIdsMapType> associate_variables_with_element_IDs(Globals::SystemType system_type) const;
    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;
};
#endif // ELEMENT_TYPES_FACTORY_H