#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"
#include "fe_data_object.h"
// ElementTypeFactory holds all ElementTypes
class ElementTypesFactory : public DataframeMapFEDataObject<int>
{
public:
    ElementTypesFactory()
        : DataframeMapFEDataObject<int>(),
          element_types_map(std::make_shared<ElementTypesMap>())
    {
    }

    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    std::unique_ptr<VariableToElementIdsMapType> associate_variables_with_element_IDs(Globals::SystemType system_type) const;

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

    void populate_data_map_from_vtk(vtkFieldData *field_data);

    void populate_element_types_map_from_dataframe();

    inline ElementTypesMap *get_element_types_map() const
    {
        return element_types_map.get();
    }

private:
    std::shared_ptr<ElementTypesMap> element_types_map;
};
#endif // ELEMENT_TYPES_FACTORY_H