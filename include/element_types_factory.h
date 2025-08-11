#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"
#include "fe_data_object.h"
// ElementTypeFactory holds all ElementTypes
class ElementTypesFactory : public FEDataObjectMap<ElementTypesMap>
{
public:
    ElementTypesFactory() = default;

    // Factory method to create an ElementType object
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    std::unique_ptr<VariableToElementIdsMapType> associate_variables_with_element_IDs(Globals::VariableCategory variable_category) const;

    inline ElementType *get_element_type(subdomain_id_type element_id)
    {
        auto it = data_map->find(element_id);
        if (it != data_map->end())
        {
            return it->second.get();
        }
        else
        {
            libmesh_error_msg("Element type not found for subdomain ID: " + std::to_string(element_id));
        }
    }

    void populate_element_types_from_vtk(vtkFieldData *field_data);
};
#endif // ELEMENT_TYPES_FACTORY_H