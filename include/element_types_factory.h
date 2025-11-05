#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"
#include "fe_data_object.h"
#include <memory>
#include <regex>
#include <vector>

// -----------------------------
// Factory to build ElementTypes from serialized columnar data
// -----------------------------
class ElementTypesFactory : public BinaryDataFEDataObject
{
public:
    explicit ElementTypesFactory() = default;

    explicit ElementTypesFactory(std::ifstream &file) { this->read_data(file); }

    void read_data(std::ifstream &file) override;

    // Accessors
    inline ElementTypesMap *get_element_types_map() const { return element_types_map.get(); }
    ElementType *get_element_type(subdomain_id_type element_id);
    std::unique_ptr<VariableToElementIdsMapType> associate_variables_with_element_IDs(Globals::SystemType system_type) const;

    // Add element type manually
    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    // Factory method
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;

    // Populates element_types_map from the copied data_map
    void populate_element_types_map();
};

#endif // ELEMENT_TYPES_FACTORY_H
