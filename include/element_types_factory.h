#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "apec_m_system.h"
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
    // Add element type manually
    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    void print_data_info() const override;

    // Factory method
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &OPTIONs);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;

    // Populates element_types_map from the copied data_map
    void populate_element_types_map();
};

#endif // ELEMENT_TYPES_FACTORY_H
