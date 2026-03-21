#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"

#include <memory>
#include <regex>
#include <vector>

#include <pybind11/pybind11.h>

// -----------------------------
// Factory to build ElementTypes from Python-side element_types attribute
// -----------------------------
class ElementTypesFactory
{
public:
    explicit ElementTypesFactory() = default;

    explicit ElementTypesFactory(const pybind11::object &element_types_attribute)
    {
        this->load_from_python(element_types_attribute);
    }
    void load_from_python(const pybind11::object &element_types_attribute);

    // Accessors
    inline ElementTypesMap *get_element_types_map() const { return element_types_map.get(); }
    ElementType *get_element_type(subdomain_id_type element_id);
    // Add element type manually
    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    // Factory method
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &OPTIONs);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;
};

#endif // ELEMENT_TYPES_FACTORY_H
