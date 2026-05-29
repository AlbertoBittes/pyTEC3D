/**
 * @file element_types_factory.h
 * @brief Factory that constructs ElementType instances from the Python model description.
 *
 * Parses the @c element_types pandas DataFrame from the APEC Python object,
 * reads the OPTION columns, and instantiates the appropriate ElementType subclass
 * for each subdomain ID via a static factory method.
 */
#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"

#include <memory>
#include <regex>
#include <vector>

#include <pybind11/pybind11.h>

/**
 * @brief Builds and owns the per-subdomain ElementType map for the current model.
 *
 * Populated from the Python @c element_types attribute by load_from_python().
 * ElementBuilder and ApecMSystem query this factory to obtain element-type
 * objects by subdomain ID.
 */
class ElementTypesFactory
{
public:
    explicit ElementTypesFactory() = default;

    explicit ElementTypesFactory(const pybind11::object &element_types_attribute)
    {
        this->load_from_python(element_types_attribute);
    }
    /** @brief Populates the element-type map from a Python @c element_types DataFrame. */
    void load_from_python(const pybind11::object &element_types_attribute);

    /** @brief Returns the raw element-type map pointer (may be null before loading). */
    inline ElementTypesMap *get_element_types_map() const { return element_types_map.get(); }

    /** @brief Looks up the ElementType for subdomain @p element_id; aborts if not found. */
    ElementType *get_element_type(libMesh::subdomain_id_type element_id);

    /** @brief Manually inserts an element type, bypassing Python loading. */
    void add_element_type(libMesh::subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    /**
     * @brief Instantiates the correct ElementType subclass for the given @p type and @p OPTIONs.
     * @param type    Enum identifying the formulation variant.
     * @param OPTIONs Option array (OPTION1..N) from the model description.
     */
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &OPTIONs);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;
};

#endif // ELEMENT_TYPES_FACTORY_H
