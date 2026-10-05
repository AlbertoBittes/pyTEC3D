/**
 * @file element_types_factory.h
 * @brief Per-subdomain element types and kernels built from the Python model description.
 *
 * Each row of the Python @c element_types DataFrame selects a base type of the JSON registry
 * (ELEMENT_TYPE picks the name) and its options (OPTION1 domain, OPTION2 order, OPTION3 axisymmetry). The factory instantiates one ElementTypeInstance and ElementKernels
 * per subdomain.
 */
#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include <memory>
#include <unordered_map>

#include <libmesh/id_types.h>
#include <pybind11/pybind11.h>

#include "element_kernels.h"
#include "element_type_instance.h"
#include "element_type_registry.h"

class ElementTypesFactory
{
public:
    struct Entry
    {
        std::unique_ptr<pytec3d::ElementTypeInstance> type;
        std::unique_ptr<pytec3d::ElementKernels> kernels;
    };
    using EntriesMap = std::unordered_map<libMesh::subdomain_id_type, Entry>;

    /** @brief Builds the entries from the Python @c element_types DataFrame; aborts on unknown types. */
    ElementTypesFactory(const pybind11::object &element_types_attribute,
                        std::shared_ptr<const pytec3d::ElementTypeRegistry> registry);

    const EntriesMap &entries() const { return entries_; }

    /** @brief Entry of subdomain @p subdomain_id; aborts if the subdomain has no element type. */
    const Entry &get(libMesh::subdomain_id_type subdomain_id) const;

private:
    std::shared_ptr<const pytec3d::ElementTypeRegistry> registry_; ///< Must outlive the entries.
    EntriesMap entries_;
};

#endif // ELEMENT_TYPES_FACTORY_H
