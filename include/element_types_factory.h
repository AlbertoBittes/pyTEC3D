/**
 * @file element_types_factory.h
 * @brief Per-subdomain element types built from the Python model description.
 *
 * Each row of the Python @c element_types DataFrame selects a type of element_types.json (ELEMENT_TYPE
 * picks the name), its physical region (OPTION1, an index of domains.json) and its options (OPTION2
 * order, OPTION3 axisymmetry). The factory has the parser build one ElementType per distinct (name,
 * options), shared by the subdomains that use it.
 */
#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>

#include <libmesh/id_types.h>
#include <pybind11/pybind11.h>

#include "element_type.h"
#include "element_type_parser.h"

class ElementTypesFactory
{
public:
    struct Entry
    {
        const pytec3d::ElementType *type = nullptr;
        unsigned int domain = 0; ///< Physical region: index into the domains of the element type.
    };
    using EntriesMap = std::unordered_map<libMesh::subdomain_id_type, Entry>;

    /** @brief Builds the entries from the Python @c element_types DataFrame; aborts on unknown types. */
    ElementTypesFactory(const pybind11::object &element_types_attribute,
                        std::shared_ptr<const pytec3d::ElementTypeParser> parser);

    const EntriesMap &entries() const { return entries_; }

    /** @brief Entry of subdomain @p subdomain_id; aborts if the subdomain has no element type. */
    const Entry &get(libMesh::subdomain_id_type subdomain_id) const;

private:
    using TypeKey = std::tuple<std::string, unsigned int, bool>;

    std::shared_ptr<const pytec3d::ElementTypeParser> parser_;
    std::map<TypeKey, std::unique_ptr<pytec3d::ElementType>> types_;
    EntriesMap entries_;
};

#endif // ELEMENT_TYPES_FACTORY_H
