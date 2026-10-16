/**
 * @file element_types_factory.h
 * @brief Per-subdomain element types built from the "element_types" section of the model data.
 *
 * Each TYPE id of the mesh selects an element type of element_types.json (element_type), its physical
 * region (option1, a domain of domains.json) and its options (option2 order, option3 axisymmetry). The factory
 * builds one ElementType per distinct (type, options), shared by the subdomains that use it.
 */
#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>

#include <libmesh/id_types.h>

#include "element_type.h"
#include "model_data.h"

class ElementTypesFactory
{
public:
    struct Entry
    {
        const pytec3d::ElementType *type = nullptr;
        pytec3d::registry::Domain domain = pytec3d::registry::Domain::COUNT; ///< Physical region of the subdomain.
    };
    using EntriesMap = std::unordered_map<libMesh::subdomain_id_type, Entry>;

    explicit ElementTypesFactory(const pytec3d::ModelData &model_data);

    const EntriesMap &entries() const { return entries_; }

    /** @brief Entry of subdomain @p subdomain_id; aborts if the subdomain has no element type. */
    const Entry &get(libMesh::subdomain_id_type subdomain_id) const;

    /** @brief Highest time order of the element types: the order of the time scheme. */
    unsigned int time_order() const;

private:
    using TypeKey = std::tuple<pytec3d::registry::ElementTypeId, unsigned int, bool>;

    std::map<TypeKey, std::unique_ptr<pytec3d::ElementType>> types_;
    EntriesMap entries_;
};

#endif // ELEMENT_TYPES_FACTORY_H
