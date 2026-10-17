/**
 * @file model_element_types.cpp
 * @brief ModelElementTypes: maps the element_types of the model data onto element types.
 */
#include "model_element_types.h"

#include <algorithm>
#include <string>

#include "logging.h"

ModelElementTypes::ModelElementTypes(const pytec3d::ModelData &model_data)
{
    for (const auto &[subdomain, spec] : model_data.element_types())
    {
        pytec3d::ElementTypeOptions options;
        options.order = spec.order;
        options.axisymmetric = spec.axisymmetric;

        const TypeKey key{spec.type, options.order, options.axisymmetric};
        auto type = types_.find(key);
        if (type == types_.end())
            type = types_.emplace(key, std::make_unique<pytec3d::ElementType>(spec.type, options)).first;

        Entry entry;
        entry.type = type->second.get();
        entry.domain = spec.domain;
        entries_[static_cast<libMesh::subdomain_id_type>(subdomain)] = entry;
    }
}

unsigned int ModelElementTypes::time_order() const
{
    unsigned int order = 0;
    for (const auto &pair : entries_)
        order = std::max(order, pair.second.type->spec().time_order);
    return order;
}

const ModelElementTypes::Entry &ModelElementTypes::get(libMesh::subdomain_id_type subdomain_id) const
{
    const auto it = entries_.find(subdomain_id);
    if (it == entries_.end())
        pytec3d::logging::error_and_abort("Element type not found for subdomain ID (TYPE): " + std::to_string(subdomain_id) +
                                          ". Define it in model.data[\"element_types\"].");
    return it->second;
}

