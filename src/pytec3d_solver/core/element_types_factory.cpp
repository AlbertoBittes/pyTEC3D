/**
 * @file element_types_factory.cpp
 * @brief ElementTypesFactory: maps the element_types of the model data onto element types.
 */
#include "element_types_factory.h"

#include <string>

#include "logging.h"

ElementTypesFactory::ElementTypesFactory(const pytec3d::ModelData &model_data,
                                         std::shared_ptr<const pytec3d::ElementTypeParser> parser)
    : parser_(std::move(parser))
{
    if (!parser_)
        pytec3d::logging::error_and_abort("ElementTypesFactory: element type parser is not available.");

    for (const auto &[subdomain, spec] : model_data.element_types())
    {
        pytec3d::ElementTypeOptions options;
        options.order = spec.order;
        options.axisymmetric = spec.axisymmetric;

        const TypeKey key{spec.type, options.order, options.axisymmetric};
        auto type = types_.find(key);
        if (type == types_.end())
            type = types_.emplace(key, parser_->build(spec.type, options)).first;

        Entry entry;
        entry.type = type->second.get();
        entry.domain = spec.domain;
        entries_[static_cast<libMesh::subdomain_id_type>(subdomain)] = entry;
    }
}

const ElementTypesFactory::Entry &ElementTypesFactory::get(libMesh::subdomain_id_type subdomain_id) const
{
    const auto it = entries_.find(subdomain_id);
    if (it == entries_.end())
        pytec3d::logging::error_and_abort("Element type not found for subdomain ID (TYPE): " + std::to_string(subdomain_id) +
                                          ". Define it in model.data[\"element_types\"].");
    return it->second;
}
