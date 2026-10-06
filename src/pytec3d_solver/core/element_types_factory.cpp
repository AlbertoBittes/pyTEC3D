/**
 * @file element_types_factory.cpp
 * @brief ElementTypesFactory: maps the Python element_types DataFrame onto element types.
 */
#include "element_types_factory.h"

#include <string>

#include "logging.h"

namespace py = pybind11;

namespace
{
    // ELEMENT_TYPE: 0 = 2D nodal, 1 = 3D nodal, 2 = 3D edge. Other formulations are other types of the JSON.
    std::string base_type_name(unsigned int element_number)
    {
        switch (element_number)
        {
        case 0:
            return "AV2D_nodal";
        case 1:
            return "AV3D_nodal";
        case 2:
            return "AV3D_edge";
        default:
            pytec3d::logging::error_and_abort("Unknown ELEMENT_TYPE " + std::to_string(element_number) + ".");
        }
    }

    unsigned int read_uint(const py::object &row, const char *column)
    {
        return static_cast<unsigned int>(py::cast<int>(row.attr("__getitem__")(column)));
    }
} // namespace

ElementTypesFactory::ElementTypesFactory(const py::object &element_types_attribute,
                                         std::shared_ptr<const pytec3d::ElementTypeParser> parser)
    : parser_(std::move(parser))
{
    if (!parser_)
        pytec3d::logging::error_and_abort("ElementTypesFactory: element type parser is not available.");
    if (element_types_attribute.is_none())
        return;

    const py::object iloc = element_types_attribute.attr("iloc");
    const std::size_t count = py::len(element_types_attribute);

    for (std::size_t i = 0; i < count; ++i)
    {
        const py::object row = iloc.attr("__getitem__")(py::int_(i));
        const unsigned int subdomain = read_uint(row, "TYPE_ID");
        const std::string name = base_type_name(read_uint(row, "ELEMENT_TYPE"));

        pytec3d::ElementTypeOptions options;
        options.order = read_uint(row, "OPTION2");
        options.axisymmetric = read_uint(row, "OPTION3") != 0;

        if (!parser_->has_type(name))
            pytec3d::logging::error_and_abort("TYPE_ID " + std::to_string(subdomain) + ": element type '" + name +
                                              "' is not defined in element_types.json.");

        const TypeKey key{name, options.order, options.axisymmetric};
        auto type = types_.find(key);
        if (type == types_.end())
            type = types_.emplace(key, parser_->build(name, options)).first;

        Entry entry;
        entry.type = type->second.get();
        entry.domain = read_uint(row, "OPTION1");
        if (entry.domain >= entry.type->domains().size())
            pytec3d::logging::error_and_abort("TYPE_ID " + std::to_string(subdomain) + ": OPTION1 " + std::to_string(entry.domain) +
                                              " is not a domain of domains.json.");
        entries_[static_cast<libMesh::subdomain_id_type>(subdomain)] = entry;
    }
}

const ElementTypesFactory::Entry &ElementTypesFactory::get(libMesh::subdomain_id_type subdomain_id) const
{
    const auto it = entries_.find(subdomain_id);
    if (it == entries_.end())
        pytec3d::logging::error_and_abort("Element type not found for subdomain ID: " + std::to_string(subdomain_id));
    return it->second;
}
