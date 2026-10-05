/**
 * @file element_types_factory.cpp
 * @brief ElementTypesFactory: maps the Python element_types DataFrame onto registry element types.
 */
#include "element_types_factory.h"

#include <string>

#include "logging.h"

namespace py = pybind11;

namespace
{
    // ELEMENT_TYPE: 0 = 2D nodal, 1 = 3D nodal, 2 = 3D edge. Other formulations are other registry types.
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
                                         std::shared_ptr<const pytec3d::ElementTypeRegistry> registry)
    : registry_(std::move(registry))
{
    if (!registry_)
        pytec3d::logging::error_and_abort("ElementTypesFactory: element type registry is not available.");
    if (element_types_attribute.is_none())
        return;

    const py::object iloc = element_types_attribute.attr("iloc");
    const std::size_t count = py::len(element_types_attribute);

    for (std::size_t i = 0; i < count; ++i)
    {
        const py::object row = iloc.attr("__getitem__")(py::int_(i));
        const unsigned int subdomain = read_uint(row, "TYPE_ID");

        pytec3d::ElementOptions options;
        options.type = base_type_name(read_uint(row, "ELEMENT_TYPE"));
        options.domain = read_uint(row, "OPTION1");
        options.order = read_uint(row, "OPTION2");
        options.axisymmetric = read_uint(row, "OPTION3") != 0;

        if (!registry_->has_type(options.type))
            pytec3d::logging::error_and_abort("TYPE_ID " + std::to_string(subdomain) + ": element type '" + options.type +
                                              "' is not defined in element_types.json.");

        Entry entry;
        entry.type = std::make_unique<pytec3d::ElementTypeInstance>(*registry_, options);
        entry.kernels = std::make_unique<pytec3d::ElementKernels>(*entry.type);
        entries_[static_cast<libMesh::subdomain_id_type>(subdomain)] = std::move(entry);
    }
}

const ElementTypesFactory::Entry &ElementTypesFactory::get(libMesh::subdomain_id_type subdomain_id) const
{
    const auto it = entries_.find(subdomain_id);
    if (it == entries_.end())
        pytec3d::logging::error_and_abort("Element type not found for subdomain ID: " + std::to_string(subdomain_id));
    return it->second;
}
