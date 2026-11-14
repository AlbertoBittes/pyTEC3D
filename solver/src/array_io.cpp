/**
 * @file array_io.cpp
 * @brief The field arrays of the element types of a model.
 */
#include "array_io.h"

#include "element_type.h"
#include "logging.h"
#include "model_element_types.h"
#include "registries.h"


#include <algorithm>
#include <array>
#include <cctype>
#include <numeric>
#include <set>
#include <tuple>

namespace pytec3d
{
    namespace
    {
        std::size_t components_of(registry::FieldId field)
        {
            return registry::is_vector(field) ? 3 : 1;
        }
    } // namespace

    std::vector<FieldArrayRequirement> required_field_arrays(const ModelElementTypes &element_types)
    {
        std::vector<FieldArrayRequirement> out;
        std::set<std::pair<bool, std::string>> seen;
        const auto add = [&](registry::FieldId field, bool point, const std::string &prefix = "")
        {
            const std::string name = prefix + std::string(registry::name(field));
            if (seen.insert({point, name}).second)
                out.push_back({name, point, components_of(field)});
        };

        std::set<const ElementType *> types;
        for (const auto &[subdomain, entry] : element_types.entries())
            types.insert(entry.type);
        for (const ElementType *type : types)
        {
            for (std::size_t index = 0; index < type->spec().fields.size(); ++index)
            {
                const formulation::FieldSpec &field = type->field(index);
                switch (field.kind)
                {
                case FieldKind::Dof:
                    add(field.id, true);
                    add(field.id, true, "_");
                    if (field.nodal_load != registry::FieldId::COUNT)
                        add(field.nodal_load, true);
                    break;
                case FieldKind::Nodal:
                    add(field.id, true);
                    break;
                case FieldKind::Element:
                    add(field.id, false);
                    break;
                }
            }
            for (std::size_t p = 0; p < registry::kNumParameterId; ++p)
            {
                const auto id = static_cast<registry::ParameterId>(p);
                if (!registry::defined(type->id(), id))
                    continue;
                const registry::ParameterInfo &info = registry::info(id);
                if (info.source == registry::ParameterSource::Cell)
                    add(info.array, false);
                else if (info.source == registry::ParameterSource::Point)
                    add(info.array, true);
            }
        }
        return out;
    }

    void check_field_arrays(const ModelArrays &arrays, const std::vector<FieldArrayRequirement> &required,
                            std::size_t n_points, std::size_t n_cells, bool warn_unused)
    {
        std::string problems;
        std::set<std::pair<bool, std::string>> used;
        for (const FieldArrayRequirement &field : required)
        {
            used.insert({field.point, field.name});
            const char *where = field.point ? "point" : "cell";
            const RealArray *view = field.point ? arrays.point_field(field.name) : arrays.cell_field(field.name);
            const std::size_t n_tuples = field.point ? n_points : n_cells;
            if (!view)
                problems += "\n  - the " + std::string(where) + " field '" + field.name + "' is missing";
            else if (view->n_tuples() != n_tuples || view->n_components() != field.n_components)
                problems += "\n  - the " + std::string(where) + " field '" + field.name + "' has shape (" +
                            std::to_string(view->n_tuples()) + ", " + std::to_string(view->n_components()) +
                            "), expected (" + std::to_string(n_tuples) + ", " + std::to_string(field.n_components) + ")";
        }
        if (!problems.empty())
            logging::error_and_abort("Model input: the field arrays of the element types of the model:" + problems);

        if (!warn_unused)
            return;
        std::string unused;
        for (const auto &[name, view] : arrays.point_fields)
            if (!used.count({true, name}))
                unused += (unused.empty() ? "" : ", ") + name + " (point)";
        for (const auto &[name, view] : arrays.cell_fields)
            if (!used.count({false, name}))
                unused += (unused.empty() ? "" : ", ") + name + " (cell)";
        if (!unused.empty())
            logging::warn("Field arrays not used by the element types of the model, ignored: {}.", unused);
    }

} // namespace pytec3d
