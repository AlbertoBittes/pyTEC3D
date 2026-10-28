/**
 * @file element_parameters.cpp
 * @brief ElementParameterValues: evaluation of the parameters on one element.
 */
#include "element_parameters.h"

#include <algorithm>
#include <cmath>

#include "logging.h"

namespace pytec3d
{
    using registry::ParameterSource;

    namespace
    {
        using libMesh::Real;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementParameterValues: " + message);
        }

        // A vector or tensor given in the axes of the element, in the global axes.
        ParameterValue to_global(const ParameterValue &value, const libMesh::TensorValue<Real> &rotation)
        {
            switch (value.type())
            {
            case ParameterType::Vector:
                return ParameterValue(libMesh::RealVectorValue(rotation * value.vector()));
            case ParameterType::Tensor:
                return ParameterValue(libMesh::RealTensorValue(rotation * value.tensor() * rotation.transpose()));
            case ParameterType::Scalar:
                break;
            }
            return value;
        }
    } // namespace

    ElementParameterValues::ElementParameterValues(const ElementType &type, const Inputs &inputs)
        : type_(type), inputs_(inputs), values_(registry::kNumParameterId), computed_(registry::kNumParameterId, 0)
    {
    }

    const ParameterValue &ElementParameterValues::symbol(registry::ParameterId id)
    {
        const std::size_t index = static_cast<std::size_t>(id);
        if (!computed_.at(index))
        {
            if (!registry::defined(type_.id(), id, inputs_.domain))
                fail("'" + std::string(registry::name(id)) + "' is not defined for element type " + std::string(type_.name()) +
                     " in domain " + std::string(inputs_.domain == registry::Domain::COUNT ? "?" : registry::name(inputs_.domain)) + ".");
            const ParameterValue computed = compute(id);
            values_[index] = computed;
            computed_[index] = 1;
        }
        return values_[index];
    }

    ParameterValue ElementParameterValues::compute(registry::ParameterId id)
    {
        const registry::ParameterInfo &info = registry::info(id);
        switch (info.source)
        {
        case ParameterSource::Material:
        {
            const ParameterValue *value = inputs_.material ? inputs_.material->find(static_cast<registry::MaterialProperty>(info.property)) : nullptr;
            if (value)
                return to_global(*value, inputs_.rotation);
            if (const auto fallback = registry::default_value(type_.id(), id))
                return to_global(*fallback, inputs_.rotation);
            undefined(info, "material", inputs_.material_id, inputs_.material != nullptr);
        }
        case ParameterSource::Constant:
        {
            const ParameterValue *value = inputs_.constants ? inputs_.constants->find(static_cast<registry::RealConstant>(info.property)) : nullptr;
            if (value)
                return to_global(*value, inputs_.rotation);
            if (const auto fallback = registry::default_value(type_.id(), id))
                return to_global(*fallback, inputs_.rotation);
            undefined(info, "real constants", inputs_.constants_id, inputs_.constants != nullptr);
        }
        case ParameterSource::Cell:
            return to_global(read_cell(info), inputs_.rotation);
        case ParameterSource::Point:
            fail("'" + std::string(registry::name(id)) + "' is point data: it can only load a DOF field.");
        case ParameterSource::Expression:
            return registry::evaluate_expression(type_.id(), id, *this);
        }
        fail("unknown parameter source.");
    }

    void ElementParameterValues::undefined(const registry::ParameterInfo &info, const char *what, unsigned int id, bool row_exists) const
    {
        const bool is_material = info.source == ParameterSource::Material;
        const std::string_view property = is_material ? registry::name(static_cast<registry::MaterialProperty>(info.property))
                                                      : registry::name(static_cast<registry::RealConstant>(info.property));
        fail("'" + std::string(property) + "' is needed by element type " + std::string(type_.name()) + " in domain " +
             std::string(registry::name(inputs_.domain)) + " but " + what + " " + std::to_string(id) +
             (row_exists ? " does not define it." : " is not in the model data."));
    }

    ParameterValue ElementParameterValues::read_cell(const registry::ParameterInfo &info) const
    {
        const std::string array_name(registry::name(info.array));
        const RealArray *array = inputs_.arrays ? inputs_.arrays->cell_field(array_name) : nullptr;
        const unsigned int n = array ? static_cast<unsigned int>(array->n_components()) : 0;
        const ParameterType type = info.type;
        if (array && n != (type == ParameterType::Scalar ? 1u : 3u))
            fail("cell array '" + array_name + "' has " + std::to_string(n) + " components; the " + type_name(type) + " needs " +
                 (type == ParameterType::Scalar ? "1." : "3."));

        // A missing array reads as zero.
        ParameterValue out = ParameterValue::zero(type);
        if (!array)
            return out;

        if (type == ParameterType::Scalar)
            out.scalar() = (*array)(inputs_.element_id);
        else
            for (unsigned int c = 0; c < 3; ++c)
                out.vector()(c) = (*array)(inputs_.element_id, c);
        return out;
    }
} // namespace pytec3d
