/**
 * @file element_parameters.cpp
 * @brief ElementParameterValues: evaluation of the parameters and laws on one element.
 */
#include "element_parameters.h"

#include <algorithm>
#include <cmath>

#include <vtkCellData.h>
#include <vtkDataArray.h>

#include "bh_curve.h"
#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using libMesh::Real;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementParameterValues: " + message);
        }

        // Isotropic H(B) and its Jacobian dH/dB from a BH curve.
        void evaluate_bh(const BHCurve &curve, const Real B[3], unsigned int dim, Real H[3], Real J[3][3])
        {
            Real b_sq = 0.0;
            for (unsigned int k = 0; k < dim; ++k)
                b_sq += B[k] * B[k];
            const Real b_mag = std::sqrt(b_sq);
            const Real nu = curve.evaluate_reluctivity_from_b(b_mag);
            const Real dhdb = curve.slope_h_vs_b(b_mag);

            for (unsigned int k = 0; k < 3; ++k)
            {
                H[k] = 0.0;
                for (unsigned int l = 0; l < 3; ++l)
                    J[k][l] = 0.0;
            }
            constexpr Real kZeroB = 1e-14;
            for (unsigned int k = 0; k < dim; ++k)
            {
                H[k] = nu * B[k];
                J[k][k] += nu;
                if (b_mag > kZeroB)
                    for (unsigned int l = 0; l < dim; ++l)
                        J[k][l] += (dhdb - nu) * B[k] * B[l] / b_sq;
            }
        }
    } // namespace

    ElementParameterValues::ElementParameterValues(const ElementType &type, const Inputs &inputs)
        : type_(type), inputs_(inputs), values_(type.parameters().size()), computed_(type.parameters().size(), 0)
    {
    }

    ParameterValue ElementParameterValues::value(const ParameterRef &reference)
    {
        switch (reference.kind)
        {
        case ParameterKind::One:
            return scalar_value(reference.sign);
        case ParameterKind::Parameter:
            return scaled(symbol(reference.index), reference.sign);
        case ParameterKind::Law:
            return scaled(symbol(static_cast<unsigned int>(type_.law(reference.law()).linear)), reference.sign);
        }
        fail("unknown parameter kind.");
    }

    bool ElementParameterValues::nonlinear(const ParameterRef &reference) const
    {
        return reference.kind == ParameterKind::Law && inputs_.bh_curve != nullptr;
    }

    void ElementParameterValues::apply_law(const ParameterRef &reference, const Real B[3], Real H[3], Real J[3][3]) const
    {
        switch (registry::kind(reference.law()))
        {
        case registry::LawKind::BHCurve:
            if (!inputs_.bh_curve)
                fail("law '" + std::string(reference.name()) + "' needs a BH curve on the element material.");
            evaluate_bh(*inputs_.bh_curve, B, type_.dimension(), H, J);
            break;
        }
        for (unsigned int k = 0; k < 3; ++k)
        {
            H[k] *= reference.sign;
            for (unsigned int l = 0; l < 3; ++l)
                J[k][l] *= reference.sign;
        }
    }

    const ParameterValue &ElementParameterValues::symbol(unsigned int index)
    {
        if (!computed_.at(index))
        {
            const ParameterValue computed = compute(type_.parameters()[index]);
            values_[index] = computed;
            computed_[index] = 1;
        }
        return values_[index];
    }

    ParameterValue ElementParameterValues::compute(const ParameterDefinition &definition)
    {
        const registry::ParameterInfo &info = definition.info();
        switch (info.source)
        {
        case ParameterSource::Material:
        {
            const ParameterValue *value = inputs_.material ? inputs_.material->find(static_cast<registry::MaterialProperty>(info.property)) : nullptr;
            if (!value)
                undefined(info, "material", inputs_.material_id, inputs_.material != nullptr);
            return *value;
        }
        case ParameterSource::Constant:
        {
            const ParameterValue *value = inputs_.constants ? inputs_.constants->find(static_cast<registry::RealConstant>(info.property)) : nullptr;
            if (!value)
                undefined(info, "real constants", inputs_.constants_id, inputs_.constants != nullptr);
            return *value;
        }
        case ParameterSource::Cell:
            return read_cell(info);
        case ParameterSource::Point:
            fail("'" + std::string(info.key) + "' is point data: it can only load a DOF field.");
        case ParameterSource::Expression:
        {
            ParameterValue result = definition.expression->evaluate(*this);
            if (result.type != info.type)
                fail("'" + std::string(info.key) + "' is defined as a " + type_name(info.type) + " but its expression is a " +
                     type_name(result.type) + ".");
            return result;
        }
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
        vtkDataArray *array = inputs_.cell_data ? inputs_.cell_data->GetArray(array_name.c_str()) : nullptr;
        const unsigned int n = array ? static_cast<unsigned int>(array->GetNumberOfComponents()) : 0;
        const ParameterType type = info.type;
        if (array && n != (type == ParameterType::Scalar ? 1u : 3u))
            fail("cell array '" + array_name + "' has " + std::to_string(n) + " components; the " + type_name(type) + " needs " +
                 (type == ParameterType::Scalar ? "1." : "3."));

        // A missing array reads as zero.
        ParameterValue out;
        out.type = type;
        if (!array)
            return out;

        if (type == ParameterType::Scalar)
            out.scalar = array->GetComponent(inputs_.element_id, 0);
        else
            for (unsigned int c = 0; c < 3; ++c)
                out.vector(c) = array->GetComponent(inputs_.element_id, static_cast<int>(c));
        return out;
    }
} // namespace pytec3d
