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
            return scaled(symbol(type_.laws().at(reference.index).linear), reference.sign);
        }
        fail("unknown parameter kind.");
    }

    bool ElementParameterValues::nonlinear(const ParameterRef &reference) const
    {
        return reference.kind == ParameterKind::Law && inputs_.bh_curve != nullptr;
    }

    void ElementParameterValues::apply_law(const ParameterRef &reference, const Real B[3], Real H[3], Real J[3][3]) const
    {
        switch (type_.laws().at(reference.index).type)
        {
        case LawType::BHCurve:
            if (!inputs_.bh_curve)
                fail("law '" + reference.name + "' needs a BH curve on the element material.");
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
        switch (definition.source)
        {
        case ParameterSource::Material:
            return read_row(inputs_.material, definition);
        case ParameterSource::Constant:
            return read_row(inputs_.constants, definition);
        case ParameterSource::Cell:
            return read_cell(definition);
        case ParameterSource::Point:
            fail("'" + definition.name + "' is point data: it can only load a DOF field.");
        case ParameterSource::Expression:
        {
            ParameterValue result = definition.expression->evaluate(*this);
            if (result.type != definition.type)
                fail("'" + definition.name + "' is defined as a " + type_name(definition.type) + " but its expression is a " +
                     type_name(result.type) + ".");
            return result;
        }
        }
        fail("unknown parameter source.");
    }

    ParameterValue ElementParameterValues::read_row(const std::vector<Real> *row, const ParameterDefinition &definition) const
    {
        if (!row)
            fail("'" + definition.name + "': the row is not available.");
        for (unsigned int column : definition.columns)
            if (column >= row->size())
                fail("'" + definition.name + "': column " + std::to_string(column) + " is not in the row.");

        ParameterValue out;
        out.type = definition.type;
        if (definition.type == ParameterType::Scalar)
            out.scalar = (*row)[definition.columns.at(0)];
        else
            for (unsigned int c = 0; c < 3; ++c)
                out.vector(c) = (*row)[definition.columns.at(c)];
        return out;
    }

    ParameterValue ElementParameterValues::read_cell(const ParameterDefinition &definition) const
    {
        vtkDataArray *array = inputs_.cell_data ? inputs_.cell_data->GetArray(definition.array.c_str()) : nullptr;
        const unsigned int n = array ? static_cast<unsigned int>(array->GetNumberOfComponents()) : 0;

        ParameterType type = definition.type;
        if (type == ParameterType::Unknown)
        {
            if (!array)
                fail("cell array '" + definition.array + "' is not in the grid, so the type of '" + definition.name + "' is unknown.");
            if (n == 1)
                type = ParameterType::Scalar;
            else if (n == 3)
                type = ParameterType::Vector;
            else if (n == 9)
                type = ParameterType::Tensor;
            else
                fail("cell array '" + definition.array + "' has " + std::to_string(n) + " components; expected 1, 3 or 9.");
        }

        // A missing array reads as zero.
        ParameterValue out;
        out.type = type;
        if (!array)
            return out;

        const auto component = [&](unsigned int c)
        { return c < n ? array->GetComponent(inputs_.element_id, static_cast<int>(c)) : 0.0; };
        if (type == ParameterType::Scalar)
            out.scalar = component(0);
        else if (type == ParameterType::Vector)
            for (unsigned int c = 0; c < 3; ++c)
                out.vector(c) = component(c);
        else
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    out.tensor(i, j) = component(3 * i + j);
        return out;
    }
} // namespace pytec3d
