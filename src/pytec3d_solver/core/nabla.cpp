/**
 * @file nabla.cpp
 * @brief Nabla: operator matrices at quadrature points and weak-form accumulation.
 */
#include "nabla.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <utility>

#include <libmesh/elem.h>
#include <libmesh/fe.h>

#include "element_type_instance.h"
#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using libMesh::DenseMatrix;
        using libMesh::Real;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Nabla: " + message);
        }

        std::string strip_spaces(const std::string &text)
        {
            std::string out;
            for (char c : text)
                if (!std::isspace(static_cast<unsigned char>(c)))
                    out.push_back(c);
            return out;
        }

        bool is_identifier(const std::string &text)
        {
            if (text.empty() || !(std::isalpha(static_cast<unsigned char>(text[0])) || text[0] == '_'))
                return false;
            return std::all_of(text.begin(), text.end(), [](char c)
                               { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
        }

        const VariableSpec &find_spec(const ElementTypeInstance &element_type, const std::string &variable)
        {
            const VariableSpec *spec = element_type.find_variable(variable);
            if (!spec)
                fail("variable '" + variable + "' is not active in element type '" + element_type.options().type + "'.");
            return *spec;
        }

        // An out-of-plane scalar is the single z component of a vector field on a 2D element.
        bool is_out_of_plane_scalar(const VariableSpec &spec, unsigned int dimension)
        {
            return dimension == 2 && spec.physical_vector && spec.active_components.size() == 1 &&
                   spec.active_components[0] == 2;
        }

        // A lumped potential is a scalar that is not a vector component, on a 2D element.
        bool is_lumped_potential(const VariableSpec &spec, unsigned int dimension)
        {
            return dimension == 2 && !spec.physical_vector;
        }

        Real out_of_plane_length(const ElementTypeInstance &element_type, Real thickness, Real x)
        {
            return element_type.options().axisymmetric ? thickness * 2.0 * libMesh::pi * std::abs(x) : thickness;
        }

        Nabla::QuadratureMatrices scalar_operator(const ElementTypeInstance &element_type,
                                                  const VariableSpec &spec,
                                                  NablaOperator op,
                                                  const libMesh::Elem &elem,
                                                  Real thickness)
        {
            const unsigned int dim = element_type.dimension();
            const bool axisymmetric = element_type.options().axisymmetric;
            const bool lumped = is_lumped_potential(spec, dim);

            if (op == NablaOperator::Curl && !is_out_of_plane_scalar(spec, dim))
                fail("curl of scalar variable '" + spec.name + "' needs an out-of-plane (z) component on a 2D element.");
            if (op == NablaOperator::Div && !spec.physical_vector)
                fail("div needs a vector field, but '" + spec.name + "' is scalar.");

            libMesh::FEBase &fe = element_type.scalar_fe(spec.name);
            // libMesh only computes what was requested before reinit().
            const auto &phi = fe.get_phi();
            const auto *dphi = (op == NablaOperator::Curl || (op == NablaOperator::Grad && !lumped)) ? &fe.get_dphi() : nullptr;
            const auto &xyz = fe.get_xyz();
            fe.reinit(&elem);

            const unsigned int nq = fe.n_quadrature_points();
            const unsigned int n = static_cast<unsigned int>(phi.size());
            Nabla::QuadratureMatrices out(nq);

            for (unsigned int qp = 0; qp < nq; ++qp)
            {
                DenseMatrix<Real> &m = out[qp];
                switch (op)
                {
                case NablaOperator::Identity:
                    m.resize(1, n);
                    for (unsigned int i = 0; i < n; ++i)
                        m(0, i) = phi[i][qp];
                    break;
                case NablaOperator::Grad:
                    if (lumped)
                    {
                        const Real length = out_of_plane_length(element_type, thickness, xyz[qp](0));
                        m.resize(1, n);
                        for (unsigned int i = 0; i < n; ++i)
                            m(0, i) = -phi[i][qp] / length;
                    }
                    else
                    {
                        m.resize(dim, n);
                        for (unsigned int a = 0; a < dim; ++a)
                            for (unsigned int i = 0; i < n; ++i)
                                m(a, i) = (*dphi)[i][qp](a);
                    }
                    break;
                case NablaOperator::Curl:
                {
                    // curl(f e_z) = (df/dy, -df/dx); axisymmetric: curl(f e_theta) = (-df/dz, df/dr + f/r).
                    const Real sign = axisymmetric ? -1.0 : 1.0;
                    const Real inv_r = axisymmetric ? 1.0 / std::abs(xyz[qp](0)) : 0.0;
                    m.resize(2, n);
                    for (unsigned int i = 0; i < n; ++i)
                    {
                        m(0, i) = sign * (*dphi)[i][qp](1);
                        m(1, i) = -sign * ((*dphi)[i][qp](0) + phi[i][qp] * inv_r);
                    }
                    break;
                }
                case NablaOperator::Div:
                    m.resize(1, n); // Out-of-plane fields are divergence free.
                    break;
                }
            }
            return out;
        }

        Nabla::QuadratureMatrices vector_operator(const ElementTypeInstance &element_type,
                                                  const VariableSpec &spec,
                                                  NablaOperator op,
                                                  const libMesh::Elem &elem)
        {
            const unsigned int dim = element_type.dimension();
            const bool axisymmetric = element_type.options().axisymmetric;

            if (op == NablaOperator::Grad)
            {
                if (spec.family != libMesh::LAGRANGE_VEC && spec.family != libMesh::L2_LAGRANGE_VEC)
                    fail("grad of vector variable '" + spec.name + "' needs a Lagrange vector family.");
                if (axisymmetric)
                    fail("grad of vector variable '" + spec.name + "' is not supported on axisymmetric elements.");
            }
            if (op == NablaOperator::Curl && dim != 3)
                fail("curl of vector variable '" + spec.name + "' is only supported on 3D elements.");

            libMesh::FEVectorBase &fe = element_type.vector_fe(spec.name);
            const auto &phi = fe.get_phi();
            const auto *dphi = (op == NablaOperator::Grad) ? &fe.get_dphi() : nullptr;
            const auto *curl_phi = (op == NablaOperator::Curl) ? &fe.get_curl_phi() : nullptr;
            const auto *div_phi = (op == NablaOperator::Div) ? &fe.get_div_phi() : nullptr;
            const auto &xyz = fe.get_xyz();
            fe.reinit(&elem);

            const unsigned int nq = fe.n_quadrature_points();
            const unsigned int n = static_cast<unsigned int>(phi.size());
            const unsigned int n_components = dim; // A vector variable always carries every component.

            // Columns stay in libMesh's order: phi_x_1 phi_y_1 phi_z_1 ... phi_x_n phi_y_n phi_z_n.
            Nabla::QuadratureMatrices out(nq);

            for (unsigned int qp = 0; qp < nq; ++qp)
            {
                DenseMatrix<Real> &m = out[qp];
                switch (op)
                {
                case NablaOperator::Identity:
                    m.resize(n_components, n);
                    for (unsigned int c = 0; c < n_components; ++c)
                        for (unsigned int i = 0; i < n; ++i)
                            m(c, i) = phi[i][qp](c);
                    break;
                case NablaOperator::Grad:
                    m.resize(n_components * dim, n);
                    for (unsigned int c = 0; c < n_components; ++c)
                        for (unsigned int b = 0; b < dim; ++b)
                            for (unsigned int i = 0; i < n; ++i)
                                m(c * dim + b, i) = (*dphi)[i][qp](c, b);
                    break;
                case NablaOperator::Curl:
                    m.resize(3, n);
                    for (unsigned int c = 0; c < 3; ++c)
                        for (unsigned int i = 0; i < n; ++i)
                            m(c, i) = (*curl_phi)[i][qp](c);
                    break;
                case NablaOperator::Div:
                {
                    // Axisymmetric (V_r, V_z): div V = dV_r/dr + V_r/r + dV_z/dz.
                    const Real inv_r = axisymmetric ? 1.0 / std::abs(xyz[qp](0)) : 0.0;
                    m.resize(1, n);
                    for (unsigned int i = 0; i < n; ++i)
                        m(0, i) = (*div_phi)[i][qp] + phi[i][qp](0) * inv_r;
                    break;
                }
                }
            }
            return out;
        }

        std::vector<unsigned int> row_components(const VariableSpec &spec,
                                                 NablaOperator op,
                                                 unsigned int dimension,
                                                 unsigned int rows)
        {
            switch (op)
            {
            case NablaOperator::Identity:
                if (spec.vector_fe)
                    break;
                return spec.active_components;
            case NablaOperator::Div:
                return {0};
            case NablaOperator::Grad:
                if (!spec.vector_fe && is_lumped_potential(spec, dimension))
                    return {2};
                break;
            case NablaOperator::Curl:
                break;
            }
            std::vector<unsigned int> components(rows);
            std::iota(components.begin(), components.end(), 0u);
            return components;
        }

    } // namespace

    OperatorExpression OperatorExpression::parse(const std::string &text)
    {
        std::string rest = strip_spaces(text);

        OperatorExpression out;
        bool have_spatial = false;
        bool have_dt = false;

        for (std::size_t open = rest.find('('); open != std::string::npos; open = rest.find('('))
        {
            if (rest.back() != ')')
                fail("malformed operator expression '" + text + "'.");
            const std::string name = rest.substr(0, open);
            std::string inner = rest.substr(open + 1, rest.size() - open - 2);

            if (name == "dt")
            {
                if (have_dt)
                    fail("repeated dt() in '" + text + "'.");
                have_dt = true;
                out.time_derivative = true;
            }
            else if (name == "grad" || name == "curl" || name == "div")
            {
                if (have_spatial)
                    fail("more than one spatial operator in '" + text + "'.");
                have_spatial = true;
                out.op = name == "grad" ? NablaOperator::Grad : name == "curl" ? NablaOperator::Curl
                                                                               : NablaOperator::Div;
            }
            else
                fail("unknown operator '" + name + "' in '" + text + "'.");
            rest = std::move(inner);
        }

        if (!is_identifier(rest))
            fail("malformed operator expression '" + text + "'.");
        out.variable = rest;
        return out;
    }

    std::vector<ExpressionTerm> ExpressionTerm::parse_sum(const std::string &text,
                                                          const std::vector<std::string> &coefficient_names)
    {
        const std::string s = strip_spaces(text);
        if (s.empty())
            fail("empty expression.");

        struct RawTerm
        {
            Real sign = 1.0;
            std::string body;
        };
        std::vector<RawTerm> raw_terms;
        RawTerm current;
        int depth = 0;
        for (char c : s)
        {
            if (c == '(')
                ++depth;
            else if (c == ')')
                --depth;
            if (depth == 0 && (c == '+' || c == '-'))
            {
                if (!current.body.empty())
                {
                    raw_terms.push_back(current);
                    current = RawTerm();
                }
                current.sign *= (c == '-') ? -1.0 : 1.0;
                continue;
            }
            current.body.push_back(c);
        }
        if (depth != 0 || current.body.empty())
            fail("malformed expression '" + text + "'.");
        raw_terms.push_back(current);

        std::vector<ExpressionTerm> terms;
        for (const RawTerm &raw : raw_terms)
        {
            std::vector<std::string> factors;
            std::string factor;
            depth = 0;
            for (char c : raw.body)
            {
                if (c == '(')
                    ++depth;
                else if (c == ')')
                    --depth;
                if (depth == 0 && c == '*')
                {
                    factors.push_back(factor);
                    factor.clear();
                    continue;
                }
                factor.push_back(c);
            }
            factors.push_back(factor);

            ExpressionTerm term;
            term.sign = raw.sign;
            bool have_operand = false;
            for (const std::string &f : factors)
            {
                if (std::find(coefficient_names.begin(), coefficient_names.end(), f) != coefficient_names.end())
                {
                    if (!term.coefficient.empty())
                        fail("more than one coefficient in a term of '" + text + "'.");
                    term.coefficient = f;
                }
                else
                {
                    if (have_operand)
                        fail("more than one operator expression in a term of '" + text + "'.");
                    term.operand = OperatorExpression::parse(f);
                    have_operand = true;
                }
            }
            if (!have_operand)
                fail("term without an operator expression in '" + text + "'.");
            terms.push_back(std::move(term));
        }
        return terms;
    }

    Nabla::Operand Nabla::operand(const ElementTypeInstance &element_type,
                                  const std::string &variable,
                                  NablaOperator op,
                                  const libMesh::Elem &elem,
                                  Real thickness)
    {
        const VariableSpec &spec = find_spec(element_type, variable);

        Operand out;
        out.values = spec.vector_fe ? vector_operator(element_type, spec, op, elem)
                                    : scalar_operator(element_type, spec, op, elem, thickness);
        const unsigned int rows = out.values.empty() ? 0 : out.values[0].m();
        out.components = row_components(spec, op, element_type.dimension(), rows);
        return out;
    }

    Nabla::QuadratureMatrices Nabla::evaluate(const ElementTypeInstance &element_type,
                                              const std::string &variable,
                                              NablaOperator op,
                                              const libMesh::Elem &elem,
                                              Real thickness)
    {
        return operand(element_type, variable, op, elem, thickness).values;
    }

    Nabla::QuadratureMatrices Nabla::evaluate(const ElementTypeInstance &element_type,
                                              const OperatorExpression &expression,
                                              const libMesh::Elem &elem,
                                              Real thickness)
    {
        return evaluate(element_type, expression.variable, expression.op, elem, thickness);
    }

    std::vector<libMesh::Real> Nabla::integration_weights(const ElementTypeInstance &element_type,
                                                          const libMesh::Elem &elem,
                                                          Real thickness)
    {
        libMesh::FEBase &fe = element_type.geometry_fe();
        const auto &JxW = fe.get_JxW();
        const auto &xyz = fe.get_xyz();
        fe.reinit(&elem);

        std::vector<Real> weights(JxW.begin(), JxW.end());
        if (element_type.dimension() == 2)
            for (std::size_t qp = 0; qp < weights.size(); ++qp)
                weights[qp] *= out_of_plane_length(element_type, thickness, xyz[qp](0));
        return weights;
    }
} // namespace pytec3d
