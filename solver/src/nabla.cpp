/**
 * @file nabla.cpp
 * @brief Nabla: the operators on the scalar and vector FEs.
 */
#include "nabla.h"

#include <cmath>
#include <string>
#include <type_traits>

#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using formulation::Operator;
        using libMesh::Real;
        using libMesh::RealTensorValue;
        using libMesh::RealVectorValue;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Nabla: " + message);
        }

        template <typename T>
        QpDofTable<T> &shaped(ValueTable &out, formulation::ValueType type, unsigned int nq, std::size_t n)
        {
            out.type = type;
            QpDofTable<T> &table = out.get<T>();
            table.assign(nq, std::vector<T>(n)); // Zero: some operators set one component only.
            return table;
        }

        // Radius of a quadrature point of an axisymmetric element, (x, y) = (r, z).
        Real radius(const libMesh::Point &p)
        {
            return std::sqrt(p(0) * p(0) + p(2) * p(2));
        }

        void evaluate_scalar(const formulation::NablaSpec &spec, libMesh::FEBase &fe, unsigned int /*dimension*/, bool axisymmetric,
                             const std::vector<Real> &length, ValueTable &out)
        {
            const auto &phi = fe.get_phi();
            const unsigned int nq = fe.n_quadrature_points();
            const std::size_t n = phi.size();

            switch (spec.op)
            {
            case Operator::Identity:
                if (spec.component < 0)
                {
                    auto &table = shaped<Real>(out, spec.output, nq, n);
                    for (unsigned int qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i] = phi[i][qp];
                }
                else
                {
                    auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                    for (unsigned int qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i](spec.component) = phi[i][qp];
                }
                return;
            case Operator::Grad:
            {
                if (spec.out_of_plane)
                {
                    auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                    for (unsigned int qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i](2) = -phi[i][qp] / length[qp];
                }
                else
                {
                    auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                    const auto &dphi = fe.get_dphi();
                    for (unsigned int qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i] = dphi[i][qp];
                }
                return;
            }
            case Operator::Curl:
            {
                const auto &dphi = fe.get_dphi();
                const auto &xyz = fe.get_xyz();
                const Real sign = axisymmetric ? -1.0 : 1.0;
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (unsigned int qp = 0; qp < nq; ++qp)
                {
                    const Real inv_r = axisymmetric ? 1.0 / radius(xyz[qp]) : 0.0;
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        table[qp][i](0) = sign * dphi[i][qp](1);
                        table[qp][i](1) = -sign * (dphi[i][qp](0) + phi[i][qp] * inv_r);
                    }
                }
                return;
            }
            case Operator::Div:
                break;
            }
            fail("operator not defined on a scalar FE.");
        }

        void evaluate_vector(const formulation::NablaSpec &spec, libMesh::FEVectorBase &fe, bool axisymmetric, ValueTable &out)
        {
            const auto &phi = fe.get_phi();
            const unsigned int nq = fe.n_quadrature_points();
            const std::size_t n = phi.size();

            switch (spec.op)
            {
            case Operator::Identity:
            {
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (unsigned int qp = 0; qp < nq; ++qp)
                    for (std::size_t i = 0; i < n; ++i)
                        table[qp][i] = phi[i][qp];
                return;
            }
            case Operator::Curl:
            {
                const auto &curl_phi = fe.get_curl_phi();
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (unsigned int qp = 0; qp < nq; ++qp)
                    for (std::size_t i = 0; i < n; ++i)
                        table[qp][i] = curl_phi[i][qp];
                return;
            }
            case Operator::Div:
            {
                const auto &div_phi = fe.get_div_phi();
                auto &table = shaped<Real>(out, spec.output, nq, n);
                const auto &xyz = fe.get_xyz();
                for (unsigned int qp = 0; qp < nq; ++qp)
                {
                    // Axisymmetric (V_r, V_z): div V = dV_r/dr + V_r / r + dV_z/dz.
                    const Real inv_r = axisymmetric ? 1.0 / radius(xyz[qp]) : 0.0;
                    for (std::size_t i = 0; i < n; ++i)
                        table[qp][i] = div_phi[i][qp] + phi[i][qp](0) * inv_r;
                }
                return;
            }
            case Operator::Grad:
            {
                const auto &dphi = fe.get_dphi();
                auto &table = shaped<RealTensorValue>(out, spec.output, nq, n);
                for (unsigned int qp = 0; qp < nq; ++qp)
                    for (std::size_t i = 0; i < n; ++i)
                        table[qp][i] = dphi[i][qp];
                return;
            }
            }
            fail("unknown operator.");
        }
    } // namespace

    std::size_t ValueTable::n_dofs() const
    {
        switch (type)
        {
        case formulation::ValueType::Scalar:
            return scalar.empty() ? 0 : scalar.front().size();
        case formulation::ValueType::Vector:
            return vector.empty() ? 0 : vector.front().size();
        case formulation::ValueType::Tensor:
            return tensor.empty() ? 0 : tensor.front().size();
        }
        return 0;
    }

    Nabla::Nabla(const formulation::NablaSpec &spec, libMesh::FEAbstract &fe, unsigned int dimension, bool axisymmetric)
        : spec_(&spec), fe_(&fe), dimension_(dimension), axisymmetric_(axisymmetric)
    {
        // The FE computes on reinit() only what was requested before.
        const bool vector = libMesh::FEInterface::field_type(fe.get_fe_type().family) == libMesh::TYPE_VECTOR;
        if (vector)
        {
            auto &vector_fe = static_cast<libMesh::FEVectorBase &>(fe);
            vector_fe.get_phi();
            if (spec.op == Operator::Curl)
                vector_fe.get_curl_phi();
            else if (spec.op == Operator::Div)
                vector_fe.get_div_phi();
            else if (spec.op == Operator::Grad)
            {
                if (axisymmetric)
                    fail("grad of a vector field is not defined on axisymmetric elements.");
                vector_fe.get_dphi();
            }
        }
        else
        {
            auto &scalar_fe = static_cast<libMesh::FEBase &>(fe);
            scalar_fe.get_phi();
            if (spec.op == Operator::Curl || (spec.op == Operator::Grad && !spec.out_of_plane))
                scalar_fe.get_dphi();
        }
        fe.get_xyz();
    }

    void Nabla::evaluate(const std::vector<Real> &length, ValueTable &out) const
    {
        if (libMesh::FEInterface::field_type(fe_->get_fe_type().family) == libMesh::TYPE_VECTOR)
            evaluate_vector(*spec_, static_cast<libMesh::FEVectorBase &>(*fe_), axisymmetric_, out);
        else
            evaluate_scalar(*spec_, static_cast<libMesh::FEBase &>(*fe_), dimension_, axisymmetric_, length, out);
    }
} // namespace pytec3d
