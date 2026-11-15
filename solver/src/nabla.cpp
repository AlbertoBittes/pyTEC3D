/**
 * @file nabla.cpp
 * @brief Nabla: the operators on the scalar, vector Lagrange and Nedelec FEs.
 */
#include "nabla.h"

#include <cmath>
#include <string>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using formulation::Operator;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Nabla: " + message);
        }

        template <typename T>
        QpDofTable<T> &shaped(ValueTable &out, formulation::ValueType type, std::size_t nq, std::size_t n)
        {
            out.type = type;
            QpDofTable<T> &table = out.get<T>();
            table.assign(nq, std::vector<T>(n)); // Zero: some operators set one component only.
            return table;
        }

        // Radius of a quadrature point of an axisymmetric element, (x, y) = (r, z).
        Real radius(const la::Point &p)
        {
            return std::sqrt(p(0) * p(0) + p(2) * p(2));
        }

        void evaluate_scalar(const formulation::NablaSpec &spec, const FEState &fe, bool axisymmetric,
                             const std::vector<la::Point> &xyz, const std::vector<Real> &length, ValueTable &out)
        {
            const auto &phi = fe.phi;
            const std::size_t nq = phi.size();
            const std::size_t n = nq ? phi[0].size() : 0;
            switch (spec.op)
            {
            case Operator::Identity:
                if (spec.component < 0)
                {
                    auto &table = shaped<Real>(out, spec.output, nq, n);
                    for (std::size_t qp = 0; qp < nq; ++qp)
                        table[qp] = phi[qp];
                }
                else
                {
                    auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                    for (std::size_t qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i](static_cast<unsigned int>(spec.component)) = phi[qp][i];
                }
                return;
            case Operator::Grad:
            {
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                if (spec.out_of_plane)
                {
                    for (std::size_t qp = 0; qp < nq; ++qp)
                        for (std::size_t i = 0; i < n; ++i)
                            table[qp][i](2) = -phi[qp][i] / length[qp];
                }
                else
                    for (std::size_t qp = 0; qp < nq; ++qp)
                        table[qp] = fe.dphi[qp];
                return;
            }
            case Operator::Curl:
            {
                const Real sign = axisymmetric ? -1.0 : 1.0;
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (std::size_t qp = 0; qp < nq; ++qp)
                {
                    const Real inv_r = axisymmetric ? 1.0 / radius(xyz[qp]) : 0.0;
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        table[qp][i](0) = sign * fe.dphi[qp][i](1);
                        table[qp][i](1) = -sign * (fe.dphi[qp][i](0) + phi[qp][i] * inv_r);
                    }
                }
                return;
            }
            case Operator::Div:
            case Operator::SymGrad:
                break;
            }
            fail("operator not defined on a scalar FE.");
        }

        void evaluate_vector_lagrange(const formulation::NablaSpec &spec, const FEState &fe, unsigned int dimension,
                                      bool axisymmetric, const std::vector<la::Point> &xyz, ValueTable &out)
        {
            const std::size_t nq = fe.phi.size();
            const std::size_t n = nq ? fe.phi[0].size() : 0;
            const unsigned int vdim = fe.vdim;
            const std::size_t total = n * vdim;
            switch (spec.op)
            {
            case Operator::Identity:
            {
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, total);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    for (unsigned int c = 0; c < vdim; ++c)
                        for (std::size_t j = 0; j < n; ++j)
                            table[qp][c * n + j](c) = fe.phi[qp][j];
                return;
            }
            case Operator::Curl:
            {
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, total);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    for (unsigned int c = 0; c < vdim; ++c)
                        for (std::size_t j = 0; j < n; ++j)
                        {
                            const RealVectorValue &g = fe.dphi[qp][j];
                            RealVectorValue &curl = table[qp][c * n + j];
                            if (dimension == 3)
                            {
                                // curl(phi e_c) = grad(phi) x e_c.
                                RealVectorValue e;
                                e(c) = 1.0;
                                curl = g.cross(e);
                            }
                            else
                            {
                                // The out-of-plane component dV_y/dx - dV_x/dy; axisymmetric the theta one (minus it).
                                const Real z = c == 0 ? -g(1) : g(0);
                                curl(2) = axisymmetric ? -z : z;
                            }
                        }
                return;
            }
            case Operator::Div:
            {
                auto &table = shaped<Real>(out, spec.output, nq, total);
                for (std::size_t qp = 0; qp < nq; ++qp)
                {
                    // Axisymmetric (V_r, V_z): div V = dV_r/dr + V_r / r + dV_z/dz.
                    const Real inv_r = axisymmetric ? 1.0 / radius(xyz[qp]) : 0.0;
                    for (unsigned int c = 0; c < vdim; ++c)
                        for (std::size_t j = 0; j < n; ++j)
                            table[qp][c * n + j] = fe.dphi[qp][j](c) + (c == 0 ? fe.phi[qp][j] * inv_r : 0.0);
                }
                return;
            }
            case Operator::Grad:
            {
                auto &table = shaped<RealTensorValue>(out, spec.output, nq, total);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    for (unsigned int c = 0; c < vdim; ++c)
                        for (std::size_t j = 0; j < n; ++j)
                            for (unsigned int k = 0; k < 3; ++k)
                                table[qp][c * n + j](c, k) = fe.dphi[qp][j](k);
                return;
            }
            case Operator::SymGrad:
            {
                // (grad u + grad u^T) / 2 of the shape function phi_j e_c.
                auto &table = shaped<RealTensorValue>(out, spec.output, nq, total);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    for (unsigned int c = 0; c < vdim; ++c)
                        for (std::size_t j = 0; j < n; ++j)
                            for (unsigned int k = 0; k < 3; ++k)
                            {
                                table[qp][c * n + j](c, k) += 0.5 * fe.dphi[qp][j](k);
                                table[qp][c * n + j](k, c) += 0.5 * fe.dphi[qp][j](k);
                            }
                return;
            }
            }
            fail("unknown operator.");
        }

        void evaluate_nedelec(const formulation::NablaSpec &spec, const FEState &fe, bool axisymmetric, ValueTable &out)
        {
            const std::size_t nq = fe.vphi.size();
            const std::size_t n = nq ? fe.vphi[0].size() : 0;
            switch (spec.op)
            {
            case Operator::Identity:
            {
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    table[qp] = fe.vphi[qp];
                return;
            }
            case Operator::Curl:
            {
                // 2D: the out-of-plane component (z). Axisymmetric, with (x, y) = (r, z) and component 2 the theta
                // one: curl(V_r, V_z) . e_theta = dV_r/dz - dV_z/dr, minus the z one.
                const Real sign = axisymmetric ? -1.0 : 1.0;
                auto &table = shaped<RealVectorValue>(out, spec.output, nq, n);
                for (std::size_t qp = 0; qp < nq; ++qp)
                    for (std::size_t i = 0; i < n; ++i)
                        table[qp][i] = sign * fe.curl[qp][i];
                return;
            }
            case Operator::Div:
            case Operator::Grad:
            case Operator::SymGrad:
                break;
            }
            fail("div, grad and sym_grad are not defined on a Nedelec FE.");
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

    Nabla::Nabla(const formulation::NablaSpec &spec, unsigned int dimension, bool axisymmetric)
        : spec_(&spec), dimension_(dimension), axisymmetric_(axisymmetric)
    {
    }

    void Nabla::request(FEState &fe) const
    {
        if (fe.nedelec())
        {
            if (spec_->op == Operator::Curl)
                fe.need_curl = true;
            return;
        }
        if ((spec_->op == Operator::Grad || spec_->op == Operator::SymGrad) && axisymmetric_ && fe.vdim > 1)
            fail("grad and sym_grad of a vector field are not defined on axisymmetric elements.");
        if (spec_->op == Operator::Curl || spec_->op == Operator::Div || spec_->op == Operator::SymGrad ||
            (spec_->op == Operator::Grad && !spec_->out_of_plane))
            fe.need_dphi = true;
    }

    void Nabla::evaluate(const FEState &fe, const std::vector<la::Point> &xyz, const std::vector<Real> &length,
                         ValueTable &out) const
    {
        if (fe.nedelec())
            evaluate_nedelec(*spec_, fe, axisymmetric_, out);
        else if (fe.vdim > 1)
            evaluate_vector_lagrange(*spec_, fe, dimension_, axisymmetric_, xyz, out);
        else
            evaluate_scalar(*spec_, fe, axisymmetric_, xyz, length, out);
    }
} // namespace pytec3d
