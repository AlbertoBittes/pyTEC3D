/**
 * @file nabla_operators.h
 * @brief The Nabla operators as function templates, fixed at compile time by output type, dimension and axisymmetry.
 *
 *     Nabla::Operand m = Gradient<libMesh::Real, 2, false>(fe, variable);
 *
 * OutputType is libMesh::Real for a scalar variable and libMesh::RealGradient for a vector variable.
 * Row layout and conventions are those of nabla.h. The operators only read the FE: the caller must have
 * requested the data they read (get_phi(), get_dphi(), get_curl_phi(), get_div_phi(), get_xyz() when
 * axisymmetric) and reinit()ed it on the element.
 */
#pragma once

#include <cmath>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/fe.h>

#include "element_field.h"
#include "logging.h"
#include "nabla.h"

namespace pytec3d
{
    namespace nabla_detail
    {
        [[noreturn]] inline void fail(const std::string &message)
        {
            logging::error_and_abort("Nabla: " + message);
        }

        template <typename OutputType>
        constexpr bool is_vector_output = std::is_same_v<OutputType, libMesh::RealGradient>;

        template <unsigned int Dim, bool Axisymmetric>
        constexpr void check_geometry()
        {
            static_assert(Dim == 2 || Dim == 3, "Nabla operators exist for 2D and 3D elements.");
            static_assert(!Axisymmetric || Dim == 2, "Axisymmetric elements are 2D.");
        }

        // The caller requests the data and reinit()s the FE; an empty array means it did not.
        template <typename Data>
        void require_computed(const Data &data, const char *what, const Field &spec)
        {
            if (data.empty())
                fail(std::string(what) + " of '" + spec.name + "' is not computed: request it and reinit() the FE first.");
        }

        inline std::vector<unsigned int> iota_components(unsigned int n)
        {
            std::vector<unsigned int> components(n);
            std::iota(components.begin(), components.end(), 0u);
            return components;
        }
    } // namespace nabla_detail

    template <typename OutputType, unsigned int Dim, bool Axisymmetric>
    Nabla::Operand Identity(libMesh::FEGenericBase<OutputType> &fe, const Field &spec)
    {
        nabla_detail::check_geometry<Dim, Axisymmetric>();

        const auto &phi = fe.get_phi();
        nabla_detail::require_computed(phi, "shape functions", spec);

        const unsigned int nq = fe.n_quadrature_points();
        const unsigned int n = static_cast<unsigned int>(phi.size());
        Nabla::Operand out;
        out.values.resize(nq);

        for (unsigned int qp = 0; qp < nq; ++qp)
        {
            libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
            if constexpr (nabla_detail::is_vector_output<OutputType>)
            {
                m.resize(Dim, n);
                for (unsigned int c = 0; c < Dim; ++c)
                    for (unsigned int i = 0; i < n; ++i)
                        m(c, i) = phi[i][qp](c);
            }
            else
            {
                m.resize(1, n);
                for (unsigned int i = 0; i < n; ++i)
                    m(0, i) = phi[i][qp];
            }
        }

        if constexpr (nabla_detail::is_vector_output<OutputType>)
            out.components = nabla_detail::iota_components(Dim);
        else
            out.components = spec.active_components;
        return out;
    }

    template <typename OutputType, unsigned int Dim, bool Axisymmetric>
    Nabla::Operand Gradient(libMesh::FEGenericBase<OutputType> &fe, const Field &spec)
    {
        nabla_detail::check_geometry<Dim, Axisymmetric>();

        const auto &phi = fe.get_phi();
        const auto &dphi = fe.get_dphi();
        nabla_detail::require_computed(phi, "shape functions", spec);
        nabla_detail::require_computed(dphi, "shape function gradients", spec);

        const unsigned int nq = fe.n_quadrature_points();
        const unsigned int n = static_cast<unsigned int>(phi.size());
        Nabla::Operand out;
        out.values.resize(nq);

        if constexpr (nabla_detail::is_vector_output<OutputType>)
        {
            if constexpr (Axisymmetric)
                nabla_detail::fail("grad of vector variable '" + spec.name + "' is not supported on axisymmetric elements.");
            if (spec.family != libMesh::LAGRANGE_VEC && spec.family != libMesh::L2_LAGRANGE_VEC)
                nabla_detail::fail("grad of vector variable '" + spec.name + "' needs a Lagrange vector family.");

            for (unsigned int qp = 0; qp < nq; ++qp)
            {
                libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
                m.resize(Dim * Dim, n);
                for (unsigned int c = 0; c < Dim; ++c)
                    for (unsigned int b = 0; b < Dim; ++b)
                        for (unsigned int i = 0; i < n; ++i)
                            m(c * Dim + b, i) = dphi[i][qp](c, b);
            }
            out.components = nabla_detail::iota_components(Dim * Dim);
        }
        else
        {
            for (unsigned int qp = 0; qp < nq; ++qp)
            {
                libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
                m.resize(Dim, n);
                for (unsigned int a = 0; a < Dim; ++a)
                    for (unsigned int i = 0; i < n; ++i)
                        m(a, i) = dphi[i][qp](a);
            }
            out.components = nabla_detail::iota_components(Dim);
        }
        return out;
    }

    template <typename OutputType, unsigned int Dim, bool Axisymmetric>
    Nabla::Operand Curl(libMesh::FEGenericBase<OutputType> &fe, const Field &spec)
    {
        nabla_detail::check_geometry<Dim, Axisymmetric>();

        Nabla::Operand out;

        if constexpr (nabla_detail::is_vector_output<OutputType>)
        {
            if constexpr (Dim != 3)
                nabla_detail::fail("curl of vector variable '" + spec.name + "' is only supported on 3D elements.");
            else
            {
                const auto &phi = fe.get_phi();
                const auto &curl_phi = fe.get_curl_phi();
                nabla_detail::require_computed(phi, "shape functions", spec);
                nabla_detail::require_computed(curl_phi, "shape function curls", spec);

                const unsigned int nq = fe.n_quadrature_points();
                const unsigned int n = static_cast<unsigned int>(phi.size());
                out.values.resize(nq);
                for (unsigned int qp = 0; qp < nq; ++qp)
                {
                    libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
                    m.resize(3, n);
                    for (unsigned int c = 0; c < 3; ++c)
                        for (unsigned int i = 0; i < n; ++i)
                            m(c, i) = curl_phi[i][qp](c);
                }
                out.components = nabla_detail::iota_components(3);
            }
        }
        else
        {
            // The curl of a scalar is that of the out-of-plane vector f e_z of a 2D element.
            if constexpr (Dim != 2)
                nabla_detail::fail("curl of scalar variable '" + spec.name + "' needs a 2D element.");
            else
            {
                const auto &phi = fe.get_phi();
                const auto &dphi = fe.get_dphi();
                nabla_detail::require_computed(phi, "shape functions", spec);
                nabla_detail::require_computed(dphi, "shape function gradients", spec);
                const std::vector<libMesh::Point> *xyz = nullptr;
                if constexpr (Axisymmetric)
                {
                    xyz = &fe.get_xyz();
                    nabla_detail::require_computed(*xyz, "quadrature points", spec);
                }

                // curl(f e_z) = (df/dy, -df/dx); axisymmetric: curl(f e_theta) = (-df/dz, df/dr + f/r).
                constexpr libMesh::Real sign = Axisymmetric ? -1.0 : 1.0;
                const unsigned int nq = fe.n_quadrature_points();
                const unsigned int n = static_cast<unsigned int>(phi.size());
                out.values.resize(nq);
                for (unsigned int qp = 0; qp < nq; ++qp)
                {
                    libMesh::Real inv_r = 0.0;
                    if constexpr (Axisymmetric)
                        inv_r = 1.0 / std::abs((*xyz)[qp](0));

                    libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
                    m.resize(2, n);
                    for (unsigned int i = 0; i < n; ++i)
                    {
                        m(0, i) = sign * dphi[i][qp](1);
                        libMesh::Real dx = dphi[i][qp](0);
                        if constexpr (Axisymmetric)
                            dx += phi[i][qp] * inv_r;
                        m(1, i) = -sign * dx;
                    }
                }
                out.components = {0, 1};
            }
        }
        return out;
    }

    template <typename OutputType, unsigned int Dim, bool Axisymmetric>
    Nabla::Operand Divergence(libMesh::FEGenericBase<OutputType> &fe, const Field &spec)
    {
        nabla_detail::check_geometry<Dim, Axisymmetric>();

        const auto &phi = fe.get_phi();
        nabla_detail::require_computed(phi, "shape functions", spec);

        const unsigned int nq = fe.n_quadrature_points();
        const unsigned int n = static_cast<unsigned int>(phi.size());
        Nabla::Operand out;
        out.values.resize(nq);

        if constexpr (nabla_detail::is_vector_output<OutputType>)
        {
            const auto &div_phi = fe.get_div_phi();
            nabla_detail::require_computed(div_phi, "shape function divergences", spec);
            const std::vector<libMesh::Point> *xyz = nullptr;
            if constexpr (Axisymmetric)
            {
                xyz = &fe.get_xyz();
                nabla_detail::require_computed(*xyz, "quadrature points", spec);
            }

            // Axisymmetric (V_r, V_z): div V = dV_r/dr + V_r/r + dV_z/dz.
            for (unsigned int qp = 0; qp < nq; ++qp)
            {
                libMesh::Real inv_r = 0.0;
                if constexpr (Axisymmetric)
                    inv_r = 1.0 / std::abs((*xyz)[qp](0));

                libMesh::DenseMatrix<libMesh::Real> &m = out.values[qp];
                m.resize(1, n);
                for (unsigned int i = 0; i < n; ++i)
                {
                    libMesh::Real value = div_phi[i][qp];
                    if constexpr (Axisymmetric)
                        value += phi[i][qp](0) * inv_r;
                    m(0, i) = value;
                }
            }
        }
        else
        {
            // The divergence of the out-of-plane vector f e_z of a 2D element is zero.
            if constexpr (Dim != 2)
                nabla_detail::fail("div of scalar variable '" + spec.name + "' needs a 2D element.");
            else
                for (libMesh::DenseMatrix<libMesh::Real> &m : out.values)
                    m.resize(1, n);
        }
        out.components = {0};
        return out;
    }
} // namespace pytec3d
