/**
 * @file fe_collections.cpp
 * @brief The serendipity elements and the second-order Lagrange collection (fe_collections.h).
 */
#include "fe_collections.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace pytec3d::fem
{
    namespace
    {
        /** @brief A number with its three first derivatives (forward-mode differentiation of the shape functions). */
        struct Dual
        {
            double v = 0.0;
            double d[3] = {0.0, 0.0, 0.0};

            Dual() = default;
            Dual(double value) : v(value) {}
            static Dual variable(double value, int i)
            {
                Dual x(value);
                x.d[i] = 1.0;
                return x;
            }
        };
        inline Dual operator+(const Dual &a, const Dual &b)
        {
            Dual r(a.v + b.v);
            for (int i = 0; i < 3; ++i)
                r.d[i] = a.d[i] + b.d[i];
            return r;
        }
        inline Dual operator-(const Dual &a, const Dual &b)
        {
            Dual r(a.v - b.v);
            for (int i = 0; i < 3; ++i)
                r.d[i] = a.d[i] - b.d[i];
            return r;
        }
        inline Dual operator-(const Dual &a)
        {
            Dual r(-a.v);
            for (int i = 0; i < 3; ++i)
                r.d[i] = -a.d[i];
            return r;
        }
        inline Dual operator*(const Dual &a, const Dual &b)
        {
            Dual r(a.v * b.v);
            for (int i = 0; i < 3; ++i)
                r.d[i] = a.d[i] * b.v + a.v * b.d[i];
            return r;
        }
        inline Dual operator/(const Dual &a, const Dual &b)
        {
            Dual r(a.v / b.v);
            for (int i = 0; i < 3; ++i)
                r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
            return r;
        }
        inline Dual operator+(double a, const Dual &b) { return Dual(a) + b; }
        inline Dual operator-(double a, const Dual &b) { return Dual(a) - b; }
        inline Dual operator*(double a, const Dual &b) { return Dual(a) * b; }
        inline Dual operator+(const Dual &a, double b) { return a + Dual(b); }
        inline Dual operator-(const Dual &a, double b) { return a - Dual(b); }
        inline Dual operator*(const Dual &a, double b) { return a * Dual(b); }
        inline Dual operator/(const Dual &a, double b) { return a / Dual(b); }

        inline double value(double x) { return x; }
        inline double value(const Dual &x) { return x.v; }

        template <typename S>
        S power(const S &x, int n)
        {
            S r(1.0);
            for (int i = 0; i < n; ++i)
                r = r * x;
            return r;
        }

        /** @brief The space of a polynomial serendipity element: exponents (x, y, z) of its monomials. */
        std::vector<std::array<int, 3>> monomials_of(mfem::Geometry::Type geometry)
        {
            switch (geometry)
            {
            case mfem::Geometry::SQUARE:
                return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {2, 0, 0}, {1, 1, 0}, {0, 2, 0}, {2, 1, 0}, {1, 2, 0}};
            case mfem::Geometry::CUBE:
                // Superlinear degree <= 2 (the degree ignoring the variables of exponent 1): 1, the linear and
                // quadratic monomials, x y z, the six x^2 y and the three x^2 y z.
                return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2},
                        {1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0},
                        {0, 2, 1}, {1, 0, 2}, {0, 1, 2}, {2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
            case mfem::Geometry::PRISM:
                return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {2, 0, 0}, {1, 1, 0}, {0, 2, 0},
                        {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {2, 0, 1}, {1, 1, 1}, {0, 2, 1},
                        {0, 0, 2}, {1, 0, 2}, {0, 1, 2}};
            default:
                return {};
            }
        }

        int n_dofs(mfem::Geometry::Type geometry)
        {
            switch (geometry)
            {
            case mfem::Geometry::SQUARE: return 8;
            case mfem::Geometry::CUBE: return 20;
            case mfem::Geometry::PRISM: return 15;
            case mfem::Geometry::PYRAMID: return 13;
            default:
                throw std::runtime_error("SerendipityElement: no serendipity element on geometry " + std::to_string(geometry));
            }
        }

        /** @brief The libMesh PYRAMID13 basis at (xi, eta, zeta) of the libMesh reference pyramid, base [-1, 1]^2. */
        template <typename S>
        void pyramid13(const S &xi, const S &eta, const S &zeta_in, S *out)
        {
            // As libMesh: zeta kept below the apex, denominators perturbed to avoid 0 / 0 there.
            const S zeta = value(zeta_in) > 1.0 - 1.e-6 ? S(1.0 - 1.e-6) : zeta_in;
            const double eps = 1.e-35;
            const S den = 1. - zeta + eps;
            out[0] = 0.25 * (-xi - eta - 1.) * ((1. - xi) * (1. - eta) - zeta + xi * eta * zeta / den);
            out[1] = 0.25 * (-eta + xi - 1.) * ((1. + xi) * (1. - eta) - zeta - xi * eta * zeta / den);
            out[2] = 0.25 * (xi + eta - 1.) * ((1. + xi) * (1. + eta) - zeta + xi * eta * zeta / den);
            out[3] = 0.25 * (eta - xi - 1.) * ((1. - xi) * (1. + eta) - zeta - xi * eta * zeta / den);
            out[4] = zeta * (2. * zeta - 1.);
            out[5] = 0.5 * (1. + xi - zeta) * (1. - xi - zeta) * (1. - eta - zeta) / den;
            out[6] = 0.5 * (1. + eta - zeta) * (1. - eta - zeta) * (1. + xi - zeta) / den;
            out[7] = 0.5 * (1. + xi - zeta) * (1. - xi - zeta) * (1. + eta - zeta) / den;
            out[8] = 0.5 * (1. + eta - zeta) * (1. - eta - zeta) * (1. - xi - zeta) / den;
            out[9] = zeta * (1. - xi - zeta) * (1. - eta - zeta) / den;
            out[10] = zeta * (1. + xi - zeta) * (1. - eta - zeta) / den;
            out[11] = zeta * (1. + eta - zeta) * (1. + xi - zeta) / den;
            out[12] = zeta * (1. - xi - zeta) * (1. + eta - zeta) / den;
        }

        const int *edges_of(mfem::Geometry::Type geometry, int &n_edges)
        {
            using G = mfem::Geometry;
            switch (geometry)
            {
            case G::SQUARE: n_edges = 4; return &G::Constants<G::SQUARE>::Edges[0][0];
            case G::CUBE: n_edges = 12; return &G::Constants<G::CUBE>::Edges[0][0];
            case G::PRISM: n_edges = 9; return &G::Constants<G::PRISM>::Edges[0][0];
            case G::PYRAMID: n_edges = 8; return &G::Constants<G::PYRAMID>::Edges[0][0];
            default: n_edges = 0; return nullptr;
            }
        }
    } // namespace

    SerendipityElement::SerendipityElement(mfem::Geometry::Type geometry)
        : mfem::NodalFiniteElement(mfem::Geometry::Dimension[geometry], geometry, n_dofs(geometry), 2,
                                   geometry == mfem::Geometry::PYRAMID ? mfem::FunctionSpace::Uk
                                                                        : mfem::FunctionSpace::Qk),
          geometry_(geometry)
    {
        // The points of the cell: vertices, then edge midpoints, in the order of the MFEM geometry.
        const mfem::IntegrationRule &vertices = *mfem::Geometries.GetVertices(geometry);
        const int n_vertices = vertices.GetNPoints();
        int n_edges = 0;
        const int *edges = edges_of(geometry, n_edges);
        if (n_vertices + n_edges != dof)
            throw std::runtime_error("SerendipityElement: inconsistent number of points");
        for (int v = 0; v < n_vertices; ++v)
            Nodes.IntPoint(v) = vertices.IntPoint(v);
        for (int e = 0; e < n_edges; ++e)
        {
            const mfem::IntegrationPoint &a = vertices.IntPoint(edges[2 * e]);
            const mfem::IntegrationPoint &b = vertices.IntPoint(edges[2 * e + 1]);
            Nodes.IntPoint(n_vertices + e).Set3(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
        }

        monomials_ = monomials_of(geometry);
        if (monomials_.empty())
            return;
        if (static_cast<int>(monomials_.size()) != dof)
            throw std::runtime_error("SerendipityElement: the space and the points differ in size");
        // V(i, j) = m_j(x_i); the nodal basis is V^-T m, coefficients_ = V^-1 (shape_i = sum_j V^-1(j, i) m_j).
        mfem::DenseMatrix vandermonde(dof);
        for (int i = 0; i < dof; ++i)
        {
            const mfem::IntegrationPoint &p = Nodes.IntPoint(i);
            for (int j = 0; j < dof; ++j)
                vandermonde(i, j) = std::pow(p.x, monomials_[j][0]) * std::pow(p.y, monomials_[j][1]) *
                                    std::pow(p.z, monomials_[j][2]);
        }
        mfem::DenseMatrixInverse inverse(vandermonde);
        coefficients_.SetSize(dof);
        inverse.GetInverseMatrix(coefficients_);
    }

    template <typename S>
    void SerendipityElement::evaluate(const std::array<S, 3> &x, S *out) const
    {
        if (geometry_ == mfem::Geometry::PYRAMID)
        {
            // libMesh's reference pyramid: base [-1, 1]^2 at zeta = 0, apex (0, 0, 1); MFEM's: base [0, 1]^2, apex
            // (0, 0, 1). The affine map xi = 2x + z - 1, eta = 2y + z - 1, zeta = z takes one onto the other, vertex
            // to vertex in the same order.
            pyramid13(2.0 * x[0] + x[2] - 1.0, 2.0 * x[1] + x[2] - 1.0, x[2], out);
            return;
        }
        std::vector<S> m(monomials_.size());
        for (std::size_t j = 0; j < monomials_.size(); ++j)
            m[j] = power(x[0], monomials_[j][0]) * power(x[1], monomials_[j][1]) * power(x[2], monomials_[j][2]);
        for (int i = 0; i < dof; ++i)
        {
            S s(0.0);
            for (int j = 0; j < dof; ++j)
                s = s + coefficients_(j, i) * m[static_cast<std::size_t>(j)];
            out[i] = s;
        }
    }

    void SerendipityElement::CalcShape(const mfem::IntegrationPoint &ip, mfem::Vector &shape) const
    {
        std::vector<double> out(static_cast<std::size_t>(dof));
        evaluate<double>({ip.x, dim > 1 ? ip.y : 0.0, dim > 2 ? ip.z : 0.0}, out.data());
        for (int i = 0; i < dof; ++i)
            shape(i) = out[static_cast<std::size_t>(i)];
    }

    void SerendipityElement::CalcDShape(const mfem::IntegrationPoint &ip, mfem::DenseMatrix &dshape) const
    {
        std::vector<Dual> out(static_cast<std::size_t>(dof));
        evaluate<Dual>({Dual::variable(ip.x, 0), Dual::variable(dim > 1 ? ip.y : 0.0, 1),
                        Dual::variable(dim > 2 ? ip.z : 0.0, 2)},
                       out.data());
        for (int i = 0; i < dof; ++i)
            for (int k = 0; k < dim; ++k)
                dshape(i, k) = out[static_cast<std::size_t>(i)].d[k];
    }

    QuadraticH1Collection::QuadraticH1Collection(int dim, bool serendipity)
        : mfem::FiniteElementCollection(2), serendipity_(serendipity), full_(2, dim)
    {
        if (serendipity)
            for (const auto geometry : {mfem::Geometry::SQUARE, mfem::Geometry::CUBE, mfem::Geometry::PRISM,
                                        mfem::Geometry::PYRAMID})
                if (mfem::Geometry::Dimension[geometry] <= dim)
                    serendipity_elements_[geometry] = std::make_unique<SerendipityElement>(geometry);
    }

    const mfem::FiniteElement *QuadraticH1Collection::FiniteElementForGeometry(mfem::Geometry::Type geometry) const
    {
        if (serendipity_ && serendipity_elements_[geometry])
            return serendipity_elements_[geometry].get();
        if (geometry == mfem::Geometry::PYRAMID)
            return nullptr; // A pyramid has no point at the centre of its base: only the serendipity collection.
        return full_.FiniteElementForGeometry(geometry);
    }

    int QuadraticH1Collection::DofForGeometry(mfem::Geometry::Type geometry) const
    {
        // Second order: one DOF per vertex and edge; per quadrilateral face and hexahedron in the full space only.
        switch (geometry)
        {
        case mfem::Geometry::POINT:
        case mfem::Geometry::SEGMENT:
            return 1;
        case mfem::Geometry::SQUARE:
        case mfem::Geometry::CUBE:
            return serendipity_ ? 0 : 1;
        default:
            return 0;
        }
    }

    const int *QuadraticH1Collection::DofOrderForOrientation(mfem::Geometry::Type /*geometry*/, int /*orientation*/) const
    {
        // At most one DOF per entity: the same in every orientation.
        static const int identity[1] = {0};
        return identity;
    }
} // namespace pytec3d::fem
