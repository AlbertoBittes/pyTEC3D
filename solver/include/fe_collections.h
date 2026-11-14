/**
 * @file fe_collections.h
 * @brief The second-order Lagrange spaces of the solver on the quadratic cells of the mesh, as MFEM collections.
 *
 * A second-order Lagrange field has one DOF on each point of its cell: on a cell with face-centre points (QUAD9,
 * HEX27, PRISM18) the tensor-product space of MFEM (H1_FECollection, order 2: nodal at the vertices, edge midpoints,
 * face and cell centres); on the cells without (QUAD8, HEX20, PRISM15, PYRAMID13) the serendipity spaces of the
 * points the cell has, which MFEM does not provide in 3D:
 *
 *   QUAD8      span{1, x, y, x^2, xy, y^2, x^2 y, x y^2}
 *   HEX20      the 20 monomials of superlinear degree <= 2
 *   PRISM15    P2(x, y) x {1, z}, and z^2, x z^2, y z^2
 *   PYRAMID13  the rational basis of libMesh (fe_lagrange_shape_3D.C)
 *
 * The polynomial ones are the nodal bases of their spaces (inverse Vandermonde), the same functions as libMesh's.
 * A mesh uses one kind of quadrilateral face, so that the space is conforming: QuadraticH1Collection(dim, serendipity).
 * Simplices (TRI6, TET10) are the same in both.
 */
#pragma once

#include <array>
#include <memory>

#include <mfem.hpp>

namespace pytec3d::fem
{
    /**
     * @brief A nodal second-order element on the points of a serendipity cell: QUAD8, HEX20, PRISM15 or PYRAMID13.
     * DOFs in the order of MFEM: vertices, then the edge midpoints in the edge order of the MFEM geometry.
     */
    class SerendipityElement : public mfem::NodalFiniteElement
    {
    public:
        explicit SerendipityElement(mfem::Geometry::Type geometry);

        using mfem::FiniteElement::CalcShape;
        using mfem::FiniteElement::CalcDShape;
        void CalcShape(const mfem::IntegrationPoint &ip, mfem::Vector &shape) const override;
        void CalcDShape(const mfem::IntegrationPoint &ip, mfem::DenseMatrix &dshape) const override;

    private:
        template <typename S>
        void evaluate(const std::array<S, 3> &x, S *out) const;

        mfem::Geometry::Type geometry_;
        /// Exponents of the monomials of the space (polynomial elements).
        std::vector<std::array<int, 3>> monomials_;
        /// Coefficients of the nodal basis in the monomials: shape_i = sum_j coefficients_(j, i) m_j.
        mfem::DenseMatrix coefficients_;
    };

    /**
     * @brief Second-order Lagrange on the cells of the mesh: one DOF per vertex and edge, and per quadrilateral face and
     * hexahedron unless @p serendipity.
     */
    class QuadraticH1Collection : public mfem::FiniteElementCollection
    {
    public:
        QuadraticH1Collection(int dim, bool serendipity);

        const mfem::FiniteElement *FiniteElementForGeometry(mfem::Geometry::Type geometry) const override;
        int DofForGeometry(mfem::Geometry::Type geometry) const override;
        const int *DofOrderForOrientation(mfem::Geometry::Type geometry, int orientation) const override;
        const char *Name() const override { return serendipity_ ? "pyTEC3D_H1_2_serendipity" : "pyTEC3D_H1_2"; }
        int GetContType() const override { return CONTINUOUS; }
        bool serendipity() const { return serendipity_; }

    private:
        bool serendipity_;
        mfem::H1_FECollection full_;
        std::array<std::unique_ptr<SerendipityElement>, mfem::Geometry::NUM_GEOMETRIES> serendipity_elements_;
    };
} // namespace pytec3d::fem
