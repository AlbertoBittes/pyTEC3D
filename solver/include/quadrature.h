/**
 * @file quadrature.h
 * @brief The quadrature of the mesh, computed for all the elements at once, and the reference shape tables.
 *
 * QuadratureData holds, for every local element, its quadrature points X, Jacobians J and determinants at the points of
 * its rule (the Gauss rule of its geometry at the quadrature order of its element type). They are computed in batches,
 * one per (geometry, order): the reference gradients of the nodal (geometry) element at the points are tabulated once,
 * and J = sum_i x_i grad(N_i)^T is a product over flat arrays of all the elements of the batch, as MFEM's
 * GeometricFactors (which does not take mixed meshes) does.
 *
 * ReferenceTables are MFEM's DofToQuad for the elements of the fields: the shape functions of a finite element and their
 * reference gradients (Lagrange) or curls (Nedelec) at the points of a rule, tabulated once per (element, rule). The
 * physical shape functions of an element follow from them and its Jacobians (ElementType::reinit):
 *   Lagrange: phi, grad phi = J^-T G;   Nedelec (covariant): phi = J^-T phi_ref, curl phi = J curl_ref / det J (3D),
 *   curl_ref / det J (2D).
 */
#pragma once

#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include <mfem.hpp>

#include "fe_mesh.h"

namespace pytec3d::fem
{
    class QuadratureData
    {
    public:
        /** @brief The quadrature of every local element, at the order @p orders[e] (by local element). */
        QuadratureData(const FEMesh &mesh, const std::vector<int> &orders);
        /** @brief The same at the points of the rule @p rules[e] of each local element (e.g. its mesh points). */
        QuadratureData(const FEMesh &mesh, const std::vector<const mfem::IntegrationRule *> &rules);

        const mfem::IntegrationRule &rule(int e) const { return *rules_[static_cast<std::size_t>(e)]; }
        int n_points(int e) const { return rule(e).GetNPoints(); }
        /** @brief Quadrature points of element @p e: 3 values per point. */
        const double *X(int e) const { return X_.data() + 3 * offsets_[static_cast<std::size_t>(e)]; }
        /** @brief Jacobians: 9 values per point, row-major (dx_i / dxi_j, i < sdim, j < dim; zero elsewhere). */
        const double *J(int e) const { return J_.data() + 9 * offsets_[static_cast<std::size_t>(e)]; }
        const double *det(int e) const { return det_.data() + offsets_[static_cast<std::size_t>(e)]; }

    private:
        void compute(const FEMesh &mesh);

        std::vector<const mfem::IntegrationRule *> rules_;
        std::vector<std::size_t> offsets_; ///< First point of each element.
        std::vector<double> X_, J_, det_;
    };

    /** @brief The reference shape functions of a finite element at the points of a rule (MFEM's DofToQuad). */
    struct ReferenceTable
    {
        int n_points = 0;
        int n_dofs = 0;
        int dim = 0;
        std::vector<double> B;  ///< Lagrange: [q][i]; Nedelec: [q][i][dim].
        std::vector<double> G;  ///< Lagrange: reference gradients [q][i][dim]; Nedelec: curls [q][i][cdim].
        int curl_dim = 0;       ///< Nedelec: 3 in 3D, 1 in 2D.
    };

    /** @brief The reference tables, computed once per (finite element, rule). Thread-safe. */
    class ReferenceTables
    {
    public:
        const ReferenceTable &get(const mfem::FiniteElement &fe, const mfem::IntegrationRule &rule) const;

    private:
        mutable std::map<std::pair<const mfem::FiniteElement *, const mfem::IntegrationRule *>, std::unique_ptr<ReferenceTable>>
            tables_;
        mutable std::mutex mutex_;
    };
} // namespace pytec3d::fem
