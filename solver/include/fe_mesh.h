/**
 * @file fe_mesh.h
 * @brief FEMesh: the mesh of the caller (VTK cells, model_input.h) as an MFEM parallel mesh, and the maps between the
 *        points of the caller and the entities of the mesh.
 *
 * Every rank reads the whole MeshView: the cells are oriented (libMesh's rule: an element with a negative corner
 * triple product is flipped, its node order permuted), their corner points become the vertices of a serial MFEM mesh,
 * which METIS partitions into the ParMesh. The other points of the quadratic cells (edge midpoints, face and cell
 * centres) are the geometry: a second-order mesh carries its points as the nodes of a QuadraticH1Collection
 * (fe_collections.h), so its elements are curved as in the caller's mesh.
 *
 * Each point of a cell is one entity of the mesh: a corner a vertex, a midpoint an edge, a face centre a face, the
 * centre the element. A DOF of a finite element space lives on an entity, hence on a point (dof_points()): this is
 * how the point arrays of the caller (Dirichlet values, nodal loads, couplings, outputs) reach the DOFs.
 *
 * The cells of a mesh have one dimension (2 or 3) and one geometric order (all linear or all quadratic); a 2D mesh lies
 * in the xy plane. Quadrilateral faces are of one kind: with centre points (QUAD9, HEX27, PRISM18) or without (QUAD8,
 * HEX20, PRISM15, PYRAMID13, the serendipity cells).
 */
#pragma once

#include <cstdint>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include <mfem.hpp>

#include "model_input.h"

namespace pytec3d::fem
{
    /** @brief The name of a VTK cell type, as the solver reports it ("HEX20", ...). */
    const char *cell_type_name(int vtk);

    class FEMesh
    {
    public:
        /** @brief Builds the mesh; collective on @p comm. @p type_ids: the TYPE of each cell. */
        FEMesh(MPI_Comm comm, const MeshView &view, const IdArray &type_ids);
        ~FEMesh();

        MPI_Comm comm() const { return comm_; }
        int rank() const { return rank_; }
        int n_ranks() const { return n_ranks_; }

        mfem::ParMesh &mesh() { return *mesh_; }
        const mfem::ParMesh &mesh() const { return *mesh_; }
        int dimension() const { return dimension_; }
        /** @brief 1: straight cells, 2: quadratic cells, p: Lagrange cells of order p. */
        int geometric_order() const { return geometric_order_; }
        /** @brief The order of the VTK Lagrange cells of the mesh (0: none): its points are the nodes of order p. */
        int lagrange_order() const { return lagrange_order_; }
        /** @brief The quadratic cells have no face-centre points (QUAD8, HEX20, PRISM15, PYRAMID13). */
        bool serendipity() const { return serendipity_; }

        std::size_t n_points() const { return n_points_; }
        std::size_t n_cells() const { return n_cells_; }
        int n_local_elements() const { return static_cast<int>(cells_.size()); }

        /** @brief The cell (caller's index) of local element @p e. */
        std::int64_t cell(int e) const { return cells_[static_cast<std::size_t>(e)]; }
        /** @brief The TYPE id of local element @p e. */
        int subdomain(int e) const { return subdomains_[static_cast<std::size_t>(e)]; }
        /** @brief The VTK cell type of local element @p e. */
        int cell_type(int e) const { return cell_types_[static_cast<std::size_t>(e)]; }
        /** @brief The points of local element @p e, in the VTK order of its cell (after orientation). */
        const std::int64_t *element_points(int e, int &n) const;

        /** @brief The reference coordinates of the points of local element @p e, in element_points() order. */
        const std::vector<mfem::IntegrationPoint> &reference_points(int e) const;
        /** @brief The same as an integration rule (no weights): the shape tables of the fields at the points. Both are
         * built for the cell types of the mesh when it is created, so concurrent calls only read. */
        const mfem::IntegrationRule &point_rule(int e) const;

        /** @brief The point of a local vertex, edge, face or element; -1 when the entity has none (a straight cell). */
        std::int64_t vertex_point(int v) const { return vertex_points_[static_cast<std::size_t>(v)]; }
        std::int64_t edge_point(int edge) const { return edge_points_[static_cast<std::size_t>(edge)]; }
        std::int64_t face_point(int face) const;
        std::int64_t element_point(int e) const { return element_points_center_[static_cast<std::size_t>(e)]; }

        /**
         * @brief The point of each DOF of the scalar space @p fes (its ndofs, not vdofs): the point of the entity the DOF
         * is on, -1 where the entity has no point or carries several DOFs.
         */
        std::vector<std::int64_t> dof_points(const mfem::ParFiniteElementSpace &fes) const;

        /** @brief The coordinates of point @p p (3 values; z = 0 on a 2D mesh). */
        const double *point(std::int64_t p) const { return coordinates_.data() + 3 * p; }
        /** @brief The end points of the edge whose midpoint is @p p (any cell of the mesh); false if @p p is none. */
        bool edge_corners(std::int64_t p, std::int64_t &a, std::int64_t &b) const;

        /** @brief The TYPE ids of the cells (all ranks). */
        const std::set<int> &subdomain_ids() const { return subdomain_ids_; }

        /** @brief Elements flipped at import (all ranks). */
        std::size_t flipped_elements() const { return flipped_elements_; }
        /** @brief Number of elements of each VTK cell type (all ranks), by cell type. */
        const std::vector<std::pair<int, std::size_t>> &cell_type_counts() const { return cell_type_counts_; }

    private:
        MPI_Comm comm_;
        int rank_ = 0;
        int n_ranks_ = 1;
        int dimension_ = 0;
        int geometric_order_ = 1;
        int lagrange_order_ = 0;
        bool serendipity_ = false;
        std::size_t n_points_ = 0;
        std::size_t n_cells_ = 0;
        std::size_t flipped_elements_ = 0;
        std::vector<std::pair<int, std::size_t>> cell_type_counts_;

        std::unique_ptr<mfem::ParMesh> mesh_;
        std::vector<std::int64_t> cells_;
        std::vector<int> subdomains_;
        std::vector<int> cell_types_;
        std::vector<std::int64_t> point_offsets_; ///< Of each local element in point_ids_.
        std::vector<std::int64_t> point_ids_;     ///< The oriented points of the local elements.
        std::vector<std::int64_t> vertex_points_;
        std::vector<std::int64_t> edge_points_;
        std::vector<std::int64_t> face_points_;
        std::vector<std::int64_t> element_points_center_;
        std::vector<double> coordinates_;
        std::set<int> subdomain_ids_;
        mutable std::unordered_map<int, std::vector<mfem::IntegrationPoint>> reference_points_; ///< By VTK cell type.
        mutable std::unordered_map<int, std::unique_ptr<mfem::IntegrationRule>> point_rules_;     ///< By VTK cell type.
        std::unordered_map<std::int64_t, std::pair<std::int64_t, std::int64_t>> edge_corners_; ///< By midpoint.
    };
} // namespace pytec3d::fem
