/**
 * @file dof_space.h
 * @brief The variables of the solution systems and the DOFs of the DOF system on the MFEM mesh.
 *
 * A variable is a field label with its finite element ("A" vector Lagrange, "V" scalar Lagrange, "A" Nedelec,
 * "A_z" one component), active on the subdomains (TYPE ids) whose element types have the field. SystemLayout lists the
 * variables of a system and the size of their blocks in the element vectors: on an element, the variables active on
 * its subdomain in variable order, each the shape functions of its finite element (vector Lagrange: component by
 * component, the scalar shape functions of each).
 *
 * DofSpace adds the DOFs of the DOF system: one MFEM parallel space per variable on the whole mesh, of which the DOFs
 * of the elements where the variable is active are the unknowns; the others do not exist. The unknowns of all the
 * variables are numbered in one global index space, in blocks of consecutive indices per rank (the rows of the PETSc
 * matrices): the DOFs a rank owns in MFEM, variable by variable.
 *
 * Two levels of vectors, as in MFEM: T-vectors (mfem::Vector of the owned unknowns of this rank, true DOFs) and
 * L-vectors (all the DOFs of the variables' spaces on this rank). The prolongation P_L (HypreParMatrix, L rows,
 * T columns, one 1 per row of an unknown) gives the L-vector of a T-vector; P_L^T assembles element contributions.
 *
 * Local indices: the DOFs of the variables' spaces on this rank (owned and shared), concatenated in variable order.
 * Each local index has its global index (-1 when not an unknown) and, for a Nedelec DOF, the orientation of the DOF
 * the caller sees: the circulation along its edge from the end point with the lexicographically smaller coordinates
 * (x, then y, then z) to the other (as libMesh, tree_gauge.py). The element DOFs carry the sign that takes the local
 * shape functions of MFEM to those global DOFs.
 *
 * Each DOF of a Lagrange or first-order Nedelec space is on one point of the mesh (FEMesh::dof_points()): the point
 * arrays of the caller address the DOFs by (variable, point, component).
 */
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <mfem.hpp>
#include <petscsys.h>

#include "fe_mesh.h"
#include "formulation_spec.h"
#include "registries.h"

namespace pytec3d::fem
{
    using formulation::FEFamily;

    /** @brief A variable of a system: one field label with its finite element. */
    struct Variable
    {
        std::string label;
        registry::VariableId id = registry::VariableId::COUNT;
        FEFamily family = FEFamily::LAGRANGE;
        int order = 1;
        /** @brief Copies of the scalar element: the components of a vector Lagrange field; 1 otherwise. */
        int vdim = 1;
        /** @brief Values per point: the components of a vector field (Lagrange vdim, Nedelec the dimension), else 1. */
        unsigned int n_components = 1;
        /** @brief Physical component (0 = x, 1 = y, 2 = z) of value c of a point (a scalar: of the field it holds). */
        std::array<unsigned int, 3> components = {0, 1, 2};
        std::set<int> subdomains; ///< TYPE ids where it is active.

        bool active_on(int subdomain) const { return subdomains.count(subdomain) != 0; }
        bool edge() const { return family == FEFamily::NEDELEC_ONE; }
        bool vector_lagrange() const
        {
            return family == FEFamily::LAGRANGE_VEC || family == FEFamily::L2_LAGRANGE_VEC;
        }
        bool monomial() const { return family == FEFamily::MONOMIAL || family == FEFamily::MONOMIAL_VEC; }
    };

    /** @brief The finite elements of the families on the cells of a mesh, by order. */
    class FECatalog
    {
    public:
        FECatalog(int dimension, bool serendipity);
        ~FECatalog();
        /** @brief The collection of a Lagrange (any of the four families) or Nedelec field; aborts otherwise. */
        const mfem::FiniteElementCollection &collection(FEFamily family, int order) const;
        /** @brief Its element on @p geometry; null for the monomial families (one value per element). */
        const mfem::FiniteElement *element(FEFamily family, int order, mfem::Geometry::Type geometry) const;

    private:
        int dimension_;
        std::unique_ptr<mfem::FiniteElementCollection> h1_[3];
        std::unique_ptr<mfem::FiniteElementCollection> nd_;
    };

    /** @brief The variables of a system and their blocks in the element vectors. */
    class SystemLayout
    {
    public:
        SystemLayout(const FEMesh &mesh, const FECatalog &catalog, std::vector<Variable> variables);
        virtual ~SystemLayout() = default;

        const FEMesh &mesh() const { return *mesh_; }
        const FECatalog &catalog() const { return *catalog_; }
        unsigned int n_vars() const { return static_cast<unsigned int>(variables_.size()); }
        const Variable &variable(unsigned int v) const { return variables_[v]; }
        /** @brief Number of the variable labelled @p label; -1 if none. */
        int variable_number(const std::string &label) const;
        /** @brief Size of the block of variable @p v in the element vectors of local element @p e (0 if inactive). */
        unsigned int n_element_dofs(unsigned int v, int e) const;
        /** @brief Size of the element vectors of local element @p e. */
        unsigned int n_element_dofs(int e) const;

    protected:
        const FEMesh *mesh_;
        const FECatalog *catalog_;
        std::vector<Variable> variables_;
    };

    /** @brief The DOFs of the DOF system (see the file comment). */
    class DofSpace : public SystemLayout
    {
    public:
        DofSpace(const FEMesh &mesh, const FECatalog &catalog, std::vector<Variable> variables);
        ~DofSpace() override;

        // ---- global indices ---------------------------------------------------------------------------------------
        PetscInt n_global() const { return n_global_; }
        PetscInt n_owned() const { return n_owned_; }
        PetscInt first_owned() const { return first_owned_; }
        /** @brief Unknowns of variable @p v, all ranks. */
        PetscInt n_global(unsigned int v) const { return n_global_by_variable_[v]; }
        /** @brief Unknowns owned by each rank (all variables). */
        const std::vector<PetscInt> &owned_by_rank() const { return owned_by_rank_; }

        // ---- local indices ----------------------------------------------------------------------------------------
        std::size_t n_local() const { return global_.size(); }
        /** @brief Global index of local index @p l; -1 when it is not an unknown. */
        PetscInt global(std::size_t l) const { return global_[l]; }
        bool owned(std::size_t l) const { return owned_[l] != 0; }
        /** @brief Variable of local index @p l. */
        unsigned int variable_of(std::size_t l) const { return variable_of_[l]; }
        /** @brief Local index of DOF (MFEM vdof) @p vdof of variable @p v. */
        std::size_t local_index(unsigned int v, int vdof) const { return offsets_[v] + static_cast<std::size_t>(vdof); }

        /**
         * @brief The local index of the DOF of variable @p v on @p point, value @p c (Variable::components); -1 when
         * the point has no unknown of the variable.
         */
        long long point_dof(unsigned int v, std::int64_t point, unsigned int c = 0) const;
        /** @brief The points with a DOF of variable @p v, and the scalar DOF (MFEM ndof) of each. */
        const std::unordered_map<std::int64_t, int> &point_dofs(unsigned int v) const { return point_dofs_[v]; }

        /** @brief The element DOFs of local element @p e: local indices and signs, in element vector order. */
        void element_dofs(int e, const int *&local, const double *&sign, unsigned int &n) const;

        const mfem::ParFiniteElementSpace &space(unsigned int v) const { return *spaces_[v]; }

        /** @brief A marker per local index, OR-ed over the ranks that share each DOF (all of them get the result). */
        void reduce_or(std::vector<int> &marker) const;

        // ---- vectors ----------------------------------------------------------------------------------------------
        /** @brief The prolongation P_L of the T-vectors (rows: local indices, columns: global indices). */
        mfem::HypreParMatrix &prolongation() const { return *prolongation_; }
        /** @brief Row partition of the T-vectors (hypre's assumed partition: first and end of this rank). */
        HYPRE_BigInt *true_partition() const { return const_cast<HYPRE_BigInt *>(true_partition_); }
        /** @brief Row partition of the L-vectors over the ranks, and their total size (a block-diagonal L matrix). */
        HYPRE_BigInt *local_partition() const { return const_cast<HYPRE_BigInt *>(local_partition_); }
        HYPRE_BigInt n_local_global() const { return n_local_global_; }
        /** @brief The T-level (true DOF) matrix P_L^T A P_L of the L-level matrix @p A (finalized, n_local square). */
        std::unique_ptr<mfem::HypreParMatrix> assemble(mfem::SparseMatrix &A) const;
        /** @brief The T-vector P_L^T b of the L-vector @p b (element contributions summed over the ranks). */
        void assemble(const mfem::Vector &b, mfem::Vector &x) const { prolongation_->MultTranspose(b, x); }
        /** @brief The L-vector of the T-vector @p x (0 where not an unknown). */
        void gather(const mfem::Vector &x, std::vector<double> &local) const;
        /** @brief Values on the local indices of one element, with the element signs: in its shape functions. */
        void element_values(int e, const std::vector<double> &local, std::vector<double> &values) const;

    private:
        void build_numbering();
        void build_edge_signs();
        void build_element_dofs();

        std::vector<std::unique_ptr<mfem::ParFiniteElementSpace>> spaces_;
        std::vector<std::size_t> offsets_;
        std::vector<PetscInt> global_;
        std::vector<char> owned_;
        std::vector<unsigned int> variable_of_;
        std::vector<double> edge_sign_; ///< By local index: 1, or -1 for an edge DOF against the caller's orientation.
        std::vector<std::unordered_map<std::int64_t, int>> point_dofs_;
        std::vector<std::size_t> element_offsets_;
        std::vector<int> element_local_;
        std::vector<double> element_sign_;

        PetscInt n_global_ = 0;
        PetscInt n_owned_ = 0;
        PetscInt first_owned_ = 0;
        std::vector<PetscInt> n_global_by_variable_;
        std::vector<PetscInt> owned_by_rank_;

        std::unique_ptr<mfem::HypreParMatrix> prolongation_;
        HYPRE_BigInt true_partition_[2] = {0, 0};
        HYPRE_BigInt local_partition_[2] = {0, 0};
        HYPRE_BigInt n_local_global_ = 0;
    };
} // namespace pytec3d::fem
