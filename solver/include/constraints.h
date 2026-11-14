/**
 * @file constraints.h
 * @brief The constraints of the DOF system: Dirichlet values, DOF couplings and the pins of static steps.
 *
 * A constraint is a row u_s = sum_m c_m u_m + value on a DOF s (the slave), the masters m other DOFs. Three sources,
 * merged in this order (a later one replaces an earlier one on the same DOF):
 *   -# DOF couplings (model data dof_couplings, the CP_<k> point arrays): on the points of a CP set, the DOFs of the
 *      coupled variable equal those of the point with the smallest id (the master); a vector Lagrange field through the
 *      node frames of the two points, an edge field with the sign of the orientations of the two edges (from the CP
 *      labels of their end points);
 *   -# Dirichlet values: the "_F" point arrays (NaN where free); a vector given in the node frame (node_angles):
 *      all its components fix the global ones, some components give rows on the global DOFs; an edge field has its
 *      value (the circulation along the edge, see dof_space.h) in component 2 of the array on the midside point;
 *   -# on a static step (no transient effects), the DOFs whose field has no row in K on any element around them are
 *      pinned to 0, unless a Dirichlet value is given (e.g. T fixed on a CUT_LAYER of T-phi).
 *
 * Every rank holds all the constraints (ConstraintSet), so that the rows of the masters, possibly constrained and
 * owned by other ranks, can be resolved (reduced_system.h).
 */
#pragma once

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include <petscsys.h>

#include "dof_space.h"
#include "la.h"
#include "model_data.h"
#include "model_input.h"

class ModelElementTypes;

namespace pytec3d
{
    /** @brief One constraint row: u = sum masters + value. */
    struct ConstraintRow
    {
        std::vector<std::pair<PetscInt, double>> masters; ///< Sorted by DOF.
        double value = 0.0;

        /** @brief Same masters and coefficients (the value may differ). */
        bool same_row(const ConstraintRow &other) const { return masters == other.masters; }
    };

    /** @brief All the constraints of the system, by global DOF, on every rank. */
    using ConstraintSet = std::map<PetscInt, ConstraintRow>;

    /** @brief The CP sets of the model: by CP id, the pairs (master point, slave point). */
    using CouplingPairs = std::map<int, std::vector<std::pair<std::int64_t, std::int64_t>>>;

    /** @brief The CP sets of the CP_<k> point arrays: master the smallest point of each set. */
    CouplingPairs coupling_pairs(const ModelArrays &arrays, std::size_t n_points);

    struct ConstraintInputs
    {
        const fem::FEMesh *mesh = nullptr;
        const fem::DofSpace *dofs = nullptr;
        const ModelArrays *arrays = nullptr;
        const ModelElementTypes *element_types = nullptr;
        const std::vector<DofCoupling> *couplings = nullptr;
        const CouplingPairs *pairs = nullptr;
        /** @brief Local-to-global node frame of each point; empty: identity. */
        const std::vector<la::RealTensorValue> *node_rotations = nullptr;
        bool static_step = false;
    };

    /** @brief Builds the constraints (collective) and logs their summary table. */
    ConstraintSet build_constraints(const ConstraintInputs &inputs);
} // namespace pytec3d
