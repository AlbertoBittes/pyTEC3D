/**
 * @file system_matrix_reduction.h
 * @brief The unconstrained local DOFs of a system: the subset its linear solve is restricted to.
 */
#pragma once

#include <vector>

#include <libmesh/parallel_object.h>
#include <libmesh/system.h>
#include <libmesh/system_subset.h>

class SystemMatrixReduction : public libMesh::SystemSubset,
                              public libMesh::ParallelObject
{
public:
    explicit SystemMatrixReduction(const libMesh::System &system);

    const std::vector<unsigned int> &dof_ids() const override;

    /** @brief Collects the local DOFs that are not constrained. */
    void init();

private:
    std::vector<unsigned int> dof_ids_;
};
