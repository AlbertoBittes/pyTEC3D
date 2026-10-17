/**
 * @file system_matrix_reduction.cpp
 * @brief SystemMatrixReduction: the unconstrained local DOFs of a system.
 */
#include "system_matrix_reduction.h"

#include <algorithm>
#include <iterator>

#include <libmesh/dof_map.h>

SystemMatrixReduction::SystemMatrixReduction(const libMesh::System &system)
    : libMesh::SystemSubset(system),
      libMesh::ParallelObject(system)
{
    init();
}

const std::vector<unsigned int> &SystemMatrixReduction::dof_ids() const
{
    return dof_ids_;
}

void SystemMatrixReduction::init()
{
    dof_ids_.clear();

    const libMesh::DofMap &dof_map = _system.get_dof_map();
    const libMesh::processor_id_type proc = this->processor_id();
    const libMesh::dof_id_type first = dof_map.first_dof(proc);
    const libMesh::dof_id_type last = dof_map.end_dof(proc);

    std::vector<libMesh::dof_id_type> constrained_dofs;
    for (auto it = dof_map.constraint_rows_begin(); it != dof_map.constraint_rows_end(); ++it)
    {
        const libMesh::dof_id_type constrained_dof = it->first;
        if (constrained_dof >= first && constrained_dof < last)
            constrained_dofs.push_back(constrained_dof);
    }
    std::sort(constrained_dofs.begin(), constrained_dofs.end());
    constrained_dofs.erase(std::unique(constrained_dofs.begin(), constrained_dofs.end()), constrained_dofs.end());

    std::vector<libMesh::dof_id_type> all_dofs(last - first);
    for (libMesh::dof_id_type dof = first; dof < last; ++dof)
        all_dofs[dof - first] = dof;

    std::vector<libMesh::dof_id_type> unconstrained_dofs;
    std::set_difference(all_dofs.begin(), all_dofs.end(),
                        constrained_dofs.begin(), constrained_dofs.end(),
                        std::back_inserter(unconstrained_dofs));

    dof_ids_.assign(unconstrained_dofs.begin(), unconstrained_dofs.end());
}
