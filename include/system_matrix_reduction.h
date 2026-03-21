#ifndef SYSTEM_MATRIX_REDUCTION_H
#define SYSTEM_MATRIX_REDUCTION_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include <libmesh/system_subset_by_subdomain.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include "apec_m_system.h" // <-- Add this line

using namespace libMesh;
// Local Includes
#include "libmesh/system_subset.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"
#include "libmesh/dof_map.h"

// C++ includes
#include <cstddef>
#include <set>

class SystemMatrixReduction : public SystemSubset,
                              public ParallelObject
{
public:
    SystemMatrixReduction(const System &system) : SystemSubset(system),
                                                  ParallelObject(system),
                                                  _dof_ids()
    {
        this->init();
    };

    virtual ~SystemMatrixReduction() = default;

    inline virtual const std::vector<unsigned int> &dof_ids() const override
    {
        return _dof_ids;
    };

    void init()
    {
        _dof_ids.clear();

        const DofMap &dof_map = _system.get_dof_map();
        const processor_id_type proc = this->processor_id();
        const dof_id_type first = dof_map.first_dof(proc);
        const dof_id_type last = dof_map.end_dof(proc);

        std::vector<dof_id_type> constrained_dofs;
        for (auto it = dof_map.constraint_rows_begin(); it != dof_map.constraint_rows_end(); ++it)
        {
            const dof_id_type constrained_dof = it->first;
            if (constrained_dof >= first && constrained_dof < last)
                constrained_dofs.push_back(constrained_dof);
        }
        std::sort(constrained_dofs.begin(), constrained_dofs.end());
        constrained_dofs.erase(std::unique(constrained_dofs.begin(), constrained_dofs.end()),
                               constrained_dofs.end());

        std::vector<dof_id_type> all_dofs(last - first);
        for (dof_id_type dof = first; dof < last; ++dof)
            all_dofs[dof - first] = dof;

        std::vector<dof_id_type> unconstrained_dofs;
        std::set_difference(
            all_dofs.begin(), all_dofs.end(),
            constrained_dofs.begin(), constrained_dofs.end(),
            std::back_inserter(unconstrained_dofs));

        _dof_ids.assign(unconstrained_dofs.begin(), unconstrained_dofs.end());
    }

protected:
    std::vector<unsigned int> _dof_ids;

}; // class SystemSubset
#endif // SYSTEM_MATRIX_REDUCTION_H
