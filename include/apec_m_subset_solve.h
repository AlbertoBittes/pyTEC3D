#ifndef APEC_M_SUBSET_SOLVE_H
#define APEC_M_SUBSET_SOLVE_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include <libmesh/system_subset_by_subdomain.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include "apec_m_system_objects.h" // <-- Add this line

#include "fem_service_abstract.h"
#include "time_step.h"

using namespace libMesh;
// Local Includes
#include "libmesh/system_subset.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <cstddef>
#include <set>

// Forward Declarations

class ApecMSubset : public SystemSubset,
                    public ParallelObject
{
public:
    ApecMSubset(const System &system) : SystemSubset(system),
                                        ParallelObject(system),
                                        _dof_ids()
    {
        this->init();
    };

    virtual ~ApecMSubset() = default;

    inline virtual const std::vector<unsigned int> &dof_ids() const override
    {
        return _dof_ids;
    };

    void init()
    {
        _dof_ids.clear();

        std::vector<std::vector<dof_id_type>> dof_ids_per_processor(this->n_processors());
        std::vector<std::vector<dof_id_type>> constrained_dof_ids_per_processor(this->n_processors());

        const DofMap &dof_map = _system.get_dof_map();
        // std::vector<dof_id_type> dof_indices;

        for (auto it = dof_map.constraint_rows_begin(); it != dof_map.constraint_rows_end(); ++it)
        {
            // 'it' is a DofConstraints::const_iterator
            // If _dof_constraints is a std::map, then:
            // *it is a std::pair<const dof_id_type, DofConstraintRow>
            dof_id_type constrained_dof = it->first;  // key
            const DofConstraintRow &row = it->second; // value

            for (auto proc : make_range(this->n_processors()))
                if ((constrained_dof >= dof_map.first_dof(proc)) && (constrained_dof < dof_map.end_dof(proc)))
                    constrained_dof_ids_per_processor[proc].push_back(constrained_dof);
        }

        // For each processor, compute unconstrained DOFs
        for (auto proc : make_range(this->n_processors()))
        {
            dof_id_type first = dof_map.first_dof(proc);
            dof_id_type last = dof_map.end_dof(proc);

            // Efficiently fill all_dofs with pre-sizing
            std::vector<dof_id_type> all_dofs(last - first);
            for (dof_id_type dof = first; dof < last; ++dof)
                all_dofs[dof - first] = dof;

            // Get the constrained DOFs for this processor
            std::vector<dof_id_type> &constrained_dofs = constrained_dof_ids_per_processor[proc];
            std::sort(constrained_dofs.begin(), constrained_dofs.end());

            // Compute unconstrained DOFs: all_dofs - constrained_dofs, move result directly
            std::vector<dof_id_type> unconstrained_dofs;
            std::set_difference(
                all_dofs.begin(), all_dofs.end(),
                constrained_dofs.begin(), constrained_dofs.end(),
                std::back_inserter(unconstrained_dofs));

            dof_ids_per_processor[proc] = std::move(unconstrained_dofs); // move instead of copy
            // Now unconstrained_dofs contains the indices of unconstrained DOFs for this processor
            // (You can store or use unconstrained_dofs as needed here)
            // Example: print the number of unconstrained DOFs for this processor
            // std::cout << "Processor " << proc << " has " << unconstrained_dofs.size() << " unconstrained DOFs." << std::endl;
        }
        /* Distribute information among processors.  */
        std::vector<Parallel::Request> request_per_processor(this->n_processors());
        for (auto proc : make_range(this->n_processors()))
            if (proc != this->processor_id())
                this->comm().send(proc, dof_ids_per_processor[proc], request_per_processor[proc]);
        for (auto proc : make_range(this->n_processors()))
        {
            std::vector<dof_id_type> received_dofs;
            if (proc == this->processor_id())
                received_dofs = dof_ids_per_processor[proc];
            else
                this->comm().receive(proc, received_dofs);
            _dof_ids.insert(_dof_ids.end(), received_dofs.begin(), received_dofs.end());
        }

        /* Sort and unique the vector (using the same mechanism as in \p
           DofMap::prepare_send_list()).  */
        std::sort(_dof_ids.begin(), _dof_ids.end());
        std::vector<unsigned int>::iterator new_end = std::unique(_dof_ids.begin(), _dof_ids.end());
        std::vector<unsigned int>(_dof_ids.begin(), new_end).swap(_dof_ids);

        /* Wait for sends to be complete.  */
        for (auto proc : make_range(this->n_processors()))
        {
            if (proc != this->processor_id())
            {
                request_per_processor[proc].wait();
            }
        }
    }

protected:
    std::vector<unsigned int> _dof_ids;

}; // class SystemSubset
#endif // APEC_M_SUBSET_SOLVE_H