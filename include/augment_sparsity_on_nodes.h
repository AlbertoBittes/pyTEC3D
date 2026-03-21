#ifndef AUGMENT_SPARSITY_ON_NODES_H
#define AUGMENT_SPARSITY_ON_NODES_H

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>
#include <pybind11/numpy.h>

using namespace libMesh;

class AugmentSparsityOnNodes : public GhostingFunctor
{
public:
    AugmentSparsityOnNodes(
        MeshBase &mesh,
        const pybind11::array &master_slave_pairs);

    virtual void operator()(const MeshBase::const_element_iterator &range_begin,
                            const MeshBase::const_element_iterator &range_end,
                            processor_id_type p,
                            map_type &coupled_elements) override;

    virtual void mesh_reinit() override;
    virtual void redistribute() override;

private:
    void rebuild_elem_neighbors();

    MeshBase &_mesh;
    pybind11::array cp_pairs;
    std::vector<std::pair<dof_id_type, dof_id_type>> elem_pairs;
    std::unordered_map<dof_id_type, dof_id_type> slave_to_master_elem;
};

#endif // AUGMENT_SPARSITY_ON_NODES_H
