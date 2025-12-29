#ifndef AUGMENT_SPARSITY_ON_NODES_H
#define AUGMENT_SPARSITY_ON_NODES_H

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>

#include "dof_couplings_factory.h"

using namespace libMesh;

class AugmentSparsityOnNodes : public GhostingFunctor
{
public:
    AugmentSparsityOnNodes(MeshBase &mesh, const DofCouplingsObject &couplings);

    void operator()(const MeshBase::const_element_iterator &range_begin,
                    const MeshBase::const_element_iterator &range_end,
                    processor_id_type p,
                    map_type &coupled_elements) override;

    void mesh_reinit() override;
    void redistribute() override;

private:
    void build_node_sets(const DofCouplingsObject &couplings);
    void rebuild_elem_neighbors();

    MeshBase &_mesh;
    std::vector<std::pair<unsigned int, unsigned int>> node_pairs;
    std::vector<std::pair<dof_id_type, dof_id_type>> elem_pairs;
};

#endif // AUGMENT_SPARSITY_ON_NODES_H
