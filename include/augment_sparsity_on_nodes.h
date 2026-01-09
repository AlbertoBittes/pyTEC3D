#ifndef AUGMENT_SPARSITY_ON_NODES_H
#define AUGMENT_SPARSITY_ON_NODES_H

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>

using namespace libMesh;

class ApecMMesh;

class AugmentSparsityOnNodes : public GhostingFunctor
{
public:
    AugmentSparsityOnNodes(ApecMMesh &mesh);

    virtual void operator()(const MeshBase::const_element_iterator &range_begin,
                            const MeshBase::const_element_iterator &range_end,
                            processor_id_type p,
                            map_type &coupled_elements) override;

    virtual void mesh_reinit() override;
    virtual void redistribute() override;

private:
    void build_node_sets();
    void rebuild_elem_neighbors();

    ApecMMesh &_mesh;
    std::vector<std::pair<unsigned int, unsigned int>> node_pairs;
    std::vector<std::pair<dof_id_type, dof_id_type>> elem_pairs;
    std::unordered_map<unsigned int, dof_id_type> node_to_elem;
    std::unordered_map<dof_id_type, dof_id_type> slave_to_master_elem;
};

#endif // AUGMENT_SPARSITY_ON_NODES_H
