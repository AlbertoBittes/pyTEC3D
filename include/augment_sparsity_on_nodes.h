#ifndef AUGMENT_SPARSITY_ON_NODES_H
#define AUGMENT_SPARSITY_ON_NODES_H

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>

class vtkUnstructuredGrid;

using namespace libMesh;

class AugmentSparsityOnNodes : public GhostingFunctor
{
public:
    using NodePairsByCp = std::unordered_map<int, std::vector<std::pair<dof_id_type, dof_id_type>>>;

    AugmentSparsityOnNodes(
        MeshBase &mesh,
        vtkUnstructuredGrid *vtk_grid);

    virtual void operator()(const MeshBase::const_element_iterator &range_begin,
                            const MeshBase::const_element_iterator &range_end,
                            processor_id_type p,
                            map_type &coupled_elements) override;

    virtual void mesh_reinit() override;
    virtual void redistribute() override;

public:
    NodePairsByCp cp_pairs_by_cp;

private:
    void rebuild_node_pairs();
    void rebuild_elem_neighbors();

    MeshBase &_mesh;
    vtkUnstructuredGrid *_vtk_grid = nullptr;
    std::vector<std::pair<dof_id_type, dof_id_type>> node_pairs;
    std::vector<std::pair<dof_id_type, dof_id_type>> elem_pairs;
    std::unordered_map<dof_id_type, dof_id_type> slave_to_master_elem;
};

#endif // AUGMENT_SPARSITY_ON_NODES_H
