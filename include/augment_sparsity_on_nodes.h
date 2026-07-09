/**
 * @file augment_sparsity_on_nodes.h
 * @brief libMesh GhostingFunctor that couples elements sharing periodic (CP) nodes.
 *
 * Reads VTK point-data arrays whose names start with @c CP_ to identify coupled-
 * node pairs, then augments the libMesh sparsity pattern so that elements sharing
 * coupled nodes are ghosted onto each other's processor partition.  This ensures
 * correct DOF-coupling constraints across periodic boundaries in distributed runs.
 */
#ifndef AUGMENT_SPARSITY_ON_NODES_H
#define AUGMENT_SPARSITY_ON_NODES_H

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>

class vtkUnstructuredGrid;

/**
 * @brief Ghosting functor that exposes element pairs sharing coupled (periodic) nodes.
 *
 * During mesh distribution libMesh queries this functor to determine which
 * remote elements must be ghosted onto each rank.  The functor inspects CP_*
 * point-data arrays from the VTK grid to build master–slave node pairs and
 * maps them to representative element pairs.
 */
class AugmentSparsityOnNodes : public libMesh::GhostingFunctor
{
public:
    struct NodePairInfo
    {
        libMesh::dof_id_type master_node = libMesh::invalid_uint;
        libMesh::dof_id_type slave_node = libMesh::invalid_uint;
        bool same_edge_orientation = true;
    };

    using NodePairsByCp = std::unordered_map<int, std::vector<NodePairInfo>>;

    /** @brief Constructs the functor and performs an initial CP-pair scan. */
    AugmentSparsityOnNodes(
        libMesh::MeshBase &mesh,
        vtkUnstructuredGrid *vtk_grid);

    /**
     * @brief Called by libMesh to extend the ghosting set for processor @p p.
     *
     * For every element in the given range that has a coupled-node master
     * element on a different processor, that master element is added to
     * @p coupled_elements so libMesh will ghost it.
     */
    virtual void operator()(const libMesh::MeshBase::const_element_iterator &range_begin,
                            const libMesh::MeshBase::const_element_iterator &range_end,
                            libMesh::processor_id_type p,
                            libMesh::GhostingFunctor::map_type &coupled_elements) override;

    /** @brief Re-scans CP pairs after the mesh has been redistributed. */
    virtual void mesh_reinit() override;

    /** @brief Satisfies the GhostingFunctor interface; delegates to mesh_reinit(). */
    virtual void redistribute() override;

    /** @brief Map from CP set ID to its master–slave node pairs. */
    NodePairsByCp cp_pairs_by_cp;

private:
    void rebuild_node_pairs();     ///< Scans VTK CP_* arrays and builds @c node_pairs.
    void rebuild_elem_neighbors(); ///< Maps node pairs to element pairs via @c slave_to_master_elem.

    libMesh::MeshBase &_mesh;                                                            ///< Reference to the distributed mesh.
    vtkUnstructuredGrid *_vtk_grid = nullptr;                                            ///< Non-owning VTK grid pointer.
    std::vector<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>> node_pairs;       ///< All master-slave node pairs.
    std::vector<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>> elem_pairs;       ///< Representative element pairs.
    std::unordered_map<libMesh::dof_id_type, libMesh::dof_id_type> slave_to_master_elem; ///< Maps each slave element ID to its master element ID.
};

#endif // AUGMENT_SPARSITY_ON_NODES_H
