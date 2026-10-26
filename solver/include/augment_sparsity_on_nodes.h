/**
 * @file augment_sparsity_on_nodes.h
 * @brief libMesh GhostingFunctor that couples elements sharing periodic (CP) nodes.
 *
 * Reads the point arrays of the model whose names start with @c CP_ to identify coupled-
 * node pairs, then augments the libMesh sparsity pattern so that elements sharing
 * coupled nodes are ghosted onto each other's processor partition.  This ensures
 * correct DOF-coupling constraints across periodic boundaries in distributed runs.
 */
#pragma once

#include <unordered_map>
#include <vector>

#include <libmesh/ghosting_functor.h>
#include <libmesh/mesh_base.h>

#include "model_input.h"

/**
 * @brief Ghosting functor that exposes element pairs sharing coupled (periodic) nodes.
 *
 * During mesh distribution libMesh queries this functor to determine which
 * remote elements must be ghosted onto each rank.  The functor inspects the CP_*
 * point arrays of the model to build master–slave node pairs and
 * maps them to representative element pairs.
 */
class AugmentSparsityOnNodes : public libMesh::GhostingFunctor
{
public:
    using NodePairsByCp = std::unordered_map<int, std::vector<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>>>;

    /** @brief Constructs the functor and scans the CP_<k> couplings of @p arrays (only read here). */
    AugmentSparsityOnNodes(
        libMesh::MeshBase &mesh,
        const pytec3d::ModelArrays &arrays);

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
    void rebuild_node_pairs(const pytec3d::ModelArrays &arrays);     ///< Scans the CP_* arrays and builds @c cp_pairs_by_cp.
    void rebuild_elem_neighbors(const pytec3d::ModelArrays &arrays); ///< Maps node pairs to element pairs via @c slave_to_master_elems.

    libMesh::MeshBase &_mesh;                                                                          ///< Reference to the distributed mesh.
    std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::dof_id_type>> slave_to_master_elems; ///< Maps each slave element ID to all of its master element IDs.
};
