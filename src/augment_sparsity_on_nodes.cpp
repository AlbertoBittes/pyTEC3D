#include "augment_sparsity_on_nodes.h"
#include "apec_m_mesh.h"

#include <algorithm>
#include <unordered_set>

#include <libmesh/elem.h>

using namespace libMesh;

namespace
{
    bool add_pair(std::vector<std::pair<dof_id_type, dof_id_type>> &pairs,
                  dof_id_type first,
                  dof_id_type second)
    {
        const auto candidate = std::make_pair(first, second);
        if (std::find(pairs.begin(), pairs.end(), candidate) != pairs.end())
            return false;
        pairs.push_back(candidate);
        return true;
    }

}

AugmentSparsityOnNodes::AugmentSparsityOnNodes(ApecMMesh &mesh)
    : _mesh(mesh)
{
    build_node_sets();
    rebuild_elem_neighbors();
}

void AugmentSparsityOnNodes::build_node_sets()
{
    node_pairs.clear();
    const auto &pairs_by_cp = _mesh.get_cp_pairs();
    for (const auto &entry : pairs_by_cp)
        node_pairs.insert(node_pairs.end(), entry.second.begin(), entry.second.end());
}

void AugmentSparsityOnNodes::rebuild_elem_neighbors()
{
    elem_pairs.clear();
    slave_to_master_elem.clear();

    std::unordered_map<unsigned int, dof_id_type> node_to_elem;
    for (const auto *elem : _mesh.element_ptr_range())
    {
        if (!elem)
            continue;
        for (unsigned int i = 0; i < elem->n_nodes(); ++i)
            node_to_elem.emplace(static_cast<unsigned int>(elem->node_id(i)), elem->id());
    }

    for (const auto &node_pair : node_pairs)
    {
        const auto master_elem_it = node_to_elem.find(node_pair.first);
        const auto slave_elem_it = node_to_elem.find(node_pair.second);
        if (master_elem_it == node_to_elem.end() || slave_elem_it == node_to_elem.end())
            continue;

        if (!add_pair(elem_pairs, master_elem_it->second, slave_elem_it->second))
            continue;

        slave_to_master_elem.emplace(slave_elem_it->second, master_elem_it->second);
    }
}

void AugmentSparsityOnNodes::operator()(const MeshBase::const_element_iterator &range_begin,
                                        const MeshBase::const_element_iterator &range_end,
                                        processor_id_type p,
                                        map_type &coupled_elements)
{
    const CouplingMatrix *const null_mat = nullptr;

    for (const auto &elem : as_range(range_begin, range_end))
    {
        if (!elem)
            continue;

        if (elem->processor_id() != p)
            continue;

        const auto master_it = slave_to_master_elem.find(elem->id());
        if (master_it == slave_to_master_elem.end())
            continue;

        const Elem *master_elem = _mesh.elem_ptr(master_it->second);
        if (master_elem && master_elem->processor_id() != p)
            coupled_elements.emplace(master_elem, null_mat);
    }
}

void AugmentSparsityOnNodes::mesh_reinit()
{
    rebuild_elem_neighbors();
}

void AugmentSparsityOnNodes::redistribute()
{
    this->mesh_reinit();
}
