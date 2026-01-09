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
    node_to_elem.clear();
    slave_to_master_elem.clear();

    for (const auto *elem : _mesh.element_ptr_range())
    {
        if (elem->type() != NODEELEM || elem->n_nodes() == 0)
            continue;
        node_to_elem.emplace(static_cast<unsigned int>(elem->node_id(0)), elem->id());
    }

    for (const auto &node_pair : node_pairs)
    {
        auto first_it = node_to_elem.find(node_pair.first);
        auto second_it = node_to_elem.find(node_pair.second);
        if (first_it == node_to_elem.end() || second_it == node_to_elem.end())
            continue;

        add_pair(elem_pairs, first_it->second, second_it->second);
        slave_to_master_elem.emplace(second_it->second, first_it->second);
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

        if (elem->processor_id() == p)
            continue;

        if (elem->type() != NODEELEM)
        {
            for (auto node_idx : elem->node_index_range())
            {
                const unsigned int node_id = static_cast<unsigned int>(elem->node_id(node_idx));
                auto node_it = node_to_elem.find(node_id);
                if (node_it == node_to_elem.end())
                    continue;

                const dof_id_type slave_elem_id = node_it->second;
                auto master_it = slave_to_master_elem.find(slave_elem_id);
                if (master_it == slave_to_master_elem.end())
                    continue;

                coupled_elements.emplace(elem, null_mat);

                const Elem *master_elem = _mesh.elem_ptr(master_it->second);
                if (master_elem && master_elem->processor_id() != p)
                    coupled_elements.emplace(master_elem, null_mat);
            }
        }

        for (const auto &pair : elem_pairs)
        {
            if (elem->id() != pair.first && elem->id() != pair.second)
                continue;

            coupled_elements.emplace(elem, null_mat);

            const dof_id_type neighbor_id = (elem->id() == pair.first) ? pair.second : pair.first;
            const Elem *neighbor = _mesh.elem_ptr(neighbor_id);
            if (neighbor && neighbor->processor_id() != p)
                coupled_elements.emplace(neighbor, null_mat);
        }
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
