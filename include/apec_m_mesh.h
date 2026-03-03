#ifndef APEC_M_MESH_H
#define APEC_M_MESH_H

#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>
#include <unordered_map>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <libmesh/vector_value.h>

#include "apec_m_tensor_value.h"
#include "coordinate_system.h"
#include "dof_couplings_factory.h"
#include "fe_data_object.h"
#include "fe_parent_object.h"

using namespace libMesh;

class FEMDataRepository;
class ElementTypesFactory;
class APEC_M;

class ApecMMesh : public libMesh::Mesh,
                  public FEDataObjectAbstract,
                  public FEParentObjectBase<APEC_M>
{
public:
    using libMesh::Mesh::Mesh;
    using FEParentObjectBase<APEC_M>::FEParentObjectBase;
    using ConstElemPtrRange = libMesh::SimpleRange<libMesh::MeshBase::const_element_iterator>;
    using ConstNodePtrRange = libMesh::SimpleRange<libMesh::MeshBase::const_node_iterator>;
    using ConstEntityPtrRange = std::variant<ConstElemPtrRange, ConstNodePtrRange>;

    struct AngleDataIndices
    {
        VectorValue<unsigned int> indices = VectorValue<unsigned int>(libMesh::invalid_uint,
                                                                      libMesh::invalid_uint,
                                                                      libMesh::invalid_uint);
        bool available = false;
    };

    static std::shared_ptr<ApecMMesh> read_mesh(libMesh::Parallel::Communicator &comm,
                                                FEMDataRepository &mesh_repository,
                                                const std::string &mesh_filepath);

    bool generate_nodal_couplings();
    const std::unordered_map<int, std::vector<std::pair<unsigned int, unsigned int>>> &get_cp_pairs() const;
    const std::map<dof_id_type, Real> &edge_midside_orientation_map();
    AngleDataIndices detect_angle_indices() const;
    CartesianCoordinateSystem node_coordinate_system(dof_id_type node_id) const;
    CartesianCoordinateSystem node_coordinate_system(const Node &node) const;
    void select(const std::string &integer_name,
                dof_id_type value,
                bool select_elements,
                std::vector<dof_id_type> &selected_ids,
                bool global = false) const;
    ConstEntityPtrRange entity_ptr_range(bool select_elements, bool global = false) const;

    template <typename Iterator>
    std::map<dof_id_type, Point> map_local_coordinates(Iterator begin,
                                                       Iterator end,
                                                       const CoordinateSystem &coordinate_system) const
    {
        std::map<dof_id_type, Point> result;
        for (auto it = begin; it != end; ++it)
        {
            const auto *entity_ptr = to_entity_ptr(*it);
            if (!entity_ptr)
                continue;
            const Point point = entity_point(*entity_ptr);
            result.emplace(entity_ptr->id(), coordinate_system.to_local_coordinates(point));
        }
        return result;
    }

    template <typename T>
    static const T *to_entity_ptr(const T *ptr)
    {
        return ptr;
    }

    template <typename T, typename std::enable_if_t<!std::is_pointer_v<T>, int> = 0>
    static const T *to_entity_ptr(const T &ref)
    {
        return &ref;
    }

    template <typename T>
    static Point entity_point(const T &entity)
    {
        if constexpr (std::is_base_of_v<Node, T>)
            return Point(entity(0), entity(1), entity(2));
        else
            return entity.vertex_average();
    }

private:
    void build_edge_midside_orientation_map();
    std::vector<std::string> collect_node_integer_arrays(const std::string &prefix) const;
    std::unordered_map<int, std::vector<std::pair<unsigned int, unsigned int>>> cp_pairs;
    std::map<dof_id_type, Real> edge_midside_orientation_cache;
    bool edge_midside_orientation_cache_built = false;
};

#endif // APEC_M_MESH_H
