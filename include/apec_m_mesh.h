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

    static std::shared_ptr<ApecMMesh> read_mesh(libMesh::Parallel::Communicator &comm,
                                                FEMDataRepository &mesh_repository,
                                                const std::string &mesh_filepath);

private:
    std::vector<std::string> collect_node_integer_arrays(const std::string &prefix) const;
};

#endif // APEC_M_MESH_H
