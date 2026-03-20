#include "apec_m_mesh.h"

#include "vtkio_updated.h"
#include "element_types_factory.h"

#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_serializer.h>
#include <libmesh/node.h>
#include <libmesh/point.h>
#include <libmesh/vector_value.h>

#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "fem_data_repository.h"

std::shared_ptr<ApecMMesh> ApecMMesh::read_mesh(libMesh::Parallel::Communicator &comm,
                                                FEMDataRepository &mesh_repository,
                                                const std::string &mesh_filepath)
{
    libMesh::out << "Reading mesh and importing field data..." << std::endl;

    if (mesh_filepath.empty())
        libmesh_error_msg("Mesh filepath not provided for ApecMMesh::read_mesh.");

    auto mesh = mesh_repository.get_object<ApecMMesh>();
    if (!mesh)
    {
        mesh = std::make_shared<ApecMMesh>(comm);
        mesh_repository.set_object<ApecMMesh>(mesh);
    }
    VTKIO_Updated vtkio(*mesh, mesh_filepath);
    MeshCommunication().broadcast(*mesh);
    MeshTools::Modification::orient_elements(*mesh);

    return mesh;
}

std::vector<std::string> ApecMMesh::collect_node_integer_arrays(const std::string &prefix) const
{
    std::vector<std::string> arrays;
    const unsigned int n_node_arrays = this->n_node_integers();
    arrays.reserve(n_node_arrays);
    for (unsigned int i = 0; i < n_node_arrays; ++i)
    {
        const std::string &name = this->get_node_integer_name(i);
        if (name.rfind(prefix, 0) != 0)
            continue;
        arrays.push_back(name);
    }
    return arrays;
}
