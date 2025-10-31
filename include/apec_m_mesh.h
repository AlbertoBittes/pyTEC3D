#ifndef APEC_M_MESH_H
#define APEC_M_MESH_H

#include <memory>
#include <set>
#include <string>

#include <libmesh/boundary_info.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>

#include "fe_data_object.h"
#include "fe_parent_object.h"

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

    static std::shared_ptr<ApecMMesh> read(libMesh::Parallel::Communicator &comm,
                                           FEMDataRepository &model_data,
                                           const std::string &mesh_filepath);

    void apply_boundary_info();
    void generate_infinite_elements(const std::shared_ptr<ElementTypesFactory> &element_types_factory);
};

#endif // APEC_M_MESH_H
