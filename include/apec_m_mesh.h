#ifndef APEC_M_MESH_H
#define APEC_M_MESH_H

#include <memory>
#include <set>
#include <string>

#include <libmesh/boundary_info.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/vector_value.h>

#include "apec_m_tensor_value.h"
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

    void apply_boundary_info();
    void generate_infinite_elements(const std::shared_ptr<ElementTypesFactory> &element_types_factory);
    bool generate_nodal_couplings(const DofCouplingsObject &couplings, const std::string &mesh_name);
    AngleDataIndices detect_angle_indices() const;
    ApecMTensorValue node_ansys_rotation_matrix(dof_id_type node_id) const;
};

#endif // APEC_M_MESH_H
