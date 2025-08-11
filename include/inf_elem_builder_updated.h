#ifndef INF_ELEM_BUILDER_UPDATED_H
#define INF_ELEM_BUILDER_UPDATED_H

#include "libmesh/inf_elem_builder.h"

using namespace libMesh;

class InfElemBuilderUpdated : public InfElemBuilder
{
public:
    InfElemBuilderUpdated(MeshBase &mesh) : InfElemBuilder(mesh), mesh_ptr(&mesh) {}

    // Method hiding (not overriding)
    void build_inf_elem(const libMesh::Point &origin,
                        const bool x_sym = false,
                        const bool y_sym = false,
                        const bool z_sym = false,
                        const bool be_verbose = false,
                        std::set<std::pair<libMesh::dof_id_type, unsigned int>> *inner_faces = nullptr);

protected:
    MeshBase *mesh_ptr;
};

#endif // INF_ELEM_BUILDER_UPDATED_H