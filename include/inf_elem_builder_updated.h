#ifndef INF_ELEM_BUILDER_UPDATED_H
#define INF_ELEM_BUILDER_UPDATED_H

#include "libmesh/inf_elem_builder.h"

using namespace libMesh;

class InfElemBuilderUpdated : public InfElemBuilder
{
public:
    InfElemBuilderUpdated(MeshBase &mesh) : InfElemBuilder(mesh), mesh_(mesh) {}

    // Method hiding (not overriding)
    void build_inf_elem(const libMesh::Point &origin,
                        const bool x_sym = false,
                        const bool y_sym = false,
                        const bool z_sym = false,
                        const bool be_verbose = false,
                        std::set<std::pair<libMesh::dof_id_type, unsigned int>> *inner_faces = nullptr);

    void set_infinite_elem_side_boundary_id(boundary_id_type id) { infinite_elem_side_boundary_id = id; }

protected:
    boundary_id_type infinite_elem_side_boundary_id;
    MeshBase &mesh_;
};

#endif // INF_ELEM_BUILDER_UPDATED_H