#ifndef APEC_M_PERIODIC_BOUNDARY_H
#define APEC_M_PERIODIC_BOUNDARY_H

#include <libmesh/libmesh.h>

#include "libmesh/periodic_boundaries.h"
#include "libmesh/equation_systems.h"
#include "fe_parent_object.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/point.h"
#include "libmesh/dense_vector.h"
#include "globals.h"
#include <memory>
#include <set>
#include <vector>
#include <string>

// Bring in everything from the libMesh namespace
using namespace libMesh;

class ApecMPeriodicBoundary : public PeriodicBoundaryBase, FEParentObjectBase<EquationSystems>
{
public:
    ApecMPeriodicBoundary(const boundary_id_type main_b, const boundary_id_type paired_b,
                          const std::vector<std::vector<Real>> &periodic_mtx,
                          TransformationType t = FORWARD,
                          const std::vector<std::string> &variable_labels = {},
                          EquationSystems *eq_sys = nullptr);

    ApecMPeriodicBoundary(const ApecMPeriodicBoundary &o, TransformationType t = FORWARD);
    virtual Point get_corresponding_pos(const Point &pt) const override;
    void set_variables();

    virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override;

    const TensorValue<Real> &get_rotation_matrix() const;
    const VectorValue<Real> &get_translation_vector() const;
    const std::vector<std::string> &get_variable_labels() const;

private:
    DenseMatrix<Real> build_transformation_matrix() const;

    TensorValue<Real> rotation_matrix;
    VectorValue<Real> translation_vector;
    std::vector<std::string> variable_labels;
};

#endif // APEC_M_PERIODIC_BOUNDARY_H
