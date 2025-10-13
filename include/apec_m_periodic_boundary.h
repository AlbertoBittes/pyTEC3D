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
                          const std::vector<std::string> &vars = {},
                          EquationSystems *eq_sys = nullptr
                          //   const std::vector<std::vector<Real>> &vars_trans_mtx,
                          )
        : PeriodicBoundaryBase(),
          FEParentObjectBase<EquationSystems>(eq_sys)
    //   periodic_matrix(periodic_mtx),
    //   variables_transformation_matrix
    {

        rotation_matrix = TensorValue<Real>(periodic_mtx[0][0], periodic_mtx[0][1], periodic_mtx[0][2],
                                            periodic_mtx[1][0], periodic_mtx[1][1], periodic_mtx[1][2],
                                            periodic_mtx[2][0], periodic_mtx[2][1], periodic_mtx[2][2]);

        translation_vector = VectorValue<Real>(periodic_mtx[0][3], periodic_mtx[1][3], periodic_mtx[2][3]);

        myboundary = main_b;
        pairedboundary = paired_b;

        if (t == INVERSE)
        {
            std::swap(myboundary, pairedboundary);
            rotation_matrix = rotation_matrix.transpose();
            translation_vector = -rotation_matrix * translation_vector;
        }

        if (vars.size() >= 0)
        {
            this->set_variables(vars);
        }

        if (variables.size() == 3)
            this->set_transformation_matrix(DenseMatrix<Real>(3, 3, {rotation_matrix(0, 0), rotation_matrix(0, 1), rotation_matrix(0, 2), rotation_matrix(1, 0), rotation_matrix(1, 1), rotation_matrix(1, 2), rotation_matrix(2, 0), rotation_matrix(2, 1), rotation_matrix(2, 2)}));
    }

    ApecMPeriodicBoundary(const ApecMPeriodicBoundary &o, TransformationType t = FORWARD)
        : ApecMPeriodicBoundary(o.myboundary, o.pairedboundary,
                                {{o.rotation_matrix(0, 0), o.rotation_matrix(0, 1), o.rotation_matrix(0, 2), o.translation_vector(0)},
                                 {o.rotation_matrix(1, 0), o.rotation_matrix(1, 1), o.rotation_matrix(1, 2), o.translation_vector(1)},
                                 {o.rotation_matrix(2, 0), o.rotation_matrix(2, 1), o.rotation_matrix(2, 2), o.translation_vector(2)},
                                 {0, 0, 0, 1}},
                                t,
                                {}, // vars
                                o.get_parent_object())
    {
    }
    // // Getters (optional, for testing)
    virtual Point get_corresponding_pos(const Point &pt) const override
    {
        // Note that since _theta defines the angle from "paired boundary" to
        // "my boundary", and we want the inverse of that here, we multiply
        // by the transpose of the transformation matrix below.
        Point pt_prime;
        pt_prime = rotation_matrix * pt + translation_vector;
        return pt_prime;
    }

    void set_variables(const std::vector<std::string> &vars)
    {
        variables.clear();
        for (const auto &var_str : vars)
        {
            unsigned int var_index = this->get_parent_object()->get_system(0).variable_number(var_str);
            set_variable(var_index);
        }
    }

    virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override
    {
        return std::make_unique<ApecMPeriodicBoundary>(*this, t);
    }

private:
    TensorValue<Real> rotation_matrix;
    VectorValue<Real> translation_vector;
};

#endif // APEC_M_PERIODIC_BOUNDARY_H