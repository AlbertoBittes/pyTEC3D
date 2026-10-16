/**
 * @file nabla.h
 * @brief Nabla: one differential operator (identity, grad, curl, div) applied to the shape functions of one FE.
 *
 * A Nabla follows what libMesh gives for its FE:
 *  - scalar FE (FEBase):       phi Real, dphi RealGradient;
 *  - vector FE (FEVectorBase): phi RealVectorValue, div_phi Real, curl_phi RealVectorValue, dphi RealTensorValue.
 * Evaluated on the element the FE was last reinit()ed on, it gives a table [quadrature point][shape function]
 * of Real, RealVectorValue or RealTensorValue (NablaSpec::output):
 *
 *   FE      operator  output  value
 *   scalar  identity  Real    phi, or phi e_c for a scalar FE holding component c of a vector field
 *   scalar  grad      Vector  dphi (the in-plane gradient on 2D elements); for an out-of-plane scalar of a 2D
 *                             element (active component z: the z component of a 3D field over the out-of-plane
 *                             length L, the declared length times 2 pi r if axisymmetric) -phi / L e_z
 *   scalar  curl      Vector  curl(phi e_z) of a 2D element: (dphi/dy, -dphi/dx); axisymmetric, with (x, y) = (r, z),
 *                             curl(phi e_theta) = (-dphi/dz, dphi/dr + phi / r)
 *   vector  identity  Vector  phi
 *   vector  curl      Vector  curl_phi (3D)
 *   vector  div       Real    div_phi; axisymmetric adds phi_r / r
 *   vector  grad      Tensor  dphi (not axisymmetric)
 *
 * The shape functions are in libMesh's order: for a vector Lagrange FE, node by node and component by
 * component within a node (shape i is the scalar shape i / dim on component i % dim), the order of the DOFs of a
 * vector variable.
 */
#pragma once

#include <vector>

#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

#include "formulation_spec.h"

namespace libMesh
{
    class FEAbstract;
}

namespace pytec3d
{
    /** @brief Values [quadrature point][shape function]. */
    template <typename T>
    using QpDofTable = std::vector<std::vector<T>>;

    /** @brief A table of one of the three value types. */
    struct ValueTable
    {
        formulation::ValueType type = formulation::ValueType::Scalar;
        QpDofTable<libMesh::Real> scalar;
        QpDofTable<libMesh::RealVectorValue> vector;
        QpDofTable<libMesh::RealTensorValue> tensor;

        template <typename T>
        const QpDofTable<T> &get() const;
        template <typename T>
        QpDofTable<T> &get();

        /** @brief Number of shape functions (columns). */
        std::size_t n_dofs() const;
    };

    template <>
    inline const QpDofTable<libMesh::Real> &ValueTable::get() const { return scalar; }
    template <>
    inline const QpDofTable<libMesh::RealVectorValue> &ValueTable::get() const { return vector; }
    template <>
    inline const QpDofTable<libMesh::RealTensorValue> &ValueTable::get() const { return tensor; }
    template <>
    inline QpDofTable<libMesh::Real> &ValueTable::get() { return scalar; }
    template <>
    inline QpDofTable<libMesh::RealVectorValue> &ValueTable::get() { return vector; }
    template <>
    inline QpDofTable<libMesh::RealTensorValue> &ValueTable::get() { return tensor; }

    class Nabla
    {
    public:
        /** @brief Binds @p spec to @p fe and requests from it the data the operator reads. */
        Nabla(const formulation::NablaSpec &spec, libMesh::FEAbstract &fe, unsigned int dimension, bool axisymmetric);

        const formulation::NablaSpec &spec() const { return *spec_; }

        /**
         * @brief The table on the element the FE was last reinit()ed on.
         * @param length Out-of-plane length at each quadrature point (2D elements; 1 on 3D ones).
         */
        void evaluate(const std::vector<libMesh::Real> &length, ValueTable &out) const;

    private:
        const formulation::NablaSpec *spec_;
        libMesh::FEAbstract *fe_;
        unsigned int dimension_;
        bool axisymmetric_;
    };
} // namespace pytec3d
