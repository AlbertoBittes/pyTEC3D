/**
 * @file nabla.h
 * @brief Nabla: one differential operator (identity, grad, curl, div) applied to the shape functions of one FE.
 *
 * The shape functions of an FE on the element are those of its MFEM finite element in physical space (FEState,
 * filled by ElementType::reinit()):
 *  - scalar FE (Lagrange):          phi Real, dphi RealGradient;
 *  - vector Lagrange (vdim copies): shape k = c n + j (c = k / n, j = k % n, n the scalar ones) is phi_j e_c;
 *  - Nedelec:                       phi RealVectorValue, curl_phi RealVectorValue (2D: its z component).
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
 *   vector  curl      Vector  curl_phi (3D); 2D: its out-of-plane component (z), axisymmetric the theta one,
 *                             dV_r/dz - dV_z/dr (minus the z one)
 *   vector  div       Real    div_phi; axisymmetric adds phi_r / r (vector Lagrange)
 *   vector  grad      Tensor  dphi, row c the gradient of component c (vector Lagrange, not axisymmetric)
 */
#pragma once

#include <vector>

#include "formulation_spec.h"
#include "la.h"

namespace pytec3d
{
    using la::Real;
    using la::RealTensorValue;
    using la::RealVectorValue;

    /** @brief Values [quadrature point][shape function]. */
    template <typename T>
    using QpDofTable = std::vector<std::vector<T>>;

    /** @brief A table of one of the three value types. */
    struct ValueTable
    {
        formulation::ValueType type = formulation::ValueType::Scalar;
        QpDofTable<Real> scalar;
        QpDofTable<RealVectorValue> vector;
        QpDofTable<RealTensorValue> tensor;

        template <typename T>
        const QpDofTable<T> &get() const;
        template <typename T>
        QpDofTable<T> &get();

        /** @brief Number of shape functions (columns). */
        std::size_t n_dofs() const;
    };

    template <>
    inline const QpDofTable<Real> &ValueTable::get() const { return scalar; }
    template <>
    inline const QpDofTable<RealVectorValue> &ValueTable::get() const { return vector; }
    template <>
    inline const QpDofTable<RealTensorValue> &ValueTable::get() const { return tensor; }
    template <>
    inline QpDofTable<Real> &ValueTable::get() { return scalar; }
    template <>
    inline QpDofTable<RealVectorValue> &ValueTable::get() { return vector; }
    template <>
    inline QpDofTable<RealTensorValue> &ValueTable::get() { return tensor; }

    /** @brief The shape functions of one FE on the current element, at its quadrature points. */
    struct FEState
    {
        formulation::FEFamily family = formulation::FEFamily::LAGRANGE;
        unsigned int vdim = 1;        ///< Copies of the scalar element (vector Lagrange), else 1.
        bool need_dphi = false;       ///< The nablas read the gradients (scalar elements).
        bool need_curl = false;       ///< The nablas read the curls (Nedelec).
        QpDofTable<Real> phi;               ///< Scalar elements.
        QpDofTable<RealVectorValue> dphi;   ///< Scalar elements.
        QpDofTable<RealVectorValue> vphi;   ///< Nedelec.
        QpDofTable<RealVectorValue> curl;   ///< Nedelec; (0, 0, curl) on 2D elements.

        bool nedelec() const { return family == formulation::FEFamily::NEDELEC_ONE; }
    };

    class Nabla
    {
    public:
        Nabla(const formulation::NablaSpec &spec, unsigned int dimension, bool axisymmetric);

        const formulation::NablaSpec &spec() const { return *spec_; }
        /** @brief What the operator reads from its FE (FEState::need_dphi, need_curl). */
        void request(FEState &fe) const;

        /**
         * @brief The table on the element of @p fe.
         * @param xyz The quadrature points (the radius of axisymmetric elements).
         * @param length Out-of-plane length at each quadrature point (2D elements; 1 on 3D ones).
         */
        void evaluate(const FEState &fe, const std::vector<la::Point> &xyz, const std::vector<Real> &length,
                      ValueTable &out) const;

    private:
        const formulation::NablaSpec *spec_;
        unsigned int dimension_;
        bool axisymmetric_;
    };
} // namespace pytec3d
