/**
 * @file element_kernels.h
 * @brief Generic element matrices and post-processed fields driven by the JSON kernels.
 *
 * ElementKernels reads the kernels, coefficients and field expressions of an ElementTypeInstance and
 * evaluates them on one element with Nabla. It replaces the hand-written per-formulation builders.
 * Nabla only builds operator matrices; ElementKernels contracts test, coefficient and trial into a
 * block and adds it at the element matrix indices of the test and trial variables (so coupling blocks
 * land in the right rows and columns):
 *
 *  - build_matrices(): every kernel active in the domain adds  integral of test^T * coeff * trial  to
 *    Ke, or to Ce when the trial carries dt().
 *  - build_rhs(): kernels without a trial are sources, integral of test^T * coeff with a vector
 *    coefficient (e.g. Js[] with test A, M[] with test curl(A)); nodal loads are added directly.
 *  - build_tangent_residual(): the nonlinear (BH) version of the kernels whose coefficient is nu[][].
 *  - nodal_fields() / element_fields(): each field is a signed sum of coefficient * operator(variable)
 *    terms (element_types.json). Every operator is evaluated from the local solution and projected onto
 *    the nodal basis (L2) or averaged over the quadrature points; the coefficients are then applied
 *    to those nodal or element values and the terms are summed.
 *
 * Operator columns follow libMesh's shape-function order of each variable (for a vector variable
 * phi_x_1 phi_y_1 phi_z_1 ... phi_x_n phi_y_n phi_z_n). Where each column sits in Ke, Ce, the RHS and the
 * local solution is given by ElementParameters::dofs, taken from the DofMap, so the variable and
 * variable-group layout of the system never has to be assumed here.
 *
 * Coefficient names understood in "coeff" and in expressions:
 *  - nu[][]         diag(1 / (mu0 * mur)), rotated to global axes.
 *  - sigma[][]      diag(1 / rho), rotated; zero when a resistivity is not positive and finite.
 *  - cbrt(abs(nu))  scalar cbrt(|det nu|), the gauge penalty.
 *  - coil_av[][], coil_va[][], coil_ea[][], coil_vv, coil_ee, coil_j[][]
 *                   stranded-coil couplings built from the real constants.
 *  - M[]            permanent magnetization of the material, rotated to global axes.
 *  - any other X[]  vector read from the element load "X" (a VTK cell array), rotated to global axes.
 */
#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>

#include "element_type_instance.h"
#include "nabla.h"

class BHCurve;

namespace pytec3d
{
    struct MaterialProperties
    {
        std::array<libMesh::Real, 3> relative_permeability{1.0, 1.0, 1.0};
        std::array<libMesh::Real, 3> resistivity{0.0, 0.0, 0.0};   ///< Not positive or not finite means non-conducting.
        std::array<libMesh::Real, 3> magnetization{0.0, 0.0, 0.0}; ///< Permanent magnetization in local axes.
    };

    struct RealConstants
    {
        libMesh::Real thickness = 1.0; ///< Out-of-plane length of 2D elements.
        libMesh::Real cross_section = 1.0;
        libMesh::Real turns = 1.0;
        libMesh::Real volume = 1.0; ///< Coil volume (3D).
        libMesh::Real radius = 1.0; ///< Coil radius (axisymmetric 2D).
        libMesh::Real resistance = 1.0;
        libMesh::Real symmetry = 1.0; ///< +1 / -1 factor symmetrising the coil equations.
        std::array<libMesh::Real, 3> current_direction{0.0, 0.0, 1.0};
    };

    /** @brief Where the columns of every DOF variable's operators sit in the element matrices, vectors and local solution. */
    struct ElementDofs
    {
        /** @brief DOF variable name -> element matrix index of each shape function, in shape-function order. */
        std::map<std::string, std::vector<unsigned int>> columns;
        unsigned int n_dofs = 0; ///< Size of the element matrices.
    };

    /** @brief Everything besides the geometry that an element needs. */
    struct ElementParameters
    {
        MaterialProperties material;
        RealConstants constants;
        libMesh::TensorValue<libMesh::Real> local_to_global{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        const BHCurve *bh_curve = nullptr; ///< Set for nonlinear elements.
        ElementDofs dofs;
    };

    /** @brief Source terms read from the model for one element; vectors are in local axes. */
    struct ElementLoads
    {
        /** @brief Per-node values of a scalar DOF variable; only nodes flagged in @c inject are loaded. */
        struct Nodal
        {
            std::vector<libMesh::Real> values;
            std::vector<char> inject;
        };

        std::map<std::string, std::array<libMesh::Real, 3>> element; ///< Load name (e.g. "Js") -> vector.
        std::map<std::string, Nodal> nodal;                          ///< DOF variable name -> nodal values.
    };

    class ElementKernels
    {
    public:
        /** @brief Reads the kernels and field expressions of @p element_type, which must outlive this object. */
        explicit ElementKernels(const ElementTypeInstance &element_type);

        /** @brief Linear stiffness Ke and time-derivative matrix Ce; kernels with the BH coefficient are skipped when nonlinear. */
        void build_matrices(const libMesh::Elem &elem,
                            const ElementParameters &parameters,
                            libMesh::DenseMatrix<libMesh::Real> &Ke,
                            libMesh::DenseMatrix<libMesh::Real> &Ce) const;

        /** @brief Source vector: source kernels plus nodal loads. */
        void build_rhs(const libMesh::Elem &elem,
                       const ElementParameters &parameters,
                       const ElementLoads &loads,
                       libMesh::DenseVector<libMesh::Real> &rhs) const;

        /** @brief Names of the element loads (VTK cell arrays) the source kernels read, e.g. "Js". */
        const std::vector<std::string> &element_load_names() const { return element_load_names_; }

        /** @brief Tangent matrix and residual of the BH kernels at solution @p u; both zero when the element is linear. */
        void build_tangent_residual(const libMesh::Elem &elem,
                                    const ElementParameters &parameters,
                                    const libMesh::DenseVector<libMesh::Real> &u,
                                    libMesh::DenseMatrix<libMesh::Real> &tangent,
                                    libMesh::DenseVector<libMesh::Real> &residual) const;

        /**
         * @brief Nodal fields, one vector per nodal variable of the element type (instance order), each in the
         *        shape-function order of that variable.
         * @param u     Local DOF solution, indexed through parameters.dofs.
         * @param u_dt  Local DOF time derivative, used by dt() operands.
         */
        std::vector<libMesh::DenseVector<libMesh::Real>> nodal_fields(const libMesh::Elem &elem,
                                                                      const ElementParameters &parameters,
                                                                      const libMesh::DenseVector<libMesh::Real> &u,
                                                                      const libMesh::DenseVector<libMesh::Real> &u_dt) const;

        /** @brief Element (cell-average) fields, one vector per element variable of the element type (instance order), one entry per component. */
        std::vector<libMesh::DenseVector<libMesh::Real>> element_fields(const libMesh::Elem &elem,
                                                                        const ElementParameters &parameters,
                                                                        const libMesh::DenseVector<libMesh::Real> &u,
                                                                        const libMesh::DenseVector<libMesh::Real> &u_dt) const;

    private:
        struct Kernel
        {
            std::string name;
            bool has_trial = false; ///< False for source kernels.
            OperatorExpression trial;
            OperatorExpression test;
            std::string coefficient;
            libMesh::Real sign = 1.0;
        };

        struct Field
        {
            const VariableSpec *variable = nullptr;
            std::vector<ExpressionTerm> terms; ///< Terms whose variable is active in the domain.
        };

        const ElementTypeInstance &element_type_;
        std::vector<Kernel> kernels_; ///< Kernels active in the instance's domain.
        std::vector<std::string> element_load_names_;
        std::vector<Field> nodal_;
        std::vector<Field> element_;
    };
} // namespace pytec3d
