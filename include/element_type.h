/**
 * @file element_type.h
 * @brief An element type: its fields, its kernels, and the finite elements that evaluate them.
 *
 * Built from element_types.json by ElementTypeParser. It holds
 *  - the fields: DOF fields (unknowns), nodal fields (L2-projected) and element fields (quadrature averages),
 *    each with the physical regions (domains) where it exists;
 *  - the kernels: DOF kernels (blocks of K and C, and sources), nodal kernels and element kernels, each tied
 *    to one or two fields and to the domains where it applies;
 *  - the FE objects. Fields of the same family and order point to the same FE, so reinit() runs once per FE;
 *  - the operators (Nabla) of the kernels, each bound to its field and FE.
 *
 * It does not evaluate parameters (element_parameters.h), which are physics, nor assemble anything: the
 * element builder reinit()s the FEs, asks for the element matrix indices of the fields and calls the kernels.
 */
#pragma once

#include <deque>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/system.h>

#include "element_field.h"
#include "expression.h"
#include "nabla.h"
#include "parameter_value.h"

namespace pytec3d
{
    struct ElementTypeOptions
    {
        unsigned int order = 1;    ///< FE order of the Lagrange families.
        bool axisymmetric = false; ///< 2D elements only.
    };

    /** @brief An operator applied to a field, as written in a kernel: curl(A), dt(V), grad(dt(V)) ... */
    struct OperandRef
    {
        const Nabla *nabla = nullptr; ///< Owned by the element type, one per (field, operator).
        bool time_derivative = false; ///< dt(...): acts on the time derivative of the solution.

        const Field &field() const { return nabla->field(); }
        NablaOperator op() const { return nabla->op(); }
    };

    /** @brief Where the value of a parameter comes from. */
    enum class ParameterSource : unsigned char
    {
        Material,   ///< Column(s) of the material row.
        Constant,   ///< Column(s) of the real-constants row.
        Cell,       ///< Cell-data array of the grid, one tuple per element.
        Point,      ///< Point-data array of the grid; only usable as the nodal load of a DOF field.
        Expression, ///< Derived from other parameters.
    };

    /** @brief A named parameter: "nu[][]" tensor, "turns_density" scalar, "Js[]" vector, "AMPS()" array typed by the grid. */
    struct ParameterDefinition
    {
        std::string name; ///< As written, with its type suffix.
        ParameterType type = ParameterType::Scalar;
        ParameterSource source = ParameterSource::Expression;
        std::string array;                 ///< Grid array name for Cell and Point.
        std::vector<unsigned int> columns; ///< Row columns for Material and Constant.
        std::shared_ptr<const Expression> expression;
    };

    /** @brief A nonlinear law, e.g. BH: acts on the trial value of a kernel instead of multiplying it. */
    enum class LawType : unsigned char
    {
        BHCurve, ///< H(B) of the BH curve of the element material.
    };

    struct NonlinearLaw
    {
        std::string name;      ///< Used in kernels as "name()".
        LawType type = LawType::BHCurve;
        unsigned int linear = 0; ///< Index of the tensor parameter the law reduces to when the element has no curve.
    };

    enum class ParameterKind : unsigned char
    {
        One,       ///< No parameter: 1.
        Parameter, ///< A ParameterDefinition.
        Law,       ///< A NonlinearLaw.
    };

    /** @brief How a kernel uses a parameter: "-nu[][]", "BH()". */
    struct ParameterRef
    {
        ParameterKind kind = ParameterKind::One;
        unsigned int index = 0; ///< Into ElementType::parameters() or laws().
        ParameterType type = ParameterType::Scalar;
        libMesh::Real sign = 1.0;
        std::string name;
    };

    /** @brief K/C block (with a trial) or right-hand-side source (without) on the DOF fields. */
    struct DofKernel
    {
        std::string name;
        bool has_trial = false;
        OperandRef test;
        OperandRef trial;
        ParameterRef parameter;
        std::vector<unsigned int> domains; ///< Empty means all.
    };

    /** @brief A term of a nodal or element field: field += parameter * operator(DOF field). */
    struct FieldKernel
    {
        std::string name;
        const Field *field = nullptr;
        OperandRef trial;
        ParameterRef parameter;
        std::vector<unsigned int> domains; ///< Empty means all.
    };

    class ElementType
    {
    public:
        ElementType(const ElementType &) = delete;
        ElementType &operator=(const ElementType &) = delete;

        const std::string &name() const { return name_; }
        unsigned int dimension() const { return dimension_; }
        const ElementTypeOptions &options() const { return options_; }
        const std::vector<std::string> &domains() const { return domains_; }

        const std::vector<Field> &fields(FieldKind kind) const;
        const Field *find_field(FieldKind kind, const std::string &name) const;

        /** @brief The operator @p op bound to @p field; the parser creates Identity for every field with an FE. */
        const Nabla &nabla(const Field &field, NablaOperator op) const;

        const std::vector<ParameterDefinition> &parameters() const { return parameters_; }
        const std::vector<NonlinearLaw> &laws() const { return laws_; }

        const std::vector<DofKernel> &dof_kernels() const { return dof_kernels_; }
        const std::vector<FieldKernel> &nodal_kernels() const { return nodal_kernels_; }
        const std::vector<FieldKernel> &element_kernels() const { return element_kernels_; }

        /** @brief True when every field of the kernel exists in @p domain and the kernel applies there. */
        bool active(const DofKernel &kernel, unsigned int domain) const;
        bool active(const FieldKernel &kernel, unsigned int domain) const;

        /** @brief reinit()s once each FE used in @p domain, and the geometry FE, on @p elem. */
        void reinit(const libMesh::Elem &elem, unsigned int domain) const;

        /** @brief JxW and quadrature points of the last reinit(). */
        const std::vector<libMesh::Real> &jxw() const { return geometry_fe_->get_JxW(); }
        const std::vector<libMesh::Point> &xyz() const { return geometry_fe_->get_xyz(); }

        /** @brief Element matrix index of every shape function of every field of one kind. */
        struct FieldIndices
        {
            std::vector<std::vector<unsigned int>> field; ///< By field index; empty for a field not active in the domain.
            unsigned int n_dofs = 0;                      ///< Size of the element matrices of the system.
        };

        /**
         * @brief Where, in the element vectors of @p system, each shape function of each field of @p kind sits.
         *
         * Follows the DofMap, so it holds for any variable and variable-group layout of the system. The
         * fields are those of @p kind active in @p domain; their columns are in shape-function order.
         */
        FieldIndices element_indices(const libMesh::System &system,
                                     const libMesh::Elem &elem,
                                     FieldKind kind,
                                     unsigned int domain) const;

    private:
        friend class ElementTypeParser;
        ElementType() = default;

        std::string name_;
        unsigned int dimension_ = 0;
        ElementTypeOptions options_;
        std::vector<std::string> domains_;

        std::vector<ParameterDefinition> parameters_;
        std::vector<NonlinearLaw> laws_;

        std::vector<Field> dof_fields_;
        std::vector<Field> nodal_fields_;
        std::vector<Field> element_fields_;
        std::vector<DofKernel> dof_kernels_;
        std::vector<FieldKernel> nodal_kernels_;
        std::vector<FieldKernel> element_kernels_;

        std::deque<Nabla> nablas_; ///< Stable addresses: kernels point into it.
        std::unique_ptr<libMesh::QBase> quadrature_;
        std::map<std::pair<int, int>, std::unique_ptr<libMesh::FEAbstract>> fes_; ///< By (family, order).
        std::unique_ptr<libMesh::FEBase> geometry_fe_;
        std::vector<std::vector<libMesh::FEAbstract *>> fes_by_domain_;
    };
} // namespace pytec3d
