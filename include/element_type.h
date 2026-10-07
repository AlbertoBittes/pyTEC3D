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
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/system.h>

#include "element_field.h"
#include "nabla.h"
#include "parameter_value.h"
#include "registries.h"

namespace pytec3d
{
    using registry::ParameterSource;

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

    /** @brief A nonlinear law, e.g. BH: acts on the trial value of a kernel instead of multiplying it. */
    struct NonlinearLaw
    {
        registry::LawId id = registry::LawId::COUNT;                 ///< Used in kernels by its name.
        registry::ParameterId linear = registry::ParameterId::COUNT; ///< Tensor parameter the law reduces to when the element has no curve.
    };

    enum class ParameterKind : unsigned char
    {
        One,       ///< No parameter: 1.
        Parameter, ///< A registry::ParameterId.
        Law,       ///< A registry::LawId.
    };

    /** @brief How a kernel uses a parameter: "-nu", "BH". */
    struct ParameterRef
    {
        ParameterKind kind = ParameterKind::One;
        unsigned int index = 0; ///< registry::ParameterId or registry::LawId, by kind.
        ParameterType type = ParameterType::Scalar;
        libMesh::Real sign = 1.0;

        registry::ParameterId parameter() const { return static_cast<registry::ParameterId>(index); }
        registry::LawId law() const { return static_cast<registry::LawId>(index); }
        std::string_view name() const
        {
            return kind == ParameterKind::Parameter ? registry::name(parameter())
                   : kind == ParameterKind::Law     ? registry::name(law())
                                                    : std::string_view("1");
        }
    };

    /** @brief K/C block (with a trial) or right-hand-side source (without) on the DOF fields. */
    struct DofKernel
    {
        std::string name;
        bool has_trial = false;
        OperandRef test;
        OperandRef trial;
        ParameterRef parameter;
        std::vector<registry::Domain> domains; ///< Empty means all.
    };

    /** @brief A term of a nodal or element field: field += parameter * operator(DOF field). */
    struct FieldKernel
    {
        std::string name;
        const Field *field = nullptr;
        OperandRef trial;
        ParameterRef parameter;
        std::vector<registry::Domain> domains; ///< Empty means all.
    };

    class ElementType
    {
    public:
        ElementType(const ElementType &) = delete;
        ElementType &operator=(const ElementType &) = delete;

        registry::ElementTypeId id() const { return id_; }
        std::string_view name() const { return registry::name(id_); }
        unsigned int dimension() const { return dimension_; }
        const ElementTypeOptions &options() const { return options_; }

        const std::vector<Field> &fields(FieldKind kind) const;
        const Field *find_field(FieldKind kind, registry::FieldId id) const;

        /** @brief The operator @p op bound to @p field; the parser creates Identity for every field with an FE. */
        const Nabla &nabla(const Field &field, NablaOperator op) const;

        /** @brief The law @p id of this element type. */
        const NonlinearLaw &law(registry::LawId id) const;

        const std::vector<DofKernel> &dof_kernels() const { return dof_kernels_; }
        const std::vector<FieldKernel> &nodal_kernels() const { return nodal_kernels_; }
        const std::vector<FieldKernel> &element_kernels() const { return element_kernels_; }

        /** @brief True when every field of the kernel exists in @p domain and the kernel applies there. */
        bool active(const DofKernel &kernel, registry::Domain domain) const;
        bool active(const FieldKernel &kernel, registry::Domain domain) const;

        /** @brief reinit()s once each FE used in @p domain, and the geometry FE, on @p elem. */
        void reinit(const libMesh::Elem &elem, registry::Domain domain) const;

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
                                     registry::Domain domain) const;

    private:
        friend class ElementTypeParser;
        ElementType() = default;

        registry::ElementTypeId id_ = registry::ElementTypeId::COUNT;
        unsigned int dimension_ = 0;
        ElementTypeOptions options_;

        std::vector<std::optional<NonlinearLaw>> laws_; ///< By registry::LawId; empty for a law the type does not define.

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
