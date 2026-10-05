/**
 * @file element_type_instance.h
 * @brief A concrete element type: a JSON base type resolved against user options.
 *
 * Resolution rules:
 *  - DOF fields are LAGRANGE (scalar) or LAGRANGE_VEC (vector) of the user order, unless the field
 *    names another family in "finite_elements" (e.g. NEDELEC_ONE, always first order).
 *  - Nodal fields are L2_LAGRANGE or L2_LAGRANGE_VEC of the user order.
 *  - Element fields are MONOMIAL / MONOMIAL_VEC of CONSTANT order.
 *  - A vector field with a single active component becomes a scalar variable of the same name.
 *  - Fields restricted to other physical domains are dropped.
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe.h>
#include <libmesh/quadrature.h>

#include "element_type_registry.h"

namespace pytec3d
{
    /** @brief User-supplied options selecting and parameterising a base element type. */
    struct ElementOptions
    {
        std::string type;          ///< Base element type name.
        unsigned int domain = 0;   ///< Option 1: index into the registry domain list.
        unsigned int order = 1;    ///< Option 2: FE order (applies to Lagrange-type families only).
        bool axisymmetric = false; ///< Option 3.
        bool nonlinear = false;    ///< Option 4.
    };

    /** @brief One libMesh variable resolved from a field. */
    struct VariableSpec
    {
        std::string name;
        FieldKind kind = FieldKind::Dof;
        libMesh::FEFamily family = libMesh::LAGRANGE;
        libMesh::Order order = libMesh::FIRST;
        std::vector<unsigned int> active_components; ///< Physical components (0=x,1=y,2=z); one entry for scalars.
        bool physical_vector = false;                ///< Declared as a vector field, even with a single active component.
        bool vector_fe = false;                      ///< Evaluated with a vector FE (FEVectorBase).
        std::string nodal_loads;                     ///< VTK point array that loads this DOF variable, if any.

        bool is_vector() const { return active_components.size() > 1; }
    };

    class ElementTypeInstance
    {
    public:
        /** @brief Resolves @p options against @p registry; aborts on invalid input. */
        ElementTypeInstance(const ElementTypeRegistry &registry, const ElementOptions &options);

        const ElementOptions &options() const { return options_; }
        const BaseElementType &base() const { return *base_; }
        const std::string &domain_name() const { return domain_name_; }
        unsigned int dimension() const { return base_->dimension; }

        /** @brief All active variables in JSON order (dof, nodal, element). */
        const std::vector<VariableSpec> &variables() const { return variables_; }
        std::vector<VariableSpec> variables_of_kind(FieldKind kind) const;
        const VariableSpec *find_variable(const std::string &name) const;

        /** @brief Gauss order shared by all variables: 2 * (highest DOF order) + 1. */
        libMesh::Order quadrature_order() const;

        /** @brief Quadrature rule shared by every FE object of this instance (built on first use). */
        libMesh::QBase &quadrature() const;

        /** @brief Scalar FE of @p variable (built on first use); aborts if the variable uses a vector FE. */
        libMesh::FEBase &scalar_fe(const std::string &variable) const;

        /** @brief Vector FE of @p variable (built on first use); aborts if the variable uses a scalar FE. */
        libMesh::FEVectorBase &vector_fe(const std::string &variable) const;

        /** @brief First-order Lagrange FE used only for geometry (JxW, xyz) on the shared quadrature. */
        libMesh::FEBase &geometry_fe() const;

    private:
        struct FEEntry
        {
            std::unique_ptr<libMesh::FEBase> scalar;
            std::unique_ptr<libMesh::FEVectorBase> vector;
        };

        std::size_t variable_index(const std::string &name) const;

        ElementOptions options_;
        const BaseElementType *base_ = nullptr; ///< Owned by the registry, which must outlive this object.
        std::string domain_name_;
        std::vector<VariableSpec> variables_;

        // FE objects keep internal reinit state, so they are cached lazily behind const accessors.
        mutable std::unique_ptr<libMesh::QBase> quadrature_;
        mutable std::vector<FEEntry> fes_; ///< Parallel to variables_.
        mutable std::unique_ptr<libMesh::FEBase> geometry_fe_;
    };
} // namespace pytec3d
