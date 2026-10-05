/**
 * @file element_type_instance.cpp
 * @brief ElementTypeInstance: resolves JSON fields into libMesh variable descriptions.
 */
#include "element_type_instance.h"

#include <algorithm>
#include <numeric>

#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        libMesh::FEFamily family_for(FieldKind kind, bool vector)
        {
            switch (kind)
            {
            case FieldKind::Dof:
                return vector ? libMesh::LAGRANGE_VEC : libMesh::LAGRANGE;
            case FieldKind::Nodal:
                return vector ? libMesh::L2_LAGRANGE_VEC : libMesh::L2_LAGRANGE;
            case FieldKind::Element:
                return vector ? libMesh::MONOMIAL_VEC : libMesh::MONOMIAL;
            }
            logging::error_and_abort("ElementTypeInstance: unknown field kind.");
        }

        bool is_lagrange_type(libMesh::FEFamily family)
        {
            return family == libMesh::LAGRANGE || family == libMesh::LAGRANGE_VEC ||
                   family == libMesh::L2_LAGRANGE || family == libMesh::L2_LAGRANGE_VEC;
        }
    } // namespace

    ElementTypeInstance::ElementTypeInstance(const ElementTypeRegistry &registry, const ElementOptions &options)
        : options_(options), base_(&registry.get_type(options.type)), domain_name_(registry.domain_name(options.domain))
    {
        const std::string context = "ElementTypeInstance '" + options_.type + "': ";

        if (options_.order < static_cast<unsigned int>(libMesh::FIRST))
            logging::error_and_abort(context + "order must be >= 1.");
        if (options_.axisymmetric && base_->dimension != 2)
            logging::error_and_abort(context + "axisymmetric requires a 2D element type.");

        for (const FieldSpec &field : base_->fields)
        {
            if (!field.domains.empty() &&
                std::find(field.domains.begin(), field.domains.end(), domain_name_) == field.domains.end())
                continue;

            VariableSpec variable;
            variable.name = field.name;
            variable.kind = field.kind;

            // A vector field has every component of the element dimension; a single component makes it a scalar.
            if (!field.is_vector)
                variable.active_components = {0};
            else if (field.active_components.size() == 1)
                variable.active_components = field.active_components;
            else
            {
                variable.active_components.resize(base_->dimension);
                std::iota(variable.active_components.begin(), variable.active_components.end(), 0u);
            }

            variable.family = family_for(field.kind, variable.is_vector());
            if (!field.finite_elements.empty())
            {
                variable.family = libMesh::Utility::string_to_enum<libMesh::FEFamily>(field.finite_elements);
                // "LAGRANGE" on a multi-component field means component-wise Lagrange.
                if (variable.is_vector() && variable.family == libMesh::LAGRANGE)
                    variable.family = libMesh::LAGRANGE_VEC;
            }
            variable.physical_vector = field.is_vector;
            variable.nodal_loads = field.nodal_loads;
            variable.vector_fe = libMesh::FEInterface::field_type(variable.family) == libMesh::TYPE_VECTOR;
            if (variable.vector_fe != variable.is_vector())
                logging::error_and_abort(context + "field '" + field.name +
                                         "': FE family does not match the number of active components.");

            if (is_lagrange_type(variable.family))
                variable.order = static_cast<libMesh::Order>(options_.order);
            else if (field.kind == FieldKind::Element)
                variable.order = libMesh::CONSTANT;
            else
                variable.order = libMesh::FIRST;

            variables_.push_back(std::move(variable));
        }
        fes_.resize(variables_.size());
    }

    std::vector<VariableSpec> ElementTypeInstance::variables_of_kind(FieldKind kind) const
    {
        std::vector<VariableSpec> out;
        for (const auto &variable : variables_)
            if (variable.kind == kind)
                out.push_back(variable);
        return out;
    }

    const VariableSpec *ElementTypeInstance::find_variable(const std::string &name) const
    {
        for (const auto &variable : variables_)
            if (variable.name == name)
                return &variable;
        return nullptr;
    }

    std::size_t ElementTypeInstance::variable_index(const std::string &name) const
    {
        for (std::size_t i = 0; i < variables_.size(); ++i)
            if (variables_[i].name == name)
                return i;
        logging::error_and_abort("ElementTypeInstance '" + options_.type + "': unknown variable '" + name + "'.");
    }

    libMesh::Order ElementTypeInstance::quadrature_order() const
    {
        unsigned int highest = libMesh::FIRST;
        for (const auto &variable : variables_)
            if (variable.kind == FieldKind::Dof)
                highest = std::max(highest, static_cast<unsigned int>(variable.order));
        return static_cast<libMesh::Order>(2 * highest + 1);
    }

    libMesh::QBase &ElementTypeInstance::quadrature() const
    {
        if (!quadrature_)
            quadrature_ = std::make_unique<libMesh::QGauss>(dimension(), quadrature_order());
        return *quadrature_;
    }

    libMesh::FEBase &ElementTypeInstance::scalar_fe(const std::string &variable) const
    {
        const std::size_t index = variable_index(variable);
        const VariableSpec &spec = variables_[index];
        if (spec.vector_fe)
            logging::error_and_abort("ElementTypeInstance: variable '" + variable + "' uses a vector FE.");

        FEEntry &entry = fes_[index];
        if (!entry.scalar)
        {
            entry.scalar = libMesh::FEBase::build(dimension(), libMesh::FEType(spec.order, spec.family));
            entry.scalar->attach_quadrature_rule(&quadrature());
        }
        return *entry.scalar;
    }

    libMesh::FEVectorBase &ElementTypeInstance::vector_fe(const std::string &variable) const
    {
        const std::size_t index = variable_index(variable);
        const VariableSpec &spec = variables_[index];
        if (!spec.vector_fe)
            logging::error_and_abort("ElementTypeInstance: variable '" + variable + "' uses a scalar FE.");

        FEEntry &entry = fes_[index];
        if (!entry.vector)
        {
            entry.vector = libMesh::FEVectorBase::build(dimension(), libMesh::FEType(spec.order, spec.family));
            entry.vector->attach_quadrature_rule(&quadrature());
        }
        return *entry.vector;
    }

    libMesh::FEBase &ElementTypeInstance::geometry_fe() const
    {
        if (!geometry_fe_)
        {
            geometry_fe_ = libMesh::FEBase::build(dimension(), libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE));
            geometry_fe_->attach_quadrature_rule(&quadrature());
        }
        return *geometry_fe_;
    }
} // namespace pytec3d
