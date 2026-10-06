/**
 * @file element_type.cpp
 * @brief ElementType: field and kernel queries, FE reinit and element matrix indices.
 */
#include "element_type.h"

#include <algorithm>
#include <unordered_map>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementType: " + message);
        }

        bool in_domains(const std::vector<unsigned int> &domains, unsigned int domain)
        {
            return domains.empty() || std::find(domains.begin(), domains.end(), domain) != domains.end();
        }
    } // namespace

    const std::vector<Field> &ElementType::fields(FieldKind kind) const
    {
        switch (kind)
        {
        case FieldKind::Dof:
            return dof_fields_;
        case FieldKind::Nodal:
            return nodal_fields_;
        case FieldKind::Element:
            return element_fields_;
        }
        fail("unknown field kind.");
    }

    const Field *ElementType::find_field(FieldKind kind, const std::string &name) const
    {
        for (const Field &field : fields(kind))
            if (field.name == name)
                return &field;
        return nullptr;
    }

    bool ElementType::active(const DofKernel &kernel, unsigned int domain) const
    {
        return in_domains(kernel.domains, domain) && kernel.test.field().active_in(domain) &&
               (!kernel.has_trial || kernel.trial.field().active_in(domain));
    }

    bool ElementType::active(const FieldKernel &kernel, unsigned int domain) const
    {
        return in_domains(kernel.domains, domain) && kernel.field->active_in(domain) && kernel.trial.field().active_in(domain);
    }

    const Nabla &ElementType::nabla(const Field &field, NablaOperator op) const
    {
        for (const Nabla &candidate : nablas_)
            if (&candidate.field() == &field && candidate.op() == op)
                return candidate;
        fail("no operator " + std::to_string(static_cast<int>(op)) + " for field '" + field.name + "'.");
    }

    void ElementType::reinit(const libMesh::Elem &elem, unsigned int domain) const
    {
        for (libMesh::FEAbstract *fe : fes_by_domain_.at(domain))
            fe->reinit(&elem);
        geometry_fe_->reinit(&elem);
    }

    ElementType::FieldIndices ElementType::element_indices(const libMesh::System &system,
                                                           const libMesh::Elem &elem,
                                                           FieldKind kind,
                                                           unsigned int domain) const
    {
        const libMesh::DofMap &dof_map = system.get_dof_map();

        std::vector<libMesh::dof_id_type> element_dofs;
        dof_map.dof_indices(&elem, element_dofs);
        std::unordered_map<libMesh::dof_id_type, unsigned int> position;
        for (unsigned int i = 0; i < element_dofs.size(); ++i)
            position[element_dofs[i]] = i;

        FieldIndices out;
        out.n_dofs = static_cast<unsigned int>(element_dofs.size());
        out.field.resize(fields(kind).size());

        std::size_t covered = 0;
        for (const Field &field : fields(kind))
        {
            if (!field.active_in(domain))
                continue;

            std::vector<std::vector<libMesh::dof_id_type>> per_variable;
            for (const LibMeshVariable &variable : field.libmesh_variables)
            {
                std::vector<libMesh::dof_id_type> variable_dofs;
                dof_map.dof_indices(&elem, variable_dofs, system.variable_number(variable.label));
                per_variable.push_back(std::move(variable_dofs));
            }

            // One libMesh variable per component: the shape functions run node by node over the components.
            std::vector<libMesh::dof_id_type> dofs;
            if (per_variable.size() == 1)
                dofs = per_variable.front();
            else
            {
                const std::size_t n_components = per_variable.size();
                const std::size_t n_nodes = per_variable.front().size();
                dofs.resize(n_components * n_nodes);
                for (std::size_t component = 0; component < n_components; ++component)
                {
                    if (per_variable[component].size() != n_nodes)
                        fail("the components of '" + field.name + "' have different DOF counts.");
                    for (std::size_t node = 0; node < n_nodes; ++node)
                        dofs[node * n_components + component] = per_variable[component][node];
                }
            }

            std::vector<unsigned int> &columns = out.field[field.index];
            for (const libMesh::dof_id_type dof : dofs)
            {
                const auto it = position.find(dof);
                if (it == position.end())
                    fail("a DOF of '" + field.name + "' is not a DOF of the element.");
                columns.push_back(it->second);
            }
            covered += columns.size();
        }

        if (covered != element_dofs.size())
            fail("the fields of kind active in the domain cover " + std::to_string(covered) + " of " +
                 std::to_string(element_dofs.size()) + " element DOFs of system '" + system.name() + "'.");
        return out;
    }
} // namespace pytec3d
