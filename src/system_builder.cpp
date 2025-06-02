#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "system_builder.h"
#include "element_types.h"
#include "libmesh/string_to_enum.h"
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
// The definition of a geometric element

using namespace libMesh;

void SystemBuilder::initialize()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;

    auto equation_systems = this->get_parent_object();
    auto element_types_collector = equation_systems->parameters.get<std::shared_ptr<ElementTypesCollector>>("element_types_collector");
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    auto &system = equation_systems->get_system(this->system_name);
    if (system.name() == "ApecMLinearSystem")
    {
        vars_set = "DOFVariables";

        std::shared_ptr<AssemblerLinear> assembler_linear = std::make_shared<AssemblerLinear>(equation_systems);
        system.attach_assemble_object(*assembler_linear);
        equation_systems->parameters.set<std::shared_ptr<AssemblerLinear>>("assembler") = assembler_linear;

        // libMesh::out << "Setting boundary conditions object..." << std::endl;
        // NodeBCObject node_dof_constraints(system);
        // system.attach_constraint_object(node_dof_constraints);
    }
    else if (system.name() == "DerivedQuantitiesSystem")
    {
        vars_set = "NodalSolutionVariables";
    }
    else if (system.name() == "DerivedQuantitiesSystemElem")
    {
        vars_set = "ElementSolutionVariables";
    }
    auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds(vars_set);
    for (const auto &pair : vars_to_subdomain_ids)
    {
        const std::string &variable = pair.first;
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->getElementTypes().at(*subdomain_ids.begin());
        // Add the variable to the system
        const auto &element_data = elem_type->get_element_type_data();
        const auto &dof_variables = element_data.at(vars_set);
        auto dof_variable_data = dof_variables.at(variable);

        // Retrieve the 6th value from the tuple
        fe_fam = std::get<5>(dof_variable_data);
        order = std::get<6>(dof_variable_data);

        system.add_variable(variable, Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_fam), &subdomain_ids);
    }
    if (system.name() == "ApecMLinearSystem")
    {
        system.add_vector("base_rhs");
    }
    // }
    // Add more conditions as needed for other combinations
}