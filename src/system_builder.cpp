#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "system_builder.h"
#include "element_types.h"
#include "assembler.h"
#include "fem_data_repository.h"
#include "apec_m_linear_system.h"
#include "libmesh/string_to_enum.h"
#include "boundary_conditions_object.h"
#include "libmesh/numeric_vector.h"
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

    System *system = this->get_parent_object();
    auto context = std::dynamic_pointer_cast<ElementPhysicsContext>(this->get_context());
    auto element_types_collector = context->get_element_types();
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    // Map system name to SystemType
    SystemType type = get_system_type_from_name(system->name());
    VariableCategory category = get_variable_category(type);
    auto vars_to_subdomain_ids = element_types_collector->associate_variables_with_element_IDs(category);

    for (const auto &pair : vars_to_subdomain_ids)
    {
        const std::string &variable_string =
            Globals::enum_to_string(pair.first, Globals::ENUM_TO_STRING);
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());

        // Get variable info for this variable
        const auto &var_info = elem_type->get_variable_information();
        auto it = var_info.find(pair.first);
        if (it == var_info.end())
            continue; // Skip if variable info not found

        const VariableData &data = it->second;
        if (!data.is_active)
            continue; // Skip inactive variables

        // Add the variable to the system, specifying subdomains
        system->add_variable(
            variable_string,
            data.order,
            data.fe_family,
            &subdomain_ids);
    }

    if (system->name() == "ApecMLinearSystem")
    {
        std::shared_ptr<AssemblerLinear> assembler_linear = std::make_shared<AssemblerLinear>(system);
        system->attach_assemble_object(*assembler_linear);

        auto main_system = dynamic_cast<ApecMLinearSystem *>(system);

        main_system->add_vector("base_rhs");

        main_system->zero_out_matrix_and_rhs = false;

        std::shared_ptr<BoundaryConditionsObject> node_bc_object = std::make_shared<BoundaryConditionsObject>();
        main_system->attach_constraint_object(*node_bc_object);
        main_system->set_boundary_conditions_object(node_bc_object);

        main_system->current_local_solution->zero();
        main_system->current_local_solution->close();
    }
    // }
    // Add more conditions as needed for other combinations
}