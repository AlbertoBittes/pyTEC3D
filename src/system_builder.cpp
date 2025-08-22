// C++ include files that we need
#include <iostream>
#include <set>

#include <libmesh/linear_implicit_system.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/transient_system.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>

#include "system_builder.h"
#include "fem_service_abstract.h"
#include "fe_parent_object.h"
#include "assembler.h"
#include "globals.h"
#include "element_types_factory.h"
#include "fem_data_repository.h"
#include "apec_m_system_objects.h"
#include "global_constraints_object.h"
// The definition of a geometric element

using namespace libMesh;

// This could actually be done in the
void SystemBuilder::initialize()
{
    System *system = this->get_parent_object();

    // Pass the same context view to the main system as used by SystemBuilder
    auto context_view = this->get_context_view();
    if (context_view)
    {
        auto apec_m_system = dynamic_cast<ApecMLinearSystem *>(system);
        if (apec_m_system)
        {
            apec_m_system->set_context_view(context_view);
        }
    }

    if (system->name() == "ApecMLinearSystem")
    {
        auto main_system = dynamic_cast<TransientLinearImplicitSystem *>(system);

        main_system->zero_out_matrix_and_rhs = false;

        main_system->current_local_solution->zero();
        main_system->current_local_solution->close();

        // main_system->add_vector("solution_derivative");
        // main_system->add_vector("solution_derivative_old");
    }
};

void SystemBuilder::init_variables()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;

    System *system = this->get_parent_object();
    auto element_types_collector = context<ElementTypesFactory>();
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    // Use new helper functions to map system name to SystemType and get category
    Globals::SystemType type = Globals::get_key_by_predicate_one(
        Globals::SYSTEM_INFO_MAP,
        [sys_name = system->name()](const Globals::SystemInfo &info)
        {
            return info.name == sys_name;
        });
    Globals::VariableCategory category = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type).category;
    auto vars_to_subdomain_ids = element_types_collector->associate_variables_with_element_IDs(category);

    for (const auto &pair : *vars_to_subdomain_ids)
    {
        const auto &variable_string = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, pair.first).global_label;
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());

        // Get variable info for this variable
        auto &var_info = elem_type->get_variable_information();
        auto it = var_info.find(pair.first);
        if (it == var_info.end())
            continue; // Skip if variable info not found

        Globals::VariableInfo &data = it->second;

        // Add the variable to the system, specifying subdomains
        system->add_variable(
            variable_string,
            data.order,
            data.fe_family,
            &subdomain_ids);
    }
    // Add more conditions as needed for other combinations
}