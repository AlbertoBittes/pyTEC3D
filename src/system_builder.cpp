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
#include "fe_parent_object.h"
#include "assembler.h"
#include "system_type.h"
#include "fem_data_repository.h"
#include "apec_m_system_objects.h"
#include "boundary_conditions_object.h"
// The definition of a geometric element

using namespace libMesh;

// This could actually be done in the
void SystemBuilder::initialize()
{
    System *system = this->get_parent_object();

    if (system->name() == "ApecMLinearSystem")
    {
        auto main_system = dynamic_cast<TransientLinearImplicitSystem *>(system);

        main_system->zero_out_matrix_and_rhs = false;

        main_system->current_local_solution->zero();
        main_system->current_local_solution->close();

        main_system->add_vector("solution_derivative");
        main_system->add_vector("solution_derivative_old");
    }
};

void SystemBuilder::init_variables()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;

    System *system = this->get_parent_object();
    auto element_types_collector = context<ElementTypeFactory>();
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    // Map system name to SystemType
    SystemType type = get_system_type_from_name(system->name());
    Globals::VariableCategory category = get_variable_category(type);
    auto vars_to_subdomain_ids = element_types_collector->associate_variables_with_element_IDs(category);

    for (const auto &pair : *vars_to_subdomain_ids)
    {
        const std::string &variable_string =
            Globals::enum_to_string(pair.first, Globals::ENUM_TO_STRING);
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());

        // Get variable info for this variable
        auto &var_info = elem_type->get_variable_information();
        auto it = var_info.find(pair.first);
        if (it == var_info.end())
            continue; // Skip if variable info not found

        VariableData &data = it->second;
        if (!data.is_active)
            continue; // Skip inactive variables

        // Add the variable to the system, specifying subdomains
        system->add_variable(
            variable_string,
            data.order,
            data.fe_family,
            &subdomain_ids);
    }

    std::vector<unsigned int> all_variable_numbers;
    system->get_all_variable_numbers(all_variable_numbers);
    for (unsigned int var_num : all_variable_numbers)
    {
        // Get the variable name from the system
        const std::string &var_name = system->variable_name(var_num);

        // Loop over all elements in the collector
        for (const auto &elem_pair : *(element_types_collector->get_data_map()))
        {
            const auto &elem_type = elem_pair.second;
            auto &var_info = elem_type->get_variable_information();

            // Find the variable enum corresponding to var_name
            for (auto &info_pair : var_info)
            {
                const std::string &enum_name = Globals::enum_to_string(info_pair.first, Globals::ENUM_TO_STRING);
                if (enum_name == var_name)
                {
                    info_pair.second.system_var_num = var_num;
                    break;
                }
            }
        }
    }
    // Add more conditions as needed for other combinations
}