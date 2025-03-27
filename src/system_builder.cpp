#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
#include "../include/assembler.h"
#include "../include/session.h"
#include "../include/element_types.h"
#include "libmesh/string_to_enum.h"
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
// The definition of a geometric element
#include "vtkDataArray.h"

using namespace libMesh;

void SystemBuilder::initialize()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;
    auto element_types = equation_systems.parameters.get<std::shared_ptr<ElementTypesMapCSV>>("element_types");
    auto element_types_collector = std::make_shared<ElementTypesCollector>(element_types);
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    auto &system = equation_systems.get_system(k);
    if (system.name() == "PlasmaArcLinearSystem")
    {
        vars_set = "DOFVariables";
    }
    else if (system.name() == "DerivedQuantitiesSystem")
    {
        vars_set = "NodalDerivedQuantitiesVariables";
    }
    else if (system.name() == "DerivedQuantitiesSystemElem")
    {
        vars_set = "ElementDerivedQuantitiesVariables";
    }
    auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds(vars_set);
    for (const auto &pair : vars_to_subdomain_ids)
    {
        const std::string &variable = pair.first;
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->getElementTypes().at(*subdomain_ids.begin());
        // Add the variable to the system
        fe_fam = elem_type.getVariableFEFamily(variable);
        order = elem_type.getVariableOrder(variable);

        system.add_variable(variable, Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_fam), &subdomain_ids);
    }
    // }
    // Add more conditions as needed for other combinations
}