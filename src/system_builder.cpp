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

void SystemBuilder::createSystem()
{
    libMesh::out << "Checking system options..." << std::endl;
    // Here I will implement various checks about the system element types, the mesh and everything
    // ferromagnetic = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    // fe_family = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["fe_family"].get<std::string>();
    // formulation = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["formulation"].get<std::string>();
    // order = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["approx_order"].get<std::string>();
    // // dim = equation_systems.get_mesh().mesh_dimension();
    // dim = 3;
    // // Check if dim is 2 or 3
    // if (dim != 2 && dim != 3)
    // {
    //     libmesh_error_msg("Error: dim must be 2 or 3.");
    // }
    // // Check if fe_family is not "LAGRANGE" and dim is not 3
    // if (fe_family != "LAGRANGE" && fe_family != "NEDELEC_ONE")
    // {
    //     libmesh_error_msg("Error: Element type must be 'LAGRANGE' (nodal formulation) or 'NEDELEC_ONE' (edge formulation).");
    // }
    // // Check if fe_family is not "LAGRANGE" and dim is not 3
    // if (fe_family != "LAGRANGE" && dim != 3)
    // {
    //     libmesh_error_msg("Error: Mesh must be 3D for edge elements.");
    // }
    // // Check if order is "FIRST" or "SECOND"
    // if (order != "FIRST" && order != "SECOND")
    // {
    //     libmesh_error_msg("Error: order must be 'FIRST' or 'SECOND'.");
    // }
    // // Check if formulation is "AV" or "Hphi"
    // if (formulation != "AV" && formulation != "Hphi")
    // {
    //     libmesh_error_msg("Error: formulation must be 'AV' or 'Hphi'.");
    // }

    // libMesh::out << "Creating system..." << std::endl;
    // if (ferromagnetic)
    // {
    //     system_name = "NonlinearTransient";
    //     // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
    //     equation_systems.add_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
    // }
    // else
    // {
    system_name = "PlasmaArcLinearSystem";
    // Add TransientLinearImplicitSystem as System and return its pointer
    equation_systems.add_system<TransientLinearImplicitSystem>(system_name);
    // }
    equation_systems.add_system<ExplicitSystem>("DerivedQuantitiesSystem");
    equation_systems.add_system<ExplicitSystem>("DerivedQuantitiesSystemElem");
};

void SystemBuilder::setVariables()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;
    auto element_types = equation_systems.parameters.get<std::shared_ptr<ElementTypesMapCSV>>("element_types");
    auto element_types_collector = std::make_shared<ElementTypesCollector>(element_types);
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    {
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

            system.add_variable(variable, Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_fam), );
        }
    }
    // Add more conditions as needed for other combinations
}