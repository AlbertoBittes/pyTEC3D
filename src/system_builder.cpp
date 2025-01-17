#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
#include "../include/compute_derived_quantities.h"
#include "../include/assembler.h"
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
    transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
    ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
    formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
    order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
    dim = equation_systems.get_mesh().mesh_dimension();

    // Check if dim is 2 or 3
    if (dim != 2 && dim != 3)
    {
        libmesh_error_msg("Error: dim must be 2 or 3.");
    }
    // Check if fe_family is not "LAGRANGE" and dim is not 3
    if (fe_family != "LAGRANGE" && fe_family != "NEDELEC_ONE")
    {
        libmesh_error_msg("Error: Element type must be 'LAGRANGE' (nodal formulation) or 'NEDELEC_ONE' (edge formulation).");
    }
    // Check if fe_family is not "LAGRANGE" and dim is not 3
    if (fe_family != "LAGRANGE" && dim != 3)
    {
        libmesh_error_msg("Error: Mesh must be 3D for edge elements.");
    }
    // Check if order is "FIRST" or "SECOND"
    if (order != "FIRST" && order != "SECOND")
    {
        libmesh_error_msg("Error: order must be 'FIRST' or 'SECOND'.");
    }
    // Check if formulation is "AV" or "Hphi"
    if (formulation != "AV" && formulation != "Hphi")
    {
        libmesh_error_msg("Error: formulation must be 'AV' or 'Hphi'.");
    }

    libMesh::out << "Creating system..." << std::endl;
    if (ferromagnetic)
    {
        system_name = "NonlinearTransient";
        // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
        equation_systems.add_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
    }
    else
    {
        system_name = "PlasmaArcLinearSystem";
        // Add TransientLinearImplicitSystem as System and return its pointer
        equation_systems.add_system<TransientLinearImplicitSystem>(system_name);
    }
    equation_systems.add_system<ExplicitSystem>("DerivedQuantitiesSystem");
    equation_systems.add_system<ExplicitSystem>("DerivedQuantitiesSystemElem");
};

void SystemBuilder::setVariables()
{
    System &system = equation_systems.get_system(system_name);
    ExplicitSystem &derived_quantities_system = equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystem");
    ExplicitSystem &derived_quantities_system_elem = equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem");
    if (formulation == "AV")
    {
        if (dim == 2)
        {
            system.add_variable("A", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
            equation_systems.parameters.set<std::string>("FE_family_A") = fe_family;
            derived_quantities_system.add_variable("J", Utility ::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>("L2_" + fe_family));
            derived_quantities_system_elem.add_variable("J_elem", CONSTANT, MONOMIAL);
        }
        if (dim == 3 && fe_family == "LAGRANGE")
        {
            system.add_variable("A", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(fe_family + "_VEC"));
            equation_systems.parameters.set<std::string>("FE_family_A") = fe_family + "_VEC";
            derived_quantities_system.add_variable("J", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>("L2_" + fe_family + "_VEC"));
            derived_quantities_system_elem.add_variable("J_elem", CONSTANT, MONOMIAL_VEC);
        }
        if (fe_family == "LAGRANGE")
        {
            derived_quantities_system.add_variable("B", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>("L2_" + fe_family + "_VEC"));
            derived_quantities_system_elem.add_variable("B_elem", CONSTANT, MONOMIAL_VEC);
            derived_quantities_system.add_variable("F", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>("L2_" + fe_family + "_VEC"));
            derived_quantities_system_elem.add_variable("F_elem", CONSTANT, MONOMIAL_VEC);
        }
        // if (transient)
        // {

        if (dim == 3)
        {
            system.add_variable("V",
                                Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family),
                                equation_systems.parameters.get<std::set<subdomain_id_type> *>("conductive_subdomain"));
        }
        assembler = std::make_shared<AssemblerLinear>(equation_systems, config);
        system.attach_assemble_object(*assembler);
        derived_quantities_calculator = std::make_shared<DerivedQuantitiesCalculator>(assembler, equation_systems);
    }
    else
    {
        libmesh_error_msg("not implemented");
    }
    // Add more conditions as needed for other combinations
}

std::shared_ptr<DerivedQuantitiesCalculator> SystemBuilder::getDerivedQuantitiesCalculator()
{
    return derived_quantities_calculator;
}