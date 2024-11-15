#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
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
};

void SystemBuilder::setVariables()
{
    System &system = equation_systems.get_system(system_name);

    if (formulation == "AV")
    {
        if (dim == 3 && fe_family == "LAGRANGE")
        {
            system.add_variable("AX", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
            system.add_variable("AY", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
        }
        system.add_variable("AZ", Utility::string_to_enum<Order>(order),
                            Utility::string_to_enum<FEFamily>(fe_family));
        if (transient)
        {
            std::set<subdomain_id_type> conductive_subdomain;

            conductive_subdomain.clear();
            conductive_subdomain.insert(1);
            system.add_variable("V",
                                Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family),
                                &conductive_subdomain);
        }
        AssemblerLinear *assembler = new AssemblerLinear(equation_systems, cell_data, config, materials);
        system.attach_assemble_object(*assembler);
    }
    else
    {
        libmesh_error_msg("not implemented");
    }
    // Add more conditions as needed for other combinations
}