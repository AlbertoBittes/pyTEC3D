// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Introduction Example 1 - Creation of a Mesh Object</h1>
// \author Benjamin S. Kirk
// \date 2003

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "simulation_controller.h"

using namespace libMesh;

// void read_parse_element_types(std::string element_types_filepath, EquationSystems &equation_systems)
// {
//     libMesh::out << "Reading element types..." << std::endl;
//     equation_systems.parameters.set<std::shared_ptr<ElementTypesCollector>>("element_types_collector") = std::make_shared<ElementTypesCollector>(element_types_filepath);
// }

// Function to execute the batch mode logic

// Main function
int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries
    LibMeshInit init(argc, argv);
    {
        if (init.comm().rank() != 0)
        {
            int i = 0;
            while (0 == i)
                sleep(1);
        }
    }
    // Check for proper calling arguments
    if (argc < 2)
    {
        libmesh_error_msg("Usage:\n"
                          << "\t" << argv[0] << " -b mesh_filename\n"
                          << "\t" << argv[0] << " -i");
    }

    // Parse command-line arguments
    std::string mode = argv[1];
    if (mode == "-b")
    {
        if (argc < 3)
        {
            libmesh_error_msg("Error: Missing mesh filename. Usage:\n"
                              << "\t" << argv[0] << " -b mesh_filename");
        }

        // Get the mesh filename
        std::string mesh_filename = argv[2];

        // Construct the SimulationController in batch mode
        SimulationController controller(init.comm(), true, mesh_filename);
        controller.run_batch_mode();
    }
    else if (mode == "-i")
    {
        libMesh::out << "Running in interactive mode..." << std::endl;

        // Construct the SimulationController in interactive mode
        SimulationController controller(init.comm(), false);
        controller.run_cli(init.comm());
    }
    else
    {
        libmesh_error_msg("Error: Invalid mode. Usage:\n"
                          << "\t" << argv[0] << " -b mesh_filename\n"
                          << "\t" << argv[0] << " -i");
    }

    // All done
    return 0;
}

// This will be done much later

// for (auto &elem : mesh.active_local_element_ptr_range())
// {
//     for ( auto &side : elem->boundary_sides_begin)
// }
// BoundaryMesh boundary_mesh(mesh.comm());

// boundaryInfo.print_summary();
// // std::map<dof_id_type, dof_id_type> node_id_map;
// // std::map<dof_id_type, unsigned char> side_id_map;
// // boundary_mesh.print_info();
// // boundaryInfo.get_side_and_node_maps(boundary_mesh, node_id_map, side_id_map);
// // iterate over the boundary mesh
// for (auto &elem : boundary_mesh.active_local_element_ptr_range())
// {

// }
// Put all this part in a function, it must also add an array of bools called conductive to make it easier to identify conductive elements
// init.comm().set_union(elem_types);
// std::map<subdomain_id_type, subdomain_id_type> subdomain_map;
// subdomain_id_type subdomain_id = 0;
// for (const auto &elem_type : elem_types)
// {
//     subdomain_map[elem_type] = subdomain_id;
//     subdomain_id++;
// }
// BoundaryInfo &boundaryInfo = mesh.get_boundary_info(); // Retrieves the BoundaryInfo object
// // boundaryInfo.sync(boundary_mesh);
// {
//     Point center(0., 0., 0.);
//     Point axis(0., 0., 1.);
//     Real angle = -90 * 2 * libMesh::pi / 360;
//     AzimuthalPeriodicBoundary periodic_bc(center, axis, angle);
//     periodic_bc.myboundary = 7;
//     periodic_bc.pairedboundary = 8;
//     equation_systems.get_system<TransientLinearImplicitSystem>(0).get_dof_map().add_periodic_boundary(periodic_bc);
// }
// vtkDataArray *orig_id = partitionedGrid->GetPointData()->GetArray("VTKorigID");
// vtkDataArray *boundary_id = partitionedGrid->GetPointData()->GetArray("boundary_id");
// Print the names of the arrays in the vtkPointData
// for (int i = 0; i < partitionedGrid->GetPointData()->GetNumberOfArrays(); i++)
// {
//     libMesh::out << "Array " << i << ": " << partitionedGrid->GetPointData()->GetArrayName(i) << std::endl;
// }
// if (!boundary_id)
// {
//     libmesh_error_msg("Error: 'boundary_id' array not found in vtkPointData.");
// }

// vtkIdType size_vtk = orig_id->GetNumberOfTuples();
// dof_id_type orig_id_num;
// const Node *node_ptr = nullptr;
// unsigned int val;
// for (vtkIdType j = 0; j < size_vtk; j++)
// {
//     orig_id_num = static_cast<dof_id_type>(orig_id->GetTuple1(j));
//     node_ptr = mesh.node_ptr(orig_id_num);
//     val = boundary_id->GetTuple1(j);
//     if (val != 0)
//         boundaryInfo.add_node(node_ptr, boundary_id->GetTuple1(j));
// }
// boundaryInfo.build_side_list_from_node_list();
// boundaryInfo.clear_boundary_node_ids();
// boundaryInfo.print_info();
// boundaryInfo.print_summary();