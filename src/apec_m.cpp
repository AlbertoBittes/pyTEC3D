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

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "simulation_controller.h"

using namespace libMesh;

// Main function
// Here the code will understand if batch mode or interactive mode was chosen
// The simulation controller will take care of managing both modes, by taking as input
// The command line arguments

// To understand the code, you can follow the order of the functions called by the batch mode,
// to follow step by step how a simulation works
int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries
    LibMeshInit init(argc, argv);
    // {
    //     if (init.comm().rank() == 1)
    //     {
    //         int i = 0;
    //         while (0 == i)
    //             sleep(1);
    //     }
    // }

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
        // If you selected batch mode, you need to provide a model file (see documentation)
        if (argc < 3)
        {
            libmesh_error_msg("Error: Missing mesh filename. Usage:\n"
                              << "\t" << argv[0] << " -b mesh_filename");
        }

        // Get the mesh filename
        std::string mesh_filename = argv[2];

        // Construct the SimulationController in batch mode
        SimulationController controller(init.comm(), true, mesh_filename);
        // Execute the batch mode
        controller.run_batch_mode();
    }
    else if (mode == "-i")
    {
        libMesh::out << "Running in interactive mode..." << std::endl;

        // Construct the SimulationController in interactive mode
        SimulationController controller(init.comm(), false);
        // Execute the interactive mode
        controller.run_cli(init.comm());
    }
    else // raise error
    {
        libmesh_error_msg("Error: Invalid mode. Usage:\n"
                          << "\t" << argv[0] << " -b mesh_filename\n"
                          << "\t" << argv[0] << " -i");
    }

    // All done
    return 0;
}
