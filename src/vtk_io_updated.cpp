// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes
#include <fstream>
#include "../include/vtk_io_updated.h"
// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"

#ifdef LIBMESH_HAVE_VTK

// I get a lot of "warning: extra ';' inside a class [-Wextra-semi]" from clang
// on VTK header files.
#include "libmesh/ignore_warnings.h"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConfigure.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkMPI.h"

#include "libmesh/restore_warnings.h"

// A convenient macro for comparing VTK versions.  Returns 1 if the
// current VTK version is < major.minor.subminor and zero otherwise.
bool ends_with(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.rfind(suffix) == (str.size() - suffix.size());
}
// It relies on the VTK version numbers detected during configure.  Note that if
// LIBMESH_HAVE_VTK is not defined, none of the LIBMESH_DETECTED_VTK_VERSION_* variables will
// be defined either.
#define VTK_VERSION_LESS_THAN(major, minor, subminor)                                                            \
    ((LIBMESH_DETECTED_VTK_VERSION_MAJOR < (major) ||                                                            \
      (LIBMESH_DETECTED_VTK_VERSION_MAJOR == (major) && (LIBMESH_DETECTED_VTK_VERSION_MINOR < (minor) ||         \
                                                         (LIBMESH_DETECTED_VTK_VERSION_MINOR == (minor) &&       \
                                                          LIBMESH_DETECTED_VTK_VERSION_SUBMINOR < (subminor))))) \
         ? 1                                                                                                     \
         : 0)

#endif // LIBMESH_HAVE_VTK

// using namespace libMesh;

// Constructor for reading
VTKIOUpdated::VTKIOUpdated(MeshBase &mesh) : MeshInput<MeshBase>(mesh, /*is_parallel_format=*/true),
                                             MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true)
#ifdef LIBMESH_HAVE_VTK
                                             ,
                                             _compress(false)
#endif
{
}

// Constructor for writing
VTKIOUpdated::VTKIOUpdated(const MeshBase &mesh) : MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true)
#ifdef LIBMESH_HAVE_VTK
                                                   ,
                                                   _compress(false)
#endif
{
}

// The rest of the file is wrapped in ifdef LIBMESH_HAVE_VTK except for
// a couple of "stub" functions at the bottom.
#ifdef LIBMESH_HAVE_VTK

void VTKIOUpdated::write_data(const MeshBase &mesh,
                              const std::string &fname,
                              const NumericVector<Number> &soln,
                              const std::vector<std::string> &names,
                              const NumericVector<Number> &soln_elem,
                              const std::vector<std::string> &names_elem)
{
    // Warn that the .pvtu file extension should be used.  Paraview
    // recognizes this, and it works in both serial and parallel.  Only
    // warn about this once.
    // if (fname.substr(fname.rfind("."), fname.size()) != ".pvtu")
    //     libmesh_do_once(libMesh::err << "The .pvtu extension should be used when writing VTK files in libMesh.");

    // If there are variable names being written, the solution vector
    // should not be empty, it should have been broadcast to all
    // processors by the MeshOutput base class, since VTK is a parallel
    // format.  Verify this before going further.
    // if (!names.empty() && soln.size() && soln_elem.size())
    // {
    //     libmesh_error_msg("Empty soln vector in VTKIOUpdated::write_nodal_data().");
    std::cout << "Size of soln: " << soln.size() << std::endl;
    std::cout << "Size of soln_elem: " << soln_elem.size() << std::endl;
// }
// we only use Unstructured grids
// _vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

// // add nodes to the grid and update _local_node_map
// _local_node_map.clear();
// this->nodes_to_vtk();
#ifdef LIBMESH_HAVE_MPI
    // Set VTK to the same communicator as libMesh
    vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New();
    MPI_Comm mpi_comm = mesh.comm().get();
    vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
    vtk_comm->InitializeExternal(&vtk_opaque_comm);

    vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New();
    vtk_mpi_ctrl->SetCommunicator(vtk_comm);

#endif
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    std::string filename = "output" + std::to_string(mesh.processor_id()) + ".vtu";
    writer->SetFileName(filename.c_str());
    writer->SetInputData(_vtk_grid);
    writer->Write();
    // // add cells to the grid
    // this->cells_to_vtk();
    // vtkSmartPointer<vtkDataArray> orig_nodes = vtkSmartPointer<vtkDataArray>::New();
    // vtkSmartPointer<vtkDataArray> orig_elems = vtkSmartPointer<vtkDataArray>::New();
    // orig_nodes = _vtk_grid->GetPointData()->GetArray("VTKorigID");
    // orig_elems = _vtk_grid->GetCellData()->GetArray("ansys_elem_num");
    // add nodal solutions to the grid, if solutions are given
    if (names.size() > 0)
    {
        std::size_t num_vars = names.size();
        dof_id_type num_nodes = mesh.n_local_nodes();

        for (std::size_t variable = 0; variable < num_vars; ++variable)
        {
            if (variable > 2)
                break;
            char n_comps = 1;
            if (ends_with(names[variable].c_str(), "_x"))
            {
                if (variable + 1 < num_vars)
                    if (ends_with(names[variable + 1].c_str(), "_y"))
                        n_comps++;
                if (variable + 2 < num_vars)
                    if (ends_with(names[variable + 2].c_str(), "_z"))
                        n_comps++;
            }
            vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();

            data->SetName(names[variable].c_str());
            std::size_t num_nodes_vtk = static_cast<std::size_t>(num_nodes);

            if (n_comps >= 2)
                data->SetNumberOfComponents(n_comps);
            // number of local and ghost nodes
            data->SetNumberOfValues(num_nodes_vtk * n_comps);
            // loop over all nodes and get the solution for the current
            // variable, if the node is in the current partition
            // Write _vtk_grid
            unsigned int j = 0;
            for (dof_id_type k = soln.first_local_index(); k < soln.last_local_index(); k++) //_vtk_grid->GetPointData()->GetArray("VTKorigID")->GetNumberOfTuples(); ++k)
            {
                // unsigned int orig_id = _vtk_grid->GetPointData()->GetArray("VTKorigID")->GetTuple1(k);
                // std::map<dof_id_type, dof_id_type>::iterator local_node_it = _local_node_map.find(k);
                // if (local_node_it == _local_node_map.end())
                //     continue; // not a local node

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                libmesh_do_once(libMesh::err << "Only writing the real part for complex numbers!\n"
                                             << "if you need this support contact " << LIBMESH_PACKAGE_BUGREPORT << std::endl);
                data->SetValue(local_node_it->second, soln[k * num_vars + variable].real());
#else
                if (n_comps == 1)
                    data->SetValue(k / num_vars, soln(k));
                else
                {
                    char i;
                    for (i = 0; i < n_comps; ++i)
                    {
                        data->SetComponent(k / num_vars, i, soln(k + i));
                    }
                    k = k + i;
                }

#endif
            }
            _vtk_grid->GetPointData()->AddArray(data);
            if (n_comps > 1)
                variable += n_comps - 1;
        }
    }
    //     if (names_elem.size() > 0)
    //     {
    //         std::size_t num_vars = names_elem.size();
    //         dof_id_type num_elems = mesh.n_local_elem();

    //         for (std::size_t variable = 0; variable < num_vars; ++variable)
    //         {
    //             char n_comps = 1;
    //             if (ends_with(names_elem[variable].c_str(), "_x"))
    //             {
    //                 if (variable + 1 < num_vars)
    //                     if (ends_with(names_elem[variable + 1].c_str(), "_y"))
    //                         n_comps++;
    //                 if (variable + 2 < num_vars)
    //                     if (ends_with(names_elem[variable + 2].c_str(), "_z"))
    //                         n_comps++;
    //             }

    //             vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    //             data->SetName(names_elem[variable].c_str());
    //             if (n_comps >= 2)
    //                 data->SetNumberOfComponents(n_comps);
    //             // number of local and ghost nodes
    //             data->SetNumberOfTuples(num_elems);

    //             // loop over all nodes and get the solution for the current
    //             // variable, if the node is in the current partition
    //             for (dof_id_type k = soln_elem.first_local_index(); k < soln_elem.last_local_index(); k++)
    //             {
    //                 // std::map<dof_id_type, dof_id_type>::iterator local_node_it = _local_node_map.find(k);
    //                 // if (local_node_it == _local_node_map.end())
    //                 //     continue; // not a local node

    // #ifdef LIBMESH_USE_COMPLEX_NUMBERS
    //                 libmesh_do_once(libMesh::err << "Only writing the real part for complex numbers!\n"
    //                                              << "if you need this support contact " << LIBMESH_PACKAGE_BUGREPORT
    //                                              << std::endl);
    //                 data->SetValue(local_node_it->second, soln[k * num_vars + variable].real());
    // #else
    //                 if (n_comps == 1)
    //                     data->SetValue(k / num_vars, soln_elem(k));
    //                 else
    //                 {
    //                     char i = 0;
    //                     for (i = 0; i < n_comps; ++i)
    //                     {
    //                         data->SetComponent(k / num_vars, i, soln_elem(k));
    //                     }
    //                     k = k + i;
    //                 }
    // #endif
    //                 if (_vtk_grid->GetCellData()->HasArray(names_elem[variable].c_str()))
    //                 {
    //                     _vtk_grid->GetCellData()->RemoveArray(names_elem[variable].c_str());
    //                 }
    //             }
    //             _vtk_grid->GetCellData()->AddArray(data);
    //             if (n_comps > 1)
    //                 variable += n_comps - 1;
    //         }
    //     }
    // Write the partitioned grids in PVTU format
    std::cout << "Writing partitioned grids in PVTU format" << std::endl;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> ugWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    ugWriter->SetInputData(_vtk_grid);
    ugWriter->SetFileName((fname + "_" + std::to_string(mesh.processor_id()) + ".vtu").c_str());
    ugWriter->Write();
    // vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New();
    // MPI_Comm mpi_comm = mesh.comm().get();
    // vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
    // vtk_comm->InitializeExternal(&vtk_opaque_comm);

    // vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New();
    // vtk_mpi_ctrl->SetCommunicator(vtk_comm);
    if (mesh.processor_id() == 0)
    {
        // vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        // writer->SetController(vtk_mpi_ctrl);
        writer->SetNumberOfPieces(mesh.n_processors());
        // writer->SetStartPiece(mesh.processor_id());
        // writer->SetEndPiece(mesh.processor_id());
        writer->SetInputData(_vtk_grid);
        writer->SetFileName((fname + ".pvtu").c_str());
        writer->SetGhostLevel(1);
        std::cout << "Ciao 0";
        writer->Write();
    }
    else
    {
        std::cout << "Ciao";
    }
}

vtkUnstructuredGrid *VTKIOUpdated::get_vtk_grid()
{
    return _vtk_grid;
}

void VTKIOUpdated::set_compression(bool b)
{
    this->_compress = b;
}

void VTKIOUpdated::write_equation_systems(const std::string &fname,
                                          const EquationSystems &es,
                                          const std::set<std::string> *system_names)
{
    LOG_SCOPE("write_equation_systems()", "MeshOutput");

    // If we're asked to write data that's associated with a different
    // mesh, output files full of garbage are the result.
    const libMesh::MeshBase &my_mesh = es.get_mesh();

    // A non-renumbered mesh may not have a contiguous numbering, and
    // that needs to be fixed before we can build a solution vector.
    // Define a constant set with "DerivedQuantitiesSystemElem" as the only element
    const std::set<std::string> system_set = {"PlasmaArcLinearSystem", "DerivedQuantitiesSystem"};
    std::set<std::string> derived_quantities_system_elem_set = {"DerivedQuantitiesSystemElem"};

    // Create a pointer to the constant set
    std::set<std::string> *const_derived_quantities_system_elem_set_ptr = &derived_quantities_system_elem_set;
    const std::set<std::string> *system_set_ptr = &system_set;
    // Build the list of variable names that will be written.
    std::vector<std::string> names;
    es.build_variable_names(names, nullptr, system_set_ptr);

    std::vector<std::string> names_elem;
    // subtract names_elem from names
    std::vector<std::string> system_names_elem_temp;
    es.build_variable_names(names_elem, nullptr, const_derived_quantities_system_elem_set_ptr);
    std::unique_ptr<NumericVector<Number>> parallel_soln =
        es.build_parallel_solution_vector(system_set_ptr);
    std::unique_ptr<NumericVector<Number>> parallel_soln_elem =
        es.build_parallel_elemental_solution_vector(system_names_elem_temp);
    // Remove all elements in names which are also in names_elem
    // names.erase(std::remove_if(names.begin(), names.end(),
    //                            [&names_elem](const std::string &name)
    //                            {
    //                                return std::find(names_elem.begin(), names_elem.end(), name) != names_elem.end();
    //                            }),
    //             names.end());
    write_data(my_mesh, fname, *parallel_soln, names, *parallel_soln_elem, names_elem);
}
#else // !LIBMESH_HAVE_VTK

void VTKIOUpdated::read(const std::string &name)
{
    libmesh_error_msg("Cannot read VTK file: " << name
                                               << "\nYou must have VTK installed and correctly configured to read VTK meshes.");
}

void VTKIOUpdated::write_nodal_data(const std::string &fname,
                                    const std::vector<Number> &,
                                    const std::vector<std::string> &)
{
    libmesh_error_msg("Cannot write VTK file: " << fname
                                                << "\nYou must have VTK installed and correctly configured to read VTK meshes.");
}

#endif // LIBMESH_HAVE_VTK
