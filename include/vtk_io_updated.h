#ifndef VTK_IO_UPDATED_H
#define VTK_IO_UPDATED_H

#include "libmesh/vtk_io.h"

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// Ignore "deprecated...header" warning from strstream
#include "libmesh/ignore_warnings.h"
#include "vtkType.h"
#include "vtkSmartPointer.h"
#include "libmesh/restore_warnings.h"
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
// C++ includes
#include <cstddef>
#include <map>

using namespace libMesh;
// Forward declarations
class vtkUnstructuredGrid;

class VTKIOUpdated : public MeshInput<MeshBase>,
                     public MeshOutput<MeshBase>
{
public:
    explicit VTKIOUpdated(MeshBase &mesh);

    explicit VTKIOUpdated(const MeshBase &mesh);

    // Implement any pure virtual functions from base classes
    virtual void read(const std::string &name) override {};
    virtual void write(const std::string &name) override {};
    // New constructor
    VTKIOUpdated(vtkSmartPointer<vtkUnstructuredGrid> vtkGrid)
        : _vtk_grid(vtkGrid), _compress(false)
    {
        // Uncomment this once you moved elsewhere the loads
        // // Remove all point data arrays
        // _vtk_grid->GetPointData()->Initialize();

        // // Remove all cell data arrays
        // _vtk_grid->GetCellData()->Initialize();
    }
    void write_data(const MeshBase &mesh,
                    const std::string &,
                    const NumericVector<Number> &,
                    const std::vector<std::string> &,
                    const NumericVector<Number> &soln_elem,
                    const std::vector<std::string> &names_elem);

    void write_equation_systems(const std::string &fname,
                                const EquationSystems &es,
                                const std::set<std::string> *system_names = nullptr);

    void set_compression(bool b);

    vtkUnstructuredGrid *get_vtk_grid();

private:
    // void system_vectors_to_vtk(const EquationSystems & es,
    //                            vtkUnstructuredGrid * & grid);

    vtkSmartPointer<vtkUnstructuredGrid> _vtk_grid;

    bool _compress;
};

#endif // VTK_IO_UPDATED_H