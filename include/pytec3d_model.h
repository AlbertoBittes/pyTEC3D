/**
 * @file pytec3d_model.h
 * @brief C++ model container for mesh/equation-system ownership in pytec3d.
 */
#ifndef pyTEC3D_MODEL_H
#define pyTEC3D_MODEL_H

#include <memory>
#include <string>

#include <libmesh/distributed_mesh.h>
#include <pybind11/pybind11.h>

#include "augment_sparsity_on_nodes.h"

class vtkUnstructuredGrid;
class pyTEC3DEquationSystems;

namespace libMesh
{
    class LibMeshInit;
}

class PYBIND11_EXPORT pyTEC3DModel
{
public:
    pyTEC3DModel() = default;
    ~pyTEC3DModel();

    void load_into_solver(const std::shared_ptr<libMesh::LibMeshInit> &libmesh_init,
                          pybind11::object pytec3dodel_object);
    void solve();
    void write_system_to_vtk_grid(const std::string &system_name);
    void clear();

    bool is_loaded() const;
    pybind11::object get_pytec3d_object() const;
    vtkUnstructuredGrid *get_active_vtk_grid() const;
    const AugmentSparsityOnNodes *get_coupling_functor() const;

private:
    std::unique_ptr<libMesh::DistributedMesh> mesh;
    std::shared_ptr<pyTEC3DEquationSystems> equation_system;
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor;
    vtkUnstructuredGrid *active_vtk_grid = nullptr;
    pybind11::object pytec3d_object;
    bool loaded = false;
    bool closed_ = false;
};

#endif // pyTEC3D_MODEL_H
