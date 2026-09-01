/**
 * @file apec_m_model.h
 * @brief C++ model container for mesh/equation-system ownership in APEC-M.
 */
#ifndef APEC_M_MODEL_H
#define APEC_M_MODEL_H

#include <memory>
#include <string>

#include <libmesh/distributed_mesh.h>
#include <pybind11/pybind11.h>

#include "apec_m_visibility.h"
#include "augment_sparsity_on_nodes.h"

class vtkUnstructuredGrid;
class ApecMEquationSystems;

namespace libMesh
{
    class LibMeshInit;
}

class APEC_M_HIDDEN ApecMModel
{
public:
    ApecMModel() = default;
    ~ApecMModel();

    void load_into_solver(const std::shared_ptr<libMesh::LibMeshInit> &libmesh_init,
                          pybind11::object apec_model_object);
    void solve();
    void write_system_to_vtk_grid(const std::string &system_name);
    void clear();

    bool is_loaded() const;
    pybind11::object get_apec_object() const;
    vtkUnstructuredGrid *get_active_vtk_grid() const;
    const AugmentSparsityOnNodes *get_coupling_functor() const;

private:
    std::unique_ptr<libMesh::DistributedMesh> mesh;
    std::shared_ptr<ApecMEquationSystems> equation_system;
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor;
    vtkUnstructuredGrid *active_vtk_grid = nullptr;
    pybind11::object apec_object;
    bool loaded = false;
    bool closed_ = false;
};

#endif // APEC_M_MODEL_H
