#ifndef APEC_M_SOLVER_H
#define APEC_M_SOLVER_H

#include <memory>
#include <string>
#include <mpi.h>
#include <pybind11/pybind11.h>

namespace libMesh
{
    class LibMeshInit;
    class DistributedMesh;
}
class vtkUnstructuredGrid;

class ApecMEquationSystems;
class AugmentSparsityOnNodes;
class ElementTypesFactory;

#if defined(__GNUC__) || defined(__clang__)
#define APEC_M_SOLVER_CLASS_HIDDEN __attribute__((visibility("hidden")))
#else
#define APEC_M_SOLVER_CLASS_HIDDEN
#endif

class APEC_M_SOLVER_CLASS_HIDDEN ApecMSolver
{
public:
    explicit ApecMSolver(pybind11::object mpi_comm = pybind11::none());
    ~ApecMSolver();

    void load_model(pybind11::object apec_object);
    void close();
    static void shutdown_runtime();

    void solve();
    void write_system_to_vtk_grid(const std::string &system_name);
    bool is_loaded() const;
    pybind11::object get_apec_object() const;

public:
    int argc = 0;
    char **argv = nullptr;
    std::shared_ptr<libMesh::LibMeshInit> libmesh_init;
    std::unique_ptr<libMesh::DistributedMesh> mesh;
    std::shared_ptr<ApecMEquationSystems> equation_system;
    std::shared_ptr<ElementTypesFactory> element_types_factory;
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor;
    vtkUnstructuredGrid *active_vtk_grid = nullptr;
    pybind11::object apec_object;
    pybind11::object mpi_comm_object;
    bool loaded = false;

private:
    void ensure_libmesh_initialized(MPI_Comm mpi_comm);
};

#undef APEC_M_SOLVER_CLASS_HIDDEN

#endif // APEC_M_SOLVER_H
