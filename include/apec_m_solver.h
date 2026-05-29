/**
 * @file apec_m_solver.h
 * @brief Top-level solver class and pybind11-facing API for the APEC-M engine.
 *
 * ApecMSolver manages the full solver lifetime:
 *   - Initialises libMesh and spdlog via the process-singleton RuntimeManager.
 *   - Imports a VTK unstructured-grid mesh from the Python-side APEC model object.
 *   - Distributes the mesh and DOFs across MPI ranks.
 *   - Drives per-time-step solves and writes results back to the VTK grid.
 *
 * The class is exposed to Python via pybind11 in main.cpp.
 */
#ifndef APEC_M_SOLVER_H
#define APEC_M_SOLVER_H

#include <memory>
#include <string>
#include <mpi.h>
#include <pybind11/pybind11.h>

#include "apec_m_visibility.h"

namespace libMesh
{
    class LibMeshInit;
    class DistributedMesh;
}
class vtkUnstructuredGrid;

class ApecMEquationSystems;
class AugmentSparsityOnNodes;
class ElementTypesFactory;

/**
 * @brief Top-level solver that owns the libMesh mesh, equation systems, and VTK grid.
 *
 * Constructed once per Python simulation object.  Call load_model() to import
 * the mesh and initialise DOFs, then call solve() for each time step.
 */
class APEC_M_HIDDEN ApecMSolver
{
public:
    /**
     * @brief Constructs the solver and initialises the libMesh / MPI runtime.
     * @param mpi_comm  An mpi4py communicator, or @c py::none() for MPI_COMM_WORLD.
     */
    explicit ApecMSolver(pybind11::object mpi_comm = pybind11::none());

    /** @brief Calls close() to free all libMesh and PETSc resources. */
    ~ApecMSolver();

    /**
     * @brief Imports the model from a Python APEC object, distributes mesh and DOFs.
     * @param apec_object  Python object exposing @c mesh, @c element_types,
     *                     @c materials, @c constants, and @c time_step attributes.
     */
    void load_model(pybind11::object apec_object);

    /** @brief Releases all libMesh, PETSc, and VTK resources.  Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-level libMesh / spdlog runtime.
     *
     * Must be called after all ApecMSolver instances have been closed.  After
     * shutdown the runtime cannot be re-initialised in the same process.
     */
    static void shutdown_runtime();

    /** @brief Executes one time-step solve and writes results back to the VTK grid. */
    void solve();

    /** @brief Localises a named system's solution onto the active VTK grid. */
    void write_system_to_vtk_grid(const std::string &system_name);

    /** @brief Returns true if load_model() has completed successfully. */
    bool is_loaded() const;

    /** @brief Returns the Python APEC model object provided to load_model(). */
    pybind11::object get_apec_object() const;

    /** @brief Returns the active VTK unstructured-grid pointer (non-owning). */
    vtkUnstructuredGrid *get_active_vtk_grid() const;

    /** @brief Returns the element-types factory built from the model description. */
    std::shared_ptr<ElementTypesFactory> get_element_types_factory() const;

    /** @brief Returns the ghosting functor used for CP-node coupling (non-owning). */
    const AugmentSparsityOnNodes *get_coupling_functor() const;

private:
    int argc = 0;                                               ///< Placeholder argc forwarded to libMesh initialisation.
    char **argv = nullptr;                                      ///< Placeholder argv forwarded to libMesh initialisation.
    std::shared_ptr<libMesh::LibMeshInit> libmesh_init;         ///< Shared handle to the process-level libMesh init.
    std::unique_ptr<libMesh::DistributedMesh> mesh;             ///< Distributed libMesh mesh.
    std::shared_ptr<ApecMEquationSystems> equation_system;      ///< Equation systems container.
    std::shared_ptr<ElementTypesFactory> element_types_factory; ///< Per-subdomain element-type registry.
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor;   ///< Ghosting functor for CP node coupling.
    vtkUnstructuredGrid *active_vtk_grid = nullptr;             ///< Non-owning pointer to the Python-side VTK grid.
    pybind11::object apec_object;                               ///< Owning reference to the Python APEC model object.
    pybind11::object mpi_comm_object;                           ///< Owning reference to the mpi4py communicator (or None).
    bool loaded = false;                                        ///< Set to true after load_model() completes.
    bool closed_ = false;                                       ///< Set to true after close() is called.

    /** @brief Initialises the libMesh runtime on @p mpi_comm if not already initialised. */
    void ensure_libmesh_initialized(MPI_Comm mpi_comm);
};

#endif // APEC_M_SOLVER_H
