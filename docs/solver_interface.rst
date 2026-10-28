Solver Interface
================

The C++ solver does not know VTK, nor the ``pytec3d`` package. It receives the model as **views** of arrays owned
by the caller (pointer, type, shape: no copy) and the model data as a JSON object, and it writes the solution into
arrays of the caller. The Python package is one such caller; any other program, in Python or in C++, can drive the
solver the same way.

The two parts are separate distributions:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Distribution
     - Content
   * - ``pytec3d``
     - Pure Python (``pip install .``, needs ``numpy`` and ``vtk`` only): the model, the selections, the meshes, the
       analysis driver, the plasma algorithm. It imports the solver only when a model is loaded.
   * - ``pytec3d-solver``
     - The compiled module ``pytec3d_solver.pytec3d_solver``, built by the conda recipe against libMesh, PETSc and MPI
       (:doc:`installation`). It holds the solver core, a C++ library (``libpytec3d_core.a``) that depends neither on
       Python nor on VTK, and the `pybind11 <https://pybind11.readthedocs.io/>`_ binding that turns NumPy arrays into
       views and dicts into the C++ settings.

The module and the package check each other with an interface version (``pytec3d_solver.pytec3d_solver.INTERFACE_VERSION``,
:data:`pytec3d.solver.SOLVER_INTERFACE_VERSION`): a solver module built for another version is refused with a message
asking to rebuild it.

What the solver receives
------------------------

The model input (``solver/include/model_input.h`` in C++, the dict of :func:`pytec3d.solver_input.model_input` in Python):

Every array has a fixed type, by role: real numbers are float64 (``double``), ids int32, mesh indices int64, cell
types uint8.

.. list-table::
   :header-rows: 1
   :widths: 27 23 50

   * - Entry
     - Shape, type
     - Content
   * - ``mesh.points``
     - (n_points, 3) float64
     - Point coordinates.
   * - ``mesh.offsets``
     - (n_cells + 1,) int64
     - Cell ``i`` has the points ``connectivity[offsets[i]:offsets[i + 1]]`` (the VTK cell array layout).
   * - ``mesh.connectivity``
     - (offsets[-1],) int64
     - Point ids, in the node order of the VTK cell types.
   * - ``mesh.cell_types``
     - (n_cells,) uint8
     - VTK cell type numbers (``VTK_QUAD`` = 9, ``VTK_QUADRATIC_HEXAHEDRON`` = 25...): the convention, no VTK code.
   * - ``arrays.ids``
     - (n_cells,) int32
     - ``TYPE``, ``MATERIAL``, ``CONSTANT``, ``COORDINATE_SYSTEM``: the rows of the model data of each cell.
   * - ``arrays.couplings``
     - (n_points,) int32
     - The ``CP_<k>`` arrays (DOF couplings).
   * - ``arrays.point_fields``
     - (n_points, 1 or 3) float64
     - By name: each DOF field ``F`` and its Dirichlet values ``_F`` (NaN where free), the nodal loads (``AMPS``,
       ``CURR``...), the nodal fields (``B``, ``H``...).
   * - ``arrays.cell_fields``
     - (n_cells, 1 or 3) float64
     - By name: the element loads (``Js``...), the element fields (``B``, ``H``, ``Jt``...).
   * - ``element_rotations``
     - (n_cells, 3, 3) float64 or none
     - Local-to-global rotation of the ``COORDINATE_SYSTEM`` of each cell at its vertex average (the mean of its corner
       points); none for identity.
   * - ``node_rotations``
     - (n_points, 9) float64 or none
     - Local-to-global rotation of the DOF frame of each point (``node_angles``); none for identity.
   * - ``model_data``
     - JSON object
     - The sections of :doc:`user_manual/model_data` (element types, materials, constants, couplings, laws).

**Only the needed fields.** The fields are exactly those of the element types named in the model data, 3 components
for a vector field and 1 for a scalar: the solver derives the list from the model data (``required_field_arrays()``,
``solver/include/array_io.h``), checks the arrays against it when the model is loaded and at every ``set_arrays``, and stops
naming any field that is missing or has the wrong shape; a field it does not use is ignored with a warning. The Python
package hands over the same list (:func:`pytec3d.element_library.grid_arrays`, through
:func:`pytec3d.solver_input.model_arrays`), which :meth:`pytec3d.Model.update_grid` creates in the mesh.

**No conversion.** Arrays must be contiguous (C order), in the native byte order and of the type of their role; an
array with several components is (n, components) or (n, 3, 3), row-major. The Python binding refuses any other type
rather than converting it (a converted array would be a copy, cut off from the mesh). :meth:`pytec3d.Model.update_grid`,
called when the model is loaded, stores the mesh with these types (:func:`pytec3d.mesh.normalize_types`): a mesh read
with float32 points, 32-bit cells or ids of another integer type is converted once, there.

**Lifetime.** The solver keeps the views of the model arrays, not the arrays: each must stay alive and keep its
memory from the load to the end of the session, or until it is replaced with ``set_arrays``; the Python binding
keeps a reference to each. The mesh and the rotations are only read by ``load_model``, which copies what it needs
into the libMesh mesh and the rotation caches of the solver: no reference to them is kept, so the rotation arrays
computed in Python are freed as soon as the load returns. :meth:`pytec3d.Solver.solve` hands the session views of the
current arrays of the mesh before every solve, so that an array replaced in the mesh after the load is the one the
solver uses.

**When it is read.** The mesh, the ids (``TYPE``, ``MATERIAL``, ``CONSTANT``, ``COORDINATE_SYSTEM``), the ``CP_<k>``
arrays, the rotations and the model data are read by ``load_model``. The loads and the Dirichlet values are read at
every ``solve``, and the output fields are written at its end. The time step and the solver settings are arguments
of ``solve``, so they can change from one solve to the next.

From Python
-----------

:class:`pytec3d.Solver` does all of this. The session can also be driven directly:

.. code-block:: python

   import pytec3d_solver.pytec3d_solver as backend
   from pytec3d import Model, solver_input

   model = Model.read_model("coil.vtu")
   model.update_grid()                                     # the arrays of the element types, with the solver's types
   session = backend.Session(1)
   session.load_model(solver_input.model_input(model))   # views of model.mesh, no copy

   session.solve({"time": 1e-3, "theta": 1.0, "transient_effects": True},
                 {"linear": {"ksp_type": "preonly", "pc_type": "lu"}})
   # model.mesh's arrays A, B, H... now hold the solution
   session.set_arrays(solver_input.model_arrays(model))  # after replacing arrays of the mesh
   session.close()

From C++
--------

Link the core (CMake target ``pytec3d_core``, include ``pytec3d_session.h``) and hand it views of your own memory:

.. code-block:: cpp

   #include "pytec3d_session.h"

   pytec3d::ModelInput input;
   input.mesh.points = pytec3d::ConstRealArray(points.data(), n_points, 3);     // std::vector<double>
   input.mesh.offsets = pytec3d::IndexArray(offsets.data(), n_cells + 1);       // std::vector<std::int64_t>
   input.mesh.connectivity = pytec3d::IndexArray(connectivity.data(), connectivity.size());
   input.mesh.cell_types = pytec3d::CellTypeArray(cell_types.data(), n_cells);  // std::vector<std::uint8_t>
   input.arrays.type = pytec3d::IdArray(type.data(), n_cells);                  // std::vector<std::int32_t>
   input.arrays.material = pytec3d::IdArray(material.data(), n_cells);          // ... CONSTANT, COORDINATE_SYSTEM
   input.arrays.point_fields = {{"A", pytec3d::RealArray(a.data(), n_points, 3)},
                                {"_A", pytec3d::RealArray(a_fixed.data(), n_points, 3)} /* ... */};
   input.arrays.cell_fields = {{"Js", pytec3d::RealArray(js.data(), n_cells, 3)} /* ... */};
   input.model_data = nlohmann::json::parse(model_data_text);

   pyTEC3DSession session(1);
   session.load_model(input);
   session.solve(time_step, settings);   // TimeStepInput, SolverSettings: the solution is now in a
   session.close();
   pyTEC3DSession::shutdown_runtime();

``solver/examples/square_conductor_2d.cpp`` is a complete program: a 2D magnetostatic model built in memory, solved
twice (the load doubled in place in between). It is built with the solver when CMake is given
``-DPYTEC3D_BUILD_CPP_EXAMPLES=ON``, e.g. in an environment with the build dependencies of the solver (those of
``pytec3d-solver`` in ``conda/recipe.yaml``, ``pytec3d-libmesh`` included):

.. code-block:: bash

   cmake -S solver -B build/solver -G Ninja -DLIBMESH_DIR=$CONDA_PREFIX -DPYTEC3D_BUILD_CPP_EXAMPLES=ON
   cmake --build build/solver && build/solver/square_conductor_2d

It runs in serial or with ``mpiexec``.
