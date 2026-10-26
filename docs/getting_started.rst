Getting Started
===============

Why pyTEC3D?
------------

Plasma disruptions induce large electromagnetic loads on tokamak components, so an accurate model of the plasma and
of the currents it induces in the surrounding structures is essential for their design. At Fusion for Energy this
analysis is done with `APEC <https://github.com/Fusion4Energy/APEC>`_ [APEC2018]_ [APEC2019]_ [APEC2020]_, a
finite-element code built on ANSYS APDL, which couples the axisymmetric plasma equilibrium and evolution equations
with the full 3D eddy-current equations of the passive structures, whereas most codes treat the whole problem as
axisymmetric.

The reliance on ANSYS brings performance issues (restarts of the iterations, repeated factorizations of the system
matrix) that limit the size of the models, and costly HPC licenses further restrict its use on large models.

pyTEC3D (python Tokamak Eddy-currents Code 3D) is an
independent, free, scalable and customizable code that reproduces the algorithm of APEC and removes these limitations.

What is pyTEC3D?
----------------

pyTEC3D (python Tokamak Eddy-currents Code 3D) solves time-domain, low-frequency (magneto-quasi-static) electromagnetic problems with the A, V-A formulation:
the magnetic vector potential **A** everywhere, an electric scalar potential **V** in the massive conductors (in 2D,
the terminal voltage of a conductor with a terminal), and the circuit unknowns of the coils (terminal voltage
**VOLT**, induced voltage **EMF**). It covers 2D planar, 2D axisymmetric and 3D models, nodal and edge (Nedelec) elements, linear and
saturating (B-H curve) materials, and static and transient analyses (see :doc:`theory/index`).

The code has two modules:

* a **Python package** (``pytec3d``) that reads the mesh, prepares the analysis (element types, materials, loads,
  boundary conditions, couplings), drives the time stepping and the plasma algorithm, and post-processes the results;
* a **C++ solver** (``pytec3d_solver``) that distributes the mesh and the degrees of freedom, builds the element
  matrices and the system in parallel with `libMesh <https://libmesh.github.io/>`_, and solves it with
  `PETSc <https://petsc.org/>`_ and `MUMPS <https://mumps-solver.org/>`_ (or an iterative solver). It depends neither
  on VTK nor on the Python package: it works on views of the arrays of the mesh, without copying them
  (:doc:`solver_interface`), and can be driven by any program.

pyTEC3D models are `VTK <https://vtk.org/>`_ unstructured grids (``.vtu``): a mesh made with any meshing tool can be
imported, and `ParaView <https://www.paraview.org/>`_ is used to inspect the model, prepare selections (with the
pyTEC3D plugin, see :doc:`user_manual/paraview`) and visualize the results.

Key concepts
------------

.. list-table::
   :widths: 22 78

   * - **Solver**
     - :class:`pytec3d.solver.Solver` is the analysis: it owns the model, the time step (``solver.time_step``), the solver
       settings (``solver.solver_settings``) and the C++ solver session. One per analysis.
   * - **Model**
     - :class:`pytec3d.model.Model`: the mesh (``model.mesh``, a ``vtkUnstructuredGrid``) and the model data
       (``model.data``, a JSON-like dict, usually read from ``model_data.json``).
   * - **Mesh arrays**
     - Everything that varies over the mesh is an array of the grid: on the cells the ids ``TYPE``, ``MATERIAL``,
       ``CONSTANT``, ``COORDINATE_SYSTEM`` and the element loads (``Js``); on the nodes the solution (``A``, ``V``...),
       the Dirichlet values (``_A``, ``_V``...: ``NaN`` = free), the nodal loads (``CURR``, ``AMPS``), the coupling ids
       (``CP_0``, ``CP_1``...) and the nodal frames (``node_angles``). See :doc:`user_manual/mesh`.
   * - **Element type**
     - A formulation of the library (``AV2D_nodal``, ``AV3D_nodal``, ``AV3D_edge`` and their time-integrated
       variants), chosen per ``TYPE`` id together with a **domain** (air, conductor, coil, magnet...), an FE order and
       the axisymmetric flag. See :doc:`user_manual/element_library`.
   * - **Selection**
     - APDL-like selection of nodes and cells (:mod:`pytec3d.selection`): select by array value, position, or
       topology, then set a value on the selection. This is how loads and boundary conditions are applied.
   * - **Session**
     - ``solver.load_model_into_solver()`` validates the model and creates the C++ session. From then on the mesh and
       the model data are locked; array values (loads, Dirichlet values) can still change between steps.
   * - **Time step**
     - ``solver.time_step["time"]`` is the end time of the next step; ``solver.solve()`` solves it (static or
       transient, linear or Newton) and writes the solution into the mesh arrays; ``solver.write_output()`` saves
       ``output/<name>_out_<n>.vtu``.

The workflow
------------

.. code-block:: text

   mesh (.vtu, from any mesher;  pytec3d-convert-cdb for ANSYS .cdb)
     |
     v
   Solver.read_model()           mesh + the arrays the element types need
   Model.read_json()             element types, materials, real constants, couplings, laws, coordinate systems
     |
     v
   selection + set_value_on_selection()   Dirichlet values (_A, _V ...), DOF couplings (CP_k), loads (Js, CURR ...)
     |
     v
   Solver.load_model_into_solver()        validation, C++ session (mesh distributed, DOFs numbered)
     |
     v
   loop over the time steps:  set time_step and loads -> Solver.solve() -> Solver.write_output()
     |
     v
   output/<name>_out_<n>.vtu, <name>_steps.csv, <name>.log   ->   ParaView, Python post-processing

A first analysis
----------------

The example below is the complete script of ``examples/coil_2D_avolt`` (with its mesh and model data):
two rings in 2D axisymmetry, the lower one a
massive conductor with a terminal (``CIRCUIT_COUPLED_CONDUCTOR``) driven first by its current and then by its voltage,
the upper one a passive ring (``SOLID_CONDUCTOR``) where eddy currents are induced. The :doc:`tutorial notebook <tutorials/tutorial_coil_2D>` runs the same analysis step by step
and looks at the results.

.. code-block:: python

   from pytec3d import selection, Solver
   import numpy as np

   with Solver() as solver:
       # Model: mesh + model data
       solver.read_model("coil_2D_avolt.vtu")
       solver.model.read_json("model_data.json")
       mesh = solver.model.mesh

       # Boundary condition: A_z = 0 on the external boundary
       selection.select_all(mesh)
       selection.nsel_ext(mesh)
       selection.set_value_on_selection(mesh, "_A", 0, 2)

       # Terminal of the conductor: one V for the whole ring (dof_couplings["1"] = ["V"])
       selection.select(mesh, "MATERIAL", 23, location="elem")
       selection.nsle(mesh)
       selection.set_value_on_selection(mesh, "CP_0", 1)
       master = selection.selected_ids(mesh, "node")[0]

       solver.load_model_into_solver()
       solver.time_step["theta"] = 1.0            # backward Euler: the loads change suddenly

       for i in (1, 2):
           solver.time_step["time"] = i * 1e-5
           selection.select(mesh, "ID", master)
           if i == 1:
               selection.set_value_on_selection(mesh, "AMPS", 100)     # 100 A through the ring
           else:
               selection.set_value_on_selection(mesh, "AMPS", 0)
               selection.set_value_on_selection(mesh, "_V", 1e4)       # then 10 kV across it
           solver.solve()
           solver.write_output()

       # Results are mesh arrays: e.g. the total current density of the cells (toroidal component)
       jt = selection.get_selection_values(mesh, "Jt", 2, location="elem")

Run it with ``python script.py`` in a folder that holds ``coil_2D_avolt.vtu`` and ``model_data.json``. It writes:

* ``output/coil_2D_avolt_out_1.vtu`` and ``_out_2.vtu``: the model with the solution of each step (open them in
  ParaView);
* ``output/coil_2D_avolt_steps.csv``: one row per step (time, iterations, convergence, residual, timing);
* ``coil_2D_avolt.log``: the run report (model tables, constraints, solver statistics, warnings).

Where next
----------

* :doc:`installation` -- build and install the code.
* :doc:`tutorials/index` -- Jupyter notebooks: complete analyses run step by step, with plots, including an ITER
  plasma equilibrium and a vertical displacement event.
* :doc:`user_manual/index` -- every step of an analysis in detail, and the element library.
* :doc:`examples/index` -- the annotated test cases: 2D and 3D coils, conductors, iron, cyclic symmetry, ITER
  equilibrium and vertical displacement.
* :doc:`theory/index` -- the equations, the discretization, the time integration and the Newton-Raphson solution.

Project links
-------------

* Source repository: https://github.com/AlbertoBittes/pyTEC3D
* Issues: the GitHub issue tracker of the repository.
* APEC, the ANSYS APDL reference of the plasma algorithm: https://github.com/Fusion4Energy/APEC
* Benchmarks and validation against ANSYS MAPDL and APEC: https://github.com/AlbertoBittes/pyTEC3D_APDL_benchmarks
  (:doc:`testing`).
