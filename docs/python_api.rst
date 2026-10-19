Python API
==========

The reference below is generated from the docstrings of the package (``sphinx.ext.autosummary`` and ``autodoc``), so
it follows the code: rebuilding the documentation updates it.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Module
     - Content
   * - :mod:`pytec3d`
     - The main objects: :class:`~pytec3d.solver.Solver`, :class:`~pytec3d.model.Model`, :class:`~pytec3d.algorithm.Algorithm`, the
       coordinate systems, the exceptions, ``pytec3d.set_log_level()``.
   * - :mod:`pytec3d.solver`
     - :class:`~pytec3d.solver.Solver`: model files, solver session, time steps, outputs, run log; the default
       settings ``DEFAULT_SOLVER_SETTINGS`` and ``DEFAULT_TIME_STEP``.
   * - :mod:`pytec3d.model`
     - :class:`~pytec3d.model.Model`: mesh and model data, JSON input / output, DOF couplings, periodic boundaries,
       nonlinear laws, coordinate systems.
   * - :mod:`pytec3d.selection`
     - Selection of nodes and cells, values on the selection.
   * - :mod:`pytec3d.mesh`
     - Reading and writing meshes, the arrays of a model, NumPy access, geometry (centroids, sizes).
   * - :mod:`pytec3d.model_edit`
     - Coordinate mappings, vector fields in local components, derived fields (poloidal flux).
   * - :mod:`pytec3d.coordinate_system`
     - Cartesian, cylindrical and spherical coordinate systems.
   * - :mod:`pytec3d.tree_gauge`
     - Tree-cotree gauge of edge fields.
   * - :mod:`pytec3d.algorithm`
     - :class:`~pytec3d.algorithm.Algorithm`: plasma equilibrium and evolution.
   * - :mod:`pytec3d.element_library`
     - The element types, domains and fields of the JSON library, and the arrays they need.
   * - :mod:`pytec3d.validation`
     - The checks of the model, of the settings and of the time step.
   * - :mod:`pytec3d.errors`
     - The exceptions.
   * - :mod:`pytec3d.pytec3d_paraview_plugin`
     - ParaView filters and console helpers.
   * - :mod:`pytec3d.mesh_converter`
     - ANSYS ``.cdb`` to ``.vtu`` conversion (``pytec3d-convert-cdb``).

Usage
-----

.. code-block:: python

   from pytec3d import Solver, selection

   with Solver() as solver:
       solver.read_model("mesh.vtu")
       solver.model.read_json("model_data.json")
       selection.select(solver.model.mesh, "MATERIAL", 2, location="elem")
       selection.nsle(solver.model.mesh)
       selection.set_value_on_selection(solver.model.mesh, "CP_0", 1)
       solver.load_model_into_solver()
       solver.solve()
       solver.write_output()

Errors derive from :class:`pytec3d.errors.PyTEC3DError` (and from the matching built-in exception): an invalid model raises
:class:`pytec3d.errors.ValidationError` listing every problem, a missing array :class:`pytec3d.errors.MissingArrayError`, a change
to a model loaded into a solver session :class:`pytec3d.errors.ModelLockedError`. The log level is set with
``pytec3d.set_log_level("DEBUG")`` or the environment variable ``PYTEC3D_LOG_LEVEL``.

Reference
---------

The objects re-exported by :mod:`pytec3d` (``from pytec3d import Solver``...) are documented in their modules:

.. autosummary::
   :toctree: generated

   pytec3d.model
   pytec3d.solver
   pytec3d.selection
   pytec3d.mesh
   pytec3d.model_edit
   pytec3d.coordinate_system
   pytec3d.tree_gauge
   pytec3d.algorithm
   pytec3d.element_library
   pytec3d.validation
   pytec3d.errors
   pytec3d.pytec3d_paraview_plugin
   pytec3d.mesh_converter

The package
~~~~~~~~~~~

.. automodule:: pytec3d
   :no-members:
