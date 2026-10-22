The Mesh
========

A pyTEC3D model lives on a `VTK <https://vtk.org/>`_ [Schroeder2006]_ **unstructured grid**
(:vtk:`vtkUnstructuredGrid`), stored in ``.vtu`` files (or ``.pvtu`` for a partitioned grid), the `VTK XML format
<https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html>`_. Everything that varies over the model is
an **array** of the grid, on the points (``location="node"``, point data) or on the cells (``location="elem"``, cell
data).

Making a mesh
-------------

Any mesher can be used, as long as the result is converted to a ``.vtu`` file with coordinates in metres:

* **ANSYS** ``.cdb`` files: ``pytec3d-convert-cdb mesh.cdb`` writes ``mesh.vtu`` next to it (with
  `ansys-mapdl-reader <https://reader.docs.pyansys.com/>`_, installed by ``build.sh``). The ANSYS element type,
  material and real constant numbers become the cell arrays ``TYPE``, ``MATERIAL`` and ``CONSTANT``, and the arrays of
  a model are added. Node and element numbers of ANSYS are kept as ``ansys_node_num`` / ``ansys_elem_num``.
* **Gmsh** (https://gmsh.info/), **Salome**, **Cubit**...: export to VTK, or convert with `meshio
  <https://github.com/nschloe/meshio>`_ (``meshio convert mesh.msh mesh.vtu``); physical groups become integer arrays
  that can be turned into ``TYPE`` / ``MATERIAL`` ids with :mod:`pytec3d.selection`.
* **ParaView / VTK**: any filter pipeline that ends in an unstructured grid can be saved as ``.vtu``.
* **Python**: build a ``vtkUnstructuredGrid`` and use :meth:`pytec3d.solver.Solver.from_grid` / :meth:`pytec3d.model.Model.from_grid`.

Supported cells
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Dimension
     - Linear
     - Quadratic
   * - 2D
     - ``VTK_TRIANGLE`` (TRI3), ``VTK_QUAD`` (QUAD4)
     - ``VTK_QUADRATIC_TRIANGLE`` (TRI6), ``VTK_QUADRATIC_QUAD`` (QUAD8), ``VTK_BIQUADRATIC_QUAD`` (QUAD9)
   * - 3D
     - ``VTK_TETRA``, ``VTK_HEXAHEDRON``, ``VTK_WEDGE``, ``VTK_PYRAMID``
     - ``VTK_QUADRATIC_TETRA`` (TET10), ``VTK_QUADRATIC_HEXAHEDRON`` (HEX20), ``VTK_TRIQUADRATIC_HEXAHEDRON`` (HEX27),
       ``VTK_QUADRATIC_WEDGE`` (PRISM15), ``VTK_BIQUADRATIC_QUADRATIC_WEDGE`` (PRISM18), ``VTK_QUADRATIC_PYRAMID``
       (PYRAMID13)

Lines and some Bezier cells are also read. Rules:

* All the cells of a model have the dimension of its element type: 2D cells for a 2D element type, 3D cells for a 3D
  one (no mixed 2D / 3D meshes).
* An FE order of 2 (``option2 = 2``) needs quadratic cells.
* **Edge elements** (``AV3D_edge*``) need quadratic cells: their DOFs are stored on the midside nodes, one per edge.
* 2D axisymmetric models use the x coordinate as the radius (``x >= 0``) and the y axis as the axis of symmetry.
* Cells whose node ordering is inverted are reoriented on import (the log counts them); degenerate and distorted
  cells are reported in the *Mesh quality* table.

The arrays of a model
---------------------

Cell arrays
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Array
     - Default
     - Meaning
   * - ``TYPE``
     - 1
     - Id of the row of ``data["element_types"]`` of the cell: its element type, domain, FE order. Also the libMesh
       subdomain id.
   * - ``MATERIAL``
     - 1
     - Id of the row of ``data["materials"]`` (and of the nonlinear laws of source ``mat``).
   * - ``CONSTANT``
     - 1
     - Id of the row of ``data["constants"]`` (real constants).
   * - ``COORDINATE_SYSTEM``
     - 0
     - Element coordinate system: 0 = global Cartesian, otherwise an id of ``data["coordinate_systems"]``. The vector
       and tensor parameters of the cell (materials, real constants, element loads) are given in it
       (:doc:`coordinate_systems`).
   * - ``Js``
     - 0
     - Element load: imposed source current density [A/m^2] of the stranded conductors (3 components).
   * - ``B``, ``H``, ``Jt``
     - 0
     - Element results: mean values over the cell (see :doc:`outputs`).
   * - ``cell_centroids``, ``cell_volume``
     - computed
     - Centroid of each cell and its area (2D) or volume (3D); used by the selection by position and by the plasma
       algorithm.
   * - ``selection_cell``
     - 0
     - The current cell selection (1 = selected).

Point arrays
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Array
     - Default
     - Meaning
   * - ``A`` / ``_A``
     - 0 / NaN
     - Magnetic vector potential [Wb/m] (3 components) and its Dirichlet values. In 2D only component 2 (A_z) is a
       DOF; with edge elements the edge values are in component 2 of ``_A`` on the midside nodes.
   * - ``V`` / ``_V``
     - 0 / NaN
     - Electric potential [V] (or its time integral [V s]): of the 3D solid conductors, and the terminal voltage of
       the 2D circuit-coupled conductors.
   * - ``VOLT`` / ``_VOLT``
     - 0 / NaN
     - Terminal voltage of a stranded coil [V] (or its time integral).
   * - ``EMF`` / ``_EMF``
     - 0 / NaN
     - Induced voltage of a stranded coil [V] (or its time integral, the linked flux [Wb]).
   * - ``CURR``
     - 0
     - Nodal load of ``VOLT``: the current of a stranded coil [A].
   * - ``AMPS``
     - 0
     - Nodal load of ``V``: current injected into a 3D conductor at the node, or total current of a 2D circuit-coupled
       conductor (on its master node) [A].
   * - ``B``, ``H``
     - 0
     - Nodal results: flux density [T] and field [A/m], projected on the nodes.
   * - ``CP_0``, ``CP_1``, ...
     - 0
     - DOF coupling ids (:doc:`couplings`). ``CP_0`` always exists; more arrays allow a node to belong to several
       couplings.
   * - ``node_angles``
     - 0
     - Orientation of the nodal frame (thxy, thyz, thzx, degrees) in which the Dirichlet values and the couplings of
       the node are given (:doc:`boundary_conditions`).
   * - ``MidsideNodes``
     - computed
     - 1 on the midside nodes of quadratic cells (computed by ``read_model``).
   * - ``selection_point``
     - 0
     - The current node selection.

Which field arrays a model has follows from the element types named in its model data
(:func:`pytec3d.element_library.grid_arrays`): a vector field has 3 components and a scalar one 1, whatever the
dimension of the model. Other arrays of the mesh (named selections such as ``X_BOUND``, component flags such as
``COILS``, results of other tools) are kept untouched and can be used in selections.

Working with the arrays in Python
---------------------------------

:mod:`pytec3d.mesh` gives NumPy access to the arrays:

.. code-block:: python

   from pytec3d import mesh as meshes

   grid = solver.model.mesh
   jt = meshes.array_values(grid, "Jt", "elem")      # NumPy VIEW, shape (cells, 3): writing into it changes the mesh
   xyz = meshes.points(grid)                          # point coordinates, shape (points, 3)
   meshes.add_array(grid, "MY_FLAG", "elem", fill=0)  # add an array (no-op if it exists)
   meshes.remove_array(grid, "MY_FLAG")

After writing into a view directly, call ``Modified()`` on the VTK array (``meshes.get_array(grid, name,
location).Modified()``) so that VTK knows it changed; the functions of the package do it for you. Point and cell
ids are the 0-based indices of the grid; the solver uses the same numbering for its nodes and elements, so the ids
printed in the log (worst cell of the mesh quality table, constrained nodes...) are grid ids.

The model data inside the mesh
------------------------------

``Model.save_model()`` (and ``Solver.save_model()`` / ``write_output()``) store the model name and the model data, as
JSON text, in the field data of the grid (``pyTEC3D_META_MESH_NAME``, ``pyTEC3D_MODEL_DATA``). A saved model is
therefore self-contained: ``Solver.read_model("saved.vtu")`` restores both. Field arrays of older tools found in a
mesh.
