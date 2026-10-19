Selections
==========

Loads, boundary conditions, couplings, element types and materials are all values of mesh arrays on groups of nodes
or cells. :mod:`pytec3d.selection` builds these groups in the style of the APDL selection commands (``NSEL``,
``ESEL``, ``NSLE``, ``ESLN``) and writes values on them.

The current selection is stored on the mesh, in the arrays ``selection_point`` and ``selection_cell`` (1 = selected):
it persists between calls, can be seen in ParaView, and is saved with the model. ``read_model()`` clears it.

The pattern
-----------

.. code-block:: python

   from pytec3d import selection

   mesh = solver.model.mesh
   selection.select(mesh, "MATERIAL", 23, location="elem")      # 1. select cells
   selection.nsle(mesh)                                         # 2. select their nodes
   selection.set_value_on_selection(mesh, "CP_0", 1)            # 3. write a value on the selected nodes

Every function takes the ``location``: ``"node"`` (default) or ``"elem"``. **Forgetting** ``location="elem"`` when
selecting by a cell array (``TYPE``, ``MATERIAL``, ``CONSTANT``...) or writing a cell array (``Js``, ``TYPE``...) is the
most common mistake: the function then looks for a point array of that name and raises
:class:`~pytec3d.errors.MissingArrayError`.

Selection modes
---------------

Every selecting function takes a ``mode``:

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Mode
     - APDL
     - Effect
   * - ``"s"``
     - ``S`` (default)
     - Replace the selection.
   * - ``"a"``
     - ``A``
     - Add to the selection.
   * - ``"r"``
     - ``R``
     - Reselect within the selection (intersection).
   * - ``"u"``
     - ``U``
     - Remove from the selection.

Functions
---------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Selects
   * - ``select(mesh, array_name, value, component=0, location="node", mode="s")``
     - Nodes / cells whose array value matches ``value``: a number (equality), a tuple ``(low, high)`` (inclusive
       range), or a list / set / 1D array of values. ``array_name="ID"`` selects by node / cell index;
       ``component="magnitude"`` compares the norm of a vector array.
   * - ``select_pos(mesh, axis, value, mode="s", coordinate_system=None, location="node")``
     - By coordinate (cells by centroid): ``axis`` 0, 1, 2 = x, y, z, or the coordinates of ``coordinate_system``
       (e.g. r, theta [degrees], axial for a cylindrical system; r, theta, phi for a spherical one). ``value`` is
       usually a range ``(low, high)``: use a tolerance, positions are never exact.
   * - ``select_all(mesh, location=None)``, ``unselect_all``, ``inv``
     - Everything, nothing, the complement; ``location=None`` acts on nodes and cells.
   * - ``nsle(mesh, mode="s")``
     - The nodes of the selected cells (APDL ``NSLE``).
   * - ``esln(mesh, all_nodes=True, mode="s")``
     - The cells whose nodes are selected: all of them (``all_nodes=True``) or any (APDL ``ESLN``).
   * - ``nsel_ext(mesh, mode="s")``
     - The nodes on the external boundary of the selected cells (surface in 3D, boundary edges in 2D). Selecting
       all the cells first gives the boundary of the model; selecting a region gives its boundary.
   * - ``apply_mask(mesh, mask, location="node", mode="s")``
     - Any boolean NumPy mask, combined with the selection by ``mode``.

Reading and writing:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Effect
   * - ``set_value_on_selection(mesh, array_name, value, component=0, location="node")``
     - Write ``value`` in one component of the array at the selected nodes / cells.
   * - ``get_selection_values(mesh, array_name, component=0, location="node")``
     - The values at the selected nodes / cells, in increasing id order (``"ID"`` gives the ids).
   * - ``selected_ids(mesh, location)``, ``selection_mask(mesh, location)``
     - The selected ids (increasing) / the boolean mask.
   * - ``external_node_adjacency(mesh)``
     - 2D only: the boundary nodes of the selected cells and their neighbours along the boundary (a graph).

Recipes
-------

Boundary of the whole model:

.. code-block:: python

   selection.select_all(mesh)
   selection.nsel_ext(mesh)
   selection.set_value_on_selection(mesh, "_A", 0, 2)

A face of a box, from a coordinate range (mind the tolerance):

.. code-block:: python

   selection.select_pos(mesh, 0, (0.5 - 1e-6, 0.5 + 1e-6))           # nodes with x = 0.5

The nodes on a cylinder or a sphere, or on a sector face:

.. code-block:: python

   from pytec3d import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
   cyl = YAxisCylindricalCoordinateSystem()
   selection.select_pos(mesh, 0, (1 - 1e-5, 1 + 1e-5), coordinate_system=cyl)                 # r = 1
   selection.select_pos(mesh, 1, (30 - 1e-5, 30 + 1e-5), coordinate_system=cyl)               # theta = 30 deg
   selection.select_pos(mesh, 0, (100 - 0.1, 100 + 0.1), coordinate_system=SphericalCoordinateSystem())

A region defined geometrically, then turned into a material and an element type (before the model is loaded):

.. code-block:: python

   selection.select_all(mesh)
   selection.select_pos(mesh, 0, (0.4, 0.6), location="elem")                 # inner window
   selection.select_pos(mesh, 1, (-0.6, 0.6), location="elem", mode="r")
   selection.inv(mesh)                                                          # everything else
   selection.select_pos(mesh, 0, (0.25, 0.75), location="elem", mode="r")     # within the outer rectangle
   selection.select_pos(mesh, 1, (-0.75, 0.75), location="elem", mode="r")
   selection.set_value_on_selection(mesh, "MATERIAL", 26, location="elem")
   selection.set_value_on_selection(mesh, "TYPE", 26, location="elem")

The midside nodes of a quadratic mesh only (edge DOFs), or only the corner nodes:

.. code-block:: python

   selection.select(mesh, "MidsideNodes", 1, mode="r")     # keep the midside nodes of the selection
   selection.select(mesh, "MidsideNodes", 0, mode="r")     # keep the corner nodes

The master node of a group (the lowest id) and a load on it:

.. code-block:: python

   selection.select(mesh, "MATERIAL", 23, location="elem")
   selection.nsle(mesh)
   master = selection.selected_ids(mesh, "node")[0]
   ...
   selection.select(mesh, "ID", master)
   selection.set_value_on_selection(mesh, "AMPS", 100.0)

Named arrays of the mesh (components exported by the mesher):

.. code-block:: python

   selection.select(mesh, "X_BOUNDS_1", 1)                                    # point array, 1 on a face
   selection.select(mesh, "COIL_NOCOND_1", 1, location="elem")               # cell array, 1 on a group

Selections in ParaView
----------------------

The pyTEC3D ParaView plugin provides the same operations as filters (*Select*, *Select Extracted Data*, *Set Value On
Selection*): a selection can be made graphically and saved with the mesh, then used from a script with
``selection.select(mesh, "selection_point", 1)``. See :doc:`paraview`.
