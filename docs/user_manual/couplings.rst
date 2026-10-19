DOF Couplings and Periodic Boundaries
=====================================

A DOF coupling makes a group of nodes share the value of some fields. It is used for:

* **terminals**: one V for the whole cross-section of a 2D circuit-coupled conductor, one VOLT and one EMF for a
  stranded coil, one V for an equipotential terminal face of a 3D conductor;
* **cyclic (periodic) symmetry**: the field on one sector face equals the field on the other face;
* any other "same value" condition between nodes.

How couplings are defined
-------------------------

Two ingredients:

#. a **point array** ``CP_k`` (``CP_0`` always exists; add ``CP_1``, ``CP_2``... when a node must belong to more than
   one coupling): the nodes with the same positive integer id in that array form a group;
#. the model data section ``dof_couplings``: the fields each id couples, ``{"<id>": ["VOLT", "EMF"]}``.

.. code-block:: python

   selection.select(mesh, "MATERIAL", 23, location="elem")
   selection.nsle(mesh)
   selection.set_value_on_selection(mesh, "CP_0", 1)          # group 1: the nodes of the coil
   solver.model.add_dof_couplings({1: ["VOLT", "EMF"]})       # or in model_data.json

In each group the node with the **lowest id is the master**; the other nodes are constrained to its values. Put the
loads and Dirichlet values of a terminal on the master node:

.. code-block:: python

   master = selection.selected_ids(mesh, "node")[0]           # the lowest id of the selection

Rules:

* every positive id of a ``CP_k`` array must have an entry in ``dof_couplings``;
* the fields must be DOF fields of the element type; a field that does not exist on a node of the group is simply
  not coupled there (e.g. ``VOLT`` on the nodes of a conductor, whose field is ``V``);
* groups are built when the model is loaded: set the ``CP_k`` arrays (and ``node_angles``) before
  ``load_model_into_solver()``; the solver stops with an error if they change afterwards;
* ids must be unique across the ``CP_k`` arrays of the model (``periodic_boundaries`` takes the next free ones).

Vector fields are coupled component by component **in the nodal frames** of the two nodes (``node_angles``): with the
frames of the second sector face rotated by the sector angle, the coupled components are the same cylindrical
components on both faces.

Periodic (cyclic) boundaries
----------------------------

:meth:`pytec3d.model.Model.periodic_boundaries` pairs the selected nodes of one face with the nodes of the image face:

.. code-block:: python

   from pytec3d import mesh as meshes, YAxisCylindricalCoordinateSystem
   import vtk

   cyl = YAxisCylindricalCoordinateSystem()
   meshes.remove_array(mesh, "CP_1")                               # a fresh array for the pairs
   meshes.add_array(mesh, "CP_1", dtype=vtk.VTK_INT)

   selection.select_pos(mesh, 1, (-1e-5, 1e-5), coordinate_system=cyl)        # nodes of the face theta = 0
   pairs = solver.model.periodic_boundaries(cyl, ("A", "V"), dy=30.0, array_name="CP_1")

   # Nodal vector potential: rotate the frames of the image face by the sector angle
   selection.select_pos(mesh, 1, (30 - 1e-5, 30 + 1e-5), coordinate_system=cyl)
   selection.set_value_on_selection(mesh, "node_angles", 30.0, 2)

Each selected node is moved by ``(dx, dy, dz)`` in the coordinates of the coordinate system (``dy = 30`` turns by 30
degrees in a cylindrical system; a Cartesian system gives a translation), matched with the nearest node (within
``tolerance``), and both nodes get a new CP id in ``array_name``; the ids are added to ``dof_couplings`` with the
given fields. Nodes already coupled in that array, and nodes without an image, are skipped (the log reports how
many). Returns ``{cp_id: fields}``.

Guidelines:

* select the whole face, including its corner and midside nodes (edge elements are coupled through their midside
  nodes; the relative orientation of the edges is handled by the solver);
* leave out the nodes on the axis if they have their own conditions, and avoid coupling nodes that are also fixed
  by a Dirichlet value with a different value;
* with **edge elements** no frame rotation is needed (edge values are scalars), but pass the CP array to
  ``tree_gauge(mesh, periodic_point_coupling="CP_1")`` so that the gauge tree accounts for the identified edges;
* couple ``V`` for conductors crossing the faces, ``A`` always.

Examples: ``tests/coil_3D_avoltemf`` (30 degrees, nodal), ``tests/coil_3D_avoltemf_edge`` (edge),
``tests/ITER_coil_3D_period`` (90 degrees), ``tests/ITER_3D_evol`` (2 degrees).
