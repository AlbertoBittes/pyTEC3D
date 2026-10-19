Boundary Conditions
===================

Dirichlet values
----------------

Every DOF field ``F`` has a point array ``_F`` of the same shape (``_A``: 3 components; ``_V``, ``_VOLT``,
``_EMF``: 1). At every solve:

* ``NaN`` (the default) leaves the component free;
* a number fixes it to that value.

.. code-block:: python

   selection.select_all(mesh)
   selection.nsel_ext(mesh)
   selection.set_value_on_selection(mesh, "_A", 0.0, 2)        # A_z = 0 on the boundary
   selection.set_value_on_selection(mesh, "_A", np.nan, 2)     # free again

A boundary without any condition is a **natural** boundary: for A it means :math:`\mathbf{n}\times\mathbf{H}=0`
(the flux crosses the boundary normally, a symmetry plane of the type "B normal"); for V it means no current through
the boundary (an insulated surface).

Dirichlet values can change between steps. Changing which components are fixed (NaN <-> number) changes the set of
constrained DOFs and rebuilds the system matrix; changing only the values reuses it.

Conditions on the vector potential
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Situation
     - Condition
     - How (nodal A)
   * - Far boundary of an open domain (truncation)
     - :math:`\mathbf{A}=0` (B tangent to the boundary)
     - all components of ``_A`` = 0
   * - Symmetry plane, current **normal** to it (B tangent)
     - :math:`\mathbf{n}\times\mathbf{A}=0`
     - the tangential components = 0
   * - Symmetry plane, current **parallel** to it (B normal)
     - natural
     - nothing
   * - Axis of a 3D sector about an axis
     - components normal to the axis = 0
     - e.g. ``_A`` x and z for the y axis
   * - 2D model, outer boundary / axis
     - :math:`A_z=0` / :math:`A_\theta=0`
     - component 2 of ``_A``

In 2D the only DOF is component 2 (A_z, or A_theta in axisymmetry): ``A_z = 0`` on a boundary makes it a flux line
(B tangent); a free boundary is a symmetry plane with B normal to it.

Planes not aligned with the axes: rotated nodal frames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The components of ``_A`` refer to the **nodal frame** of each node, defined by the point array ``node_angles``: three
angles in degrees (``thxy``, ``thyz``, ``thzx``: rotations about z, then about the new x, then about the new y),
zero by default (global axes). The rotations are computed **when the model is loaded**: set ``node_angles`` before
``load_model_into_solver()`` (a later change stops the solver with an error, instead of being silently ignored). On a
plane at 30 degrees, rotate the frames of its nodes so that one local axis is normal to the plane, then fix the local
tangential components:

.. code-block:: python

   cyl = YAxisCylindricalCoordinateSystem()
   selection.select_pos(mesh, 1, (30 - 1e-5, 30 + 1e-5), coordinate_system=cyl)    # nodes at theta = 30 deg
   selection.set_value_on_selection(mesh, "node_angles", 30.0, 2)                   # rotate by thzx = 30 deg
   # local x = radial, local y = global y, local z = normal to the plane:
   selection.set_value_on_selection(mesh, "_A", 0.0, 0)
   selection.set_value_on_selection(mesh, "_A", 0.0, 1)

The rotation of a set of angles can be checked in Python:
``CoordinateSystem.build_local_to_global_matrix(thxy, thyz, thzx)`` returns the matrix whose columns are the local
axes in global coordinates. The nodal frames also apply to the DOF couplings (cyclic symmetry, :doc:`couplings`).
The solution arrays are always written in global components.

Edge elements
~~~~~~~~~~~~~

With ``AV3D_edge*`` element types, the DOFs of A are the circulations along the edges, stored on the **midside
nodes** in **component 2** of ``_A``. Fixing it on the midside nodes of a face fixes the tangential component of A on
that face, :math:`\mathbf{n}\times\mathbf{A}=0` (far field and "B tangent" symmetry planes). The values of ``_A`` on
corner nodes and the other components are not used. No nodal frame rotation is needed. Restrict a selection to the
midside nodes with ``selection.select(mesh, "MidsideNodes", 1, mode="r")``.

Gauging
-------

The vector potential is defined up to a gradient, :math:`\mathbf{A}+\nabla\psi` gives the same B. In the conductors
the eddy-current term fixes it; in the non-conducting regions the curl-curl matrix is singular. pyTEC3D handles this
as follows (:ref:`theory-gauge`):

* **nodal elements** (``AV2D_*``, ``AV3D_nodal*``): the 3D nodal types include a Coulomb-gauge penalty term
  :math:`\int\nu_g(\nabla\cdot\mathbf{w})(\nabla\cdot\mathbf{A})`, nothing to do (2D has no gauge freedom);
* **edge elements**: build a tree-cotree gauge with :func:`pytec3d.tree_gauge.tree_gauge` after the Dirichlet values
  of A are set. It builds a spanning tree of the edges of the **selected cells** (the region to gauge), grown from the
  **selected nodes** (the start surface, usually the far-field boundary where A is already fixed), and sets ``_A = 0``
  on the midside nodes of the tree edges:

  .. code-block:: python

     selection.select_all(mesh)                 # cells: the region whose edges are gauged
     selection.nsel_ext(mesh)                   # nodes: the start surface (the boundary)
     selection.set_value_on_selection(mesh, "_A", 0.0, 2)
     tree_gauge.tree_gauge(mesh, periodic_point_coupling="CP_1")   # tell it about periodic pairs, if any

  The log reports the tree: ``num_tree`` must be positive; ``tree_identity_ok`` checks the graph. With periodic
  couplings, pass their CP array so that the slave edges are left to the periodic constraint.

Electric potentials
-------------------

* The potential V of every 3D conductor (``SOLID_CONDUCTOR``) is defined up to a
  constant: fix it on at least one node of each conductor (``_V = 0``), unless another condition already does (a
  terminal at a fixed potential, a symmetry plane with ``V = 0``). A floating potential makes the system singular; the
  direct solver often still solves it, but the factorization can fail or lose accuracy (enable the MUMPS null pivot
  check, :doc:`solver_settings`).
* Terminals: fix V on the nodes of a terminal face to impose a potential; fix it on two faces to impose a voltage.
* ``_V`` of a 2D ``CIRCUIT_COUPLED_CONDUCTOR`` is its terminal voltage (fixed on the master node of its coupling);
  a 2D ``SOLID_CONDUCTOR`` has no electric DOF.
* ``_VOLT`` applies to the VOLT field of the stranded coils only: it has **no effect** on conductors, whose field is
  V.
* With the ``_timint`` element types the electric Dirichlet values are time integrals ([V s]).

Constraints table of the log
----------------------------

Each step prints a *Constraints* table: the number of static-step constraints (DOFs without stiffness fixed to zero
on a static step), Dirichlet constraints and coupling constraints, by variable. Check it: a boundary condition that
does not appear there was written on nodes that do not carry that field.
