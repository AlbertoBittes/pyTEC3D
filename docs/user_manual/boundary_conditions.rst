Boundary Conditions
===================

Dirichlet values
----------------

Every DOF field ``F`` has a point array ``_F`` of the same shape (``_A``, ``_T``, ``_Hs``: 3 components; ``_V``,
``_VOLT``, ``_EMF``, ``_PHI``: 1; the edge fields use component 2 on the midside nodes). At every solve:

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
axes in global coordinates. The nodal frames also apply to the DOF couplings (cyclic symmetry, :doc:`couplings`) and
to the vector nodal loads (:doc:`loads`). The solution arrays are always written in global components.

A field with a single DOF component (A_z of a 2D element) ignores the nodal frames: its Dirichlet value is component 2
of ``_A`` (the component of the field) as given, and its couplings and nodal loads are not rotated either.

Edge elements
~~~~~~~~~~~~~

With ``AVA3D_edge*`` element types, the DOFs of A are the circulations along the edges, stored on the **midside
nodes** in **component 2** of ``_A``. Fixing it on the midside nodes of a face fixes the tangential component of A on
that face, :math:`\mathbf{n}\times\mathbf{A}=0` (far field and "B tangent" symmetry planes). The values of ``_A`` on
corner nodes and the other components are not used. No nodal frame rotation is needed. Restrict a selection to the
midside nodes with ``selection.select(mesh, "MidsideNodes", 1, mode="r")``.

Gauging
-------

The vector potential is defined up to a gradient, :math:`\mathbf{A}+\nabla\psi` gives the same B. In the conductors
the eddy-current term fixes it; in the non-conducting regions the curl-curl matrix is singular. pyTEC3D handles this
as follows (:ref:`theory-gauge`):

* **nodal elements** (``AVA2D_*``, ``AVA3D_nodal*``): the 3D nodal types include a Coulomb-gauge penalty term
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

Scalar potentials: PHI, T, Hs
------------------------------

The scalar-potential types (``PHI3D_nodal``, ``TPHI3D_edge``, ``TPHI2D_edge``; :ref:`el-phi3d`, :ref:`el-tphi3d`)
have :math:`\mathbf{H}=-\nabla\phi` outside the conductors (plus :math:`\mathbf{H}_s` or :math:`\mathbf{T}`), so the
conditions are the dual of those on A:

.. list-table::
   :header-rows: 1
   :widths: 32 30 38

   * - Situation
     - Condition
     - How
   * - Free outer boundary (truncation)
     - :math:`\mathbf{n}\cdot\mathbf{B}=0` (flux wall: B tangent)
     - nothing on ``_PHI`` (natural), but **fix** ``_PHI`` **on one node**: :math:`\phi` is defined up to a
       constant
   * - Applied field :math:`\mathbf{H}_0` (uniform, or any known exterior field)
     - the total potential on the boundary
     - ``_PHI`` = :math:`-\mathbf{H}_0\cdot\mathbf{x}` (+ the exterior potential of the model's sources when known),
       one value per boundary node (``set_value_on_selection`` with an array)
   * - Symmetry plane, B **tangent** to it
     - natural
     - nothing
   * - Symmetry plane, B **normal** to it
     - :math:`\phi` constant
     - ``_PHI`` = 0 (or the value of the applied potential) on the plane
   * - 2D axisymmetric axis
     - natural (:math:`B_r=0`)
     - nothing

*In the conductors* (T-phi): :math:`\mathbf{n}\times\mathbf{T}=0` on the conductor surface, ``_T`` = 0 on the
midside nodes of its surface edges (component 2), and a gauge inside: ``_PHI`` = 0 at the corner nodes strictly inside
the conductor (first-order :math:`\phi`: its gradients are in the edge space of T), or a tree gauge of T rooted on the
surface (``tree_gauge.tree_gauge(..., array_name="_T")`` over the conductor cells) with :math:`\phi` free. A conductor
with a net current needs a cut (:ref:`el-cuts`); a coil is ``_T`` fixed on the edges of its winding and hole
(:ref:`el-tphi3d`).

*Static steps* (T-phi): at rest :math:`\mathbf{H}=-\nabla\phi` in the conductors too. On a static step set ``_T`` = 0
on every conductor edge and ``_PHI`` = NaN at the inner nodes, then restore the transient gauges before the next
transient step (example ``bulk_2D_tphi_axisym``). The solver does not do it for you.

*The source field* (``PHI3D_nodal``): a tree gauge of ``_Hs`` over all the cells, rooted on the boundary, and
:math:`\mathbf{n}\times\mathbf{H}_s=0` on the boundary, always (also without sources):

.. code-block:: python

   selection.select_all(mesh)
   selection.nsel_ext(mesh)                                   # the boundary nodes: the roots of the tree
   tree_gauge.tree_gauge(mesh, array_name="_Hs")              # _Hs = 0 on the tree edges
   selection.set_value_on_selection(mesh, "_Hs", 0.0, 2)      # n x Hs = 0 on the boundary edges

or ``_Hs`` given on every edge (the circulations of a known source field, :func:`pytec3d.tree_gauge.edges`,
:func:`pytec3d.tree_gauge.edge_circulations`), which removes the Hs unknowns.

The value of an edge DOF is a circulation along the edge in the solver's orientation (from the end point with the
lexicographically smaller coordinates to the larger); :func:`pytec3d.tree_gauge.edges` returns each edge's midside
point and its end points in that orientation.

Constraints table of the log
----------------------------

Each step prints a *Constraints* table: the number of static-step constraints (DOFs without stiffness fixed to zero
on a static step), Dirichlet constraints and coupling constraints, by variable. Check it: a boundary condition that
does not appear there was written on nodes that do not carry that field.
