.. _element-library:

Element Library
===============

An **element type** is a complete finite-element formulation: its unknowns (DOF fields) and their finite elements,
the parameters it reads from the model data and the mesh, the terms of its weak form, and the fields it computes after
each solve. A **domain** (``option1``) selects, inside an element type, the physical region a cell represents: the
same element type gives air, a magnet, a coil or a massive conductor depending on the domain.

The library is defined in ``element_library/pytec3d_element_library/data/element_types.json`` (with ``fields.json``, ``domains.json`` and
``nonlinear_laws.json``) and compiled into the solver: this page describes its content. The equations are derived in
:doc:`../theory/index`: :doc:`../theory/domains` for the domains, :ref:`theory-2d` and :ref:`theory-3d` for the
formulations, :doc:`../theory/circuits` for coils and terminals; :ref:`theory-kernels` explains how the JSON terms map
to the weak form.

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 8 20 8 30 34

   * - Index
     - Name
     - Dim.
     - DOF fields (finite elements)
     - Use
   * - 0
     - ``AVA2D_nodal``
     - 2
     - A_z (Lagrange); V, VOLT, EMF (Lagrange, out of plane)
     - 2D planar or axisymmetric eddy currents (A, V-A), conductors with a terminal, stranded coils.
   * - 1
     - ``AVA3D_nodal``
     - 3
     - A (Lagrange vector); V, VOLT, EMF (Lagrange)
     - 3D eddy currents (A, V-A), nodal vector potential with a Coulomb gauge (penalty) term.
   * - 2
     - ``AVA3D_edge``
     - 3
     - A (Nedelec, first kind); V, VOLT, EMF (Lagrange)
     - 3D eddy currents (A, V-A) with edge elements: tangential continuity only, tree-cotree gauge.
   * - 3
     - ``AVA2D_nodal_timint``
     - 2
     - as 0, V, VOLT and EMF time-integrated
     - Symmetric coupling in 2D.
   * - 4
     - ``AVA3D_nodal_timint``
     - 3
     - as 1, V, VOLT and EMF time-integrated
     - Symmetric coupling, nodal.
   * - 5
     - ``AVA3D_edge_timint``
     - 3
     - as 2, V, VOLT and EMF time-integrated
     - Symmetric coupling, edge.
   * - 6
     - ``TPHI3D_edge``
     - 3
     - PHI (Lagrange), T (Nedelec)
     - 3D eddy currents (and magnetostatics) in conductors and superconductors, T-phi formulation (also T-Omega,
       H-phi): the scalar potential everywhere, the current vector potential T in the conductors only. No source
       current: a coil is T fixed on its edges.
   * - 7
     - ``PHI3D_nodal``
     - 3
     - PHI (Lagrange), Hs (Nedelec)
     - 3D magnetostatics with the scalar potential: H = Hs - grad(PHI), the source field Hs of the currents
       (curl Hs = Js) solved in the same system.
   * - 8
     - ``TPHI2D_edge``
     - 2
     - PHI (Lagrange), T (Nedelec, in-plane)
     - 2D T-phi (the 2D H-formulation of superconductor cross-sections): H in the plane, J and E out of plane;
       planar or axisymmetric.

The A, V-A types, 6 and 8 are first order in time (``time_order = 1``: the theta method, see
:ref:`theory-time-integration`); 7 is static (``time_order = 0``). A model uses one element type for all its cells, with
any of its domains. 2D magnetostatics uses ``AVA2D_nodal``. Theory of 6 and 7: :doc:`../theory/formulations_potentials`.

.. _el-domains:

Domains
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Domain (``option1``)
     - Physics
     - Fields and inputs
   * - ``NON_CONDUCTING`` (0)
     - Magnetostatics only, :math:`\mathbf{J}=0`: air, vacuum, insulators.
     - A. Material: ``mur``.
   * - ``SOFT_MAGNETIC`` (1)
     - As NON_CONDUCTING: iron without current, usually with a ``bh_curve`` law on its material.
     - A. Material: ``mur`` (or the B-H curve).
   * - ``HARD_MAGNETIC`` (2)
     - Permanent magnet: :math:`\mathbf{H}=\nu\mathbf{B}-\mathbf{M}`, the magnetization adds the source
       :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{M}` (T-phi: :math:`\int\nabla q\cdot\mathbf{B}_r`,
       :math:`\mathbf{B}_r=\mu_0\mu_r\mathbf{M}`). With a B-H curve see :doc:`../theory/constitutive`.
     - A. Material: ``mur``, ``M``.
   * - ``STRANDED_CONDUCTOR`` (3)
     - Imposed current density, no eddy currents: :math:`\mathbf{J}=\mathbf{J}_s`.
     - A. Material: ``mur`` (``M``). Element load ``Js``. ``PHI3D_nodal``: the source of Hs.
   * - ``CIRCUIT_COUPLED_CONDUCTOR`` (4)
     - **2D only.** Massive conductor with eddy currents and a **terminal**: the voltage V across its out-of-plane
       length, one value per conductor (coupled). Driven by its total current (``AMPS``) or by its voltage (``_V``).
     - A + V. Material: ``mur``, ``rho`` (``M``). Nodal load ``AMPS``.
   * - ``SOLID_CONDUCTOR`` (5)
     - Massive conductor with eddy currents. 2D: no terminal, :math:`J_z=-\sigma\partial_tA_z` (a closed ring, an
       infinitely long bar). 3D: with the electric potential V (charge conservation), currents injected with ``AMPS``
       and potentials fixed with ``_V``.
     - 2D: A. 3D: A + V, nodal load ``AMPS``. Material: ``mur``, ``rho`` (``M``).
   * - ``STRANDED_COIL`` (6)
     - Multi-turn winding with uniform current density :math:`\mathbf{J}=(N/S)\,i\,\mathbf{T}` and circuit equation
       :math:`\mathrm{VOLT}=R\,i+\mathrm{EMF}`, :math:`\mathrm{EMF}=d\Psi/dt`. No eddy currents.
     - A + VOLT + EMF (coupled per coil). Material: ``mur`` (``M``). Real constants ``cross_section``, ``turns``,
       ``resistance``, ``symmetry``, ``T``, and ``radius`` (2D) or ``volume`` (3D). Nodal load ``CURR``. Not in the
       scalar potential types (T-phi: a coil is T fixed on its edges, :ref:`el-tphi3d`).
   * - ``CUT_LAYER`` (7)
     - ``TPHI3D_edge``: cells with T but no conduction: a layer of one cell spanning the hole of a ring (its net
       current), or the hole of a coil given as T fixed on its edges (:ref:`el-tphi3d`).
     - PHI + T (collapsed in a cut, fixed in a coil). Material: ``mur``.

* The **magnetization** ``M`` is read in every domain except ``NON_CONDUCTING`` and ``SOFT_MAGNETIC`` (default 0): a
  magnetized region can also carry currents.
* ``VOLT`` and ``EMF`` belong to the stranded coils; the potential of the massive conductors is ``V``.
* ``CIRCUIT_COUPLED_CONDUCTOR`` is not a domain of the 3D element types: in 3D every solid conductor has V. Each
  element type lists its domains (``"domains"`` in the JSON); the model check rejects a row with another one.

Notation
~~~~~~~~

On this page :math:`\mathbf{w}` is the test function of A, :math:`q` the test function of the scalar field of the
row, :math:`\nu=(\mu_0\,\mathrm{mur})^{-1}` the reluctivity tensor, :math:`\sigma=\mathrm{rho}^{-1}` the conductivity
tensor, :math:`\partial_t` the time derivative, and the integrals run over the cells of the domain where the term is
defined. A term is listed with the matrix it goes into: **K** (stiffness), **C** (damping: terms with a time
derivative of the unknown) or **F** (right-hand side).

.. _el-av2d:

``AVA2D_nodal`` (index 0)
-------------------------

2D formulation, planar (``option3 = 0``) or axisymmetric about the y axis (``option3 = 1``). The vector potential has
only its out-of-plane component, :math:`\mathbf{A}=A_z\mathbf{e}_z` (planar) or :math:`A_\theta\mathbf{e}_\theta`
(axisymmetric), stored as component 2 of ``A``. Theory: :ref:`theory-2d`.

Out-of-plane length
   Every integral is over the 3D body: the integration weight of a 2D cell is multiplied by the out-of-plane length
   :math:`L=\mathrm{thickness}` (planar) or :math:`L=2\pi r\,\mathrm{thickness}` (axisymmetric), with ``thickness``
   the real constant (default 1). The terminal voltage V of a circuit-coupled conductor acts across that length: its
   gradient is :math:`\nabla V=-V/L\;\mathbf{e}_z`.

Fields
   .. list-table::
      :header-rows: 1
      :widths: 15 20 25 40

      * - Field
        - Finite elements
        - Domains
        - Notes
      * - ``A`` (A_z)
        - Lagrange, order ``option2``
        - all
        - Dirichlet values: component 2 of ``_A``.
      * - ``V``
        - Lagrange, out of plane
        - CIRCUIT_COUPLED_CONDUCTOR
        - Terminal voltage. Couple it over each conductor (one value). Nodal load: ``AMPS``.
      * - ``VOLT``
        - Lagrange, out of plane
        - STRANDED_COIL
        - Coil terminal voltage. Couple it with EMF over each coil. Nodal load: ``CURR``.
      * - ``EMF``
        - Lagrange, out of plane
        - STRANDED_COIL
        - Induced voltage of the coil.

Parameters
   .. list-table::
      :header-rows: 1
      :widths: 22 30 48

      * - Parameter
        - Domains
        - Meaning
      * - ``mat.mur`` (tensor)
        - all
        - Relative permeability.
      * - ``mat.rho`` (tensor)
        - SOLID_CONDUCTOR, CIRCUIT_COUPLED_CONDUCTOR
        - Resistivity [Ohm m].
      * - ``mat.M`` (vector, default 0)
        - all but NON_CONDUCTING, SOFT_MAGNETIC
        - Magnetization [A/m] (in-plane components).
      * - ``me.Js`` (element load)
        - STRANDED_CONDUCTOR
        - Source current density [A/m^2]; component 2 (out of plane / toroidal).
      * - ``mn.AMPS`` (nodal load)
        - CIRCUIT_COUPLED_CONDUCTOR
        - Load of the V equation: total current through the conductor [A], on the master node.
      * - ``mn.CURR`` (nodal load)
        - STRANDED_COIL
        - Load of the VOLT equation: coil current [A], on the master node.
      * - ``rc.thickness`` (default 1)
        - all
        - Out-of-plane length factor.
      * - ``rc.cross_section``, ``rc.turns``, ``rc.radius``, ``rc.resistance``, ``rc.symmetry``, ``rc.T``
        - STRANDED_COIL
        - Winding data (see :doc:`model_data`).

Derived coefficients (expressions)
   :math:`\nu=(\mu_0\,\mathrm{mur})^{-1}`, :math:`\sigma=\mathrm{rho}^{-1}`; for stranded coils
   :math:`n=N/(S\,R)`, coil volume :math:`V_c=(\text{axisym}\;?\;2\pi\,\mathrm{radius}:1)\,S\,\mathrm{thickness}`,
   :math:`k=1/(R\,V_c)`; the circuit row is scaled by ``symmetry`` (:ref:`theory-stranded-coil`, :doc:`model_data`).

Weak form
   .. list-table::
      :header-rows: 1
      :widths: 30 10 60

      * - Domain
        - Matrix
        - Term
      * - all
        - K
        - :math:`\int\nabla\times\mathbf{w}\cdot\nu\,\nabla\times\mathbf{A}`
      * - STRANDED_CONDUCTOR
        - F
        - :math:`\int\mathbf{w}\cdot\mathbf{J}_s`
      * - all but NON_CONDUCTING, SOFT_MAGNETIC
        - F
        - :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{M}`
      * - SOLID, CIRCUIT_COUPLED
        - C
        - :math:`\int\mathbf{w}\cdot\sigma\,\partial_t\mathbf{A}`
      * - CIRCUIT_COUPLED
        - K
        - :math:`\int\mathbf{w}\cdot\sigma\nabla V`, :math:`\int\nabla q\cdot\sigma\nabla V`
      * - CIRCUIT_COUPLED
        - C
        - :math:`\int\nabla q\cdot\sigma\,\partial_t\mathbf{A}`
      * - STRANDED_COIL
        - K
        - A row: :math:`-\int n\,\mathbf{w}\cdot\mathbf{T}\,(\mathrm{VOLT}-\mathrm{EMF})`; VOLT row:
          :math:`\int q\,s\,k\,\mathrm{VOLT}`; EMF row: :math:`-\int q\,k\,\mathrm{EMF}`
      * - STRANDED_COIL
        - C
        - VOLT row: :math:`-\int q\,s\,n\,\mathbf{T}\cdot\partial_t\mathbf{A}`; EMF row:
          :math:`\int q\,n\,\mathbf{T}\cdot\partial_t\mathbf{A}`

   with :math:`s` = ``symmetry``. With VOLT and EMF coupled over the coil these rows are the circuit equations
   :math:`\mathrm{EMF}=d\Psi/dt` and :math:`\mathrm{CURR}=(\mathrm{VOLT}-\mathrm{EMF})/R`
   (:ref:`theory-stranded-coil`); with V coupled over a circuit-coupled conductor its row states that the total
   current through the cross-section equals ``AMPS`` (:ref:`theory-solid-2d`).

Outputs
   .. list-table::
      :header-rows: 1
      :widths: 15 25 60

      * - Field
        - Where
        - Value
      * - ``B``
        - node and element
        - :math:`\nabla\times\mathbf{A}` (in-plane components)
      * - ``H``
        - node and element
        - :math:`\nu\,\nabla\times\mathbf{A}` (the magnetization M is not subtracted)
      * - ``Jt``
        - element (component 2)
        - STRANDED_CONDUCTOR: :math:`J_s`; SOLID: :math:`-\sigma\partial_tA`; CIRCUIT_COUPLED:
          :math:`\sigma(V/L-\partial_tA)`; STRANDED_COIL: :math:`n\,T_z(\mathrm{VOLT}-\mathrm{EMF})`; 0 elsewhere.

Boundary conditions and setup
   * ``_A`` component 2 on the external boundary (A = 0: flux parallel to the boundary; on the axis of an
     axisymmetric model A_theta = 0 by symmetry). A free boundary (no condition) is a flux-normal (symmetry) plane.
   * A solid conductor needs nothing: it is a closed conductor carrying the induced current.
   * Every circuit-coupled conductor needs a coupling of its V over its nodes (:doc:`couplings`), and then either a
     current (``AMPS`` on the master node), a voltage (``_V``) or ``_V = 0`` (shorted, as a solid conductor).
   * Every coil needs a coupling of VOLT and EMF, and a current (``CURR``) or a voltage (``_VOLT``).
   * Examples: ``examples/coil_2D_avolt*`` (circuit-coupled conductor and a solid ring), ``examples/coil_2D_avoltemf*``
     (stranded coil and a solid ring), ``pyTEC3D_APDL_benchmarks/suites/simple/cases/two_wires_2D``, ``examples/ITER_2D*`` (solid vessel structures).

.. _el-av3d-nodal:

``AVA3D_nodal`` (index 1)
-------------------------

3D formulation with a nodal (Lagrange) vector potential, three components per node, and a Coulomb-gauge penalty term
that makes the curl-curl operator definite. Theory: :ref:`theory-3d` and :ref:`theory-gauge`.

Fields
   .. list-table::
      :header-rows: 1
      :widths: 15 22 28 35

      * - Field
        - Finite elements
        - Domains
        - Notes
      * - ``A``
        - Lagrange vector, order ``option2``
        - all
        - Dirichlet values: ``_A`` components 0, 1, 2 (in the nodal frame, ``node_angles``).
      * - ``V``
        - Lagrange
        - SOLID_CONDUCTOR
        - Electric potential. Nodal load ``AMPS``.
      * - ``VOLT``
        - Lagrange
        - STRANDED_COIL
        - Coil terminal voltage. Nodal load ``CURR``.
      * - ``EMF``
        - Lagrange
        - STRANDED_COIL
        - Induced voltage.

Parameters
   As ``AVA2D_nodal``, with these differences: no ``thickness``; ``mn.AMPS`` is the nodal load of V in the solid
   conductors (current injected at the node [A]); the coil volume is the real constant ``rc.volume`` (volume of the
   whole coil) instead of ``radius``. ``Js``, ``M`` and ``T`` are 3-vectors.

Derived coefficients
   :math:`\nu`, :math:`\sigma`, :math:`n=N/(S\,R)`, :math:`k=1/(R\,\mathrm{volume})`, and the gauge coefficient
   :math:`\nu_g=\sqrt[3]{|\det\nu|}` (a scalar reluctivity of the cell).

Weak form
   .. list-table::
      :header-rows: 1
      :widths: 30 10 60

      * - Domain
        - Matrix
        - Term
      * - all
        - K
        - :math:`\int\nabla\times\mathbf{w}\cdot\nu\nabla\times\mathbf{A}+\int\nu_g\,(\nabla\cdot\mathbf{w})(\nabla\cdot\mathbf{A})`
      * - STRANDED_CONDUCTOR
        - F
        - :math:`\int\mathbf{w}\cdot\mathbf{J}_s`
      * - all but NON_CONDUCTING, SOFT_MAGNETIC
        - F
        - :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{M}`
      * - SOLID
        - C
        - :math:`\int\mathbf{w}\cdot\sigma\partial_t\mathbf{A}`, :math:`\int\nabla q\cdot\sigma\partial_t\mathbf{A}`
      * - SOLID
        - K
        - :math:`\int\mathbf{w}\cdot\sigma\nabla V`, :math:`\int\nabla q\cdot\sigma\nabla V`
      * - STRANDED_COIL
        - K, C
        - as in ``AVA2D_nodal``, with 3D vectors.

   The V row is the charge conservation :math:`\nabla\cdot\mathbf{J}=0` in weak form; ``AMPS`` on a node is a
   current source at that node (:ref:`theory-circuit-conductor`).

Outputs
   ``B`` and ``H`` (nodes and elements, 3 components); ``Jt`` (elements): STRANDED_CONDUCTOR :math:`\mathbf{J}_s`;
   SOLID :math:`-\sigma(\partial_t\mathbf{A}+\nabla V)`; STRANDED_COIL :math:`n\,\mathbf{T}(\mathrm{VOLT}-\mathrm{EMF})`.

Boundary conditions and setup
   * Far field: A = 0 (all components). Symmetry plane where B is tangent (current normal to the plane): the
     tangential components of A are 0. Symmetry plane where B is normal (current parallel to the plane): no condition
     on A (natural).
   * Axis of a sector model: the components of A normal to the axis are 0.
   * Planes that are not aligned with the global axes: rotate the nodal frames with ``node_angles`` (before the model
     is loaded) and fix the local components (:doc:`boundary_conditions`).
   * The potential V of each conductor is defined up to a constant: fix it on one node of every conductor that has
     no other condition on V (``_V = 0``).
   * Cyclic symmetry: :meth:`~pytec3d.model.Model.periodic_boundaries` with ``("A", "V")`` and rotated frames on the
     second face (:doc:`couplings`).
   * Examples: ``examples/coil_3D_avoltemf*``, ``pyTEC3D_APDL_benchmarks/suites/simple/cases/sphere_coil_*``,
     ``examples/ITER_3D_evol_coils_edge``.

.. _el-av3d-edge:

``AVA3D_edge`` (index 2)
------------------------

``AVA3D_nodal`` with the vector potential discretized by first-order Nedelec (edge) elements of the first kind
[Nedelec1980]_, ``NEDELEC_ONE``: one DOF per edge, the circulation of A along it; A is tangentially continuous only,
which represents
fields at material interfaces correctly. Theory: :ref:`theory-3d`.

Differences from ``AVA3D_nodal``
   * The mesh must be **quadratic** (TET10, HEX20/27, PRISM15/18, PYRAMID13): the edge DOFs are stored on the midside
     nodes. ``option2`` is the order of the Lagrange fields V, VOLT and EMF (2 as MAPDL SOLID236, whose V is on all 20
     nodes); A is first-order Nedelec whatever ``option2`` (there is no second-order Nedelec element): refine the mesh
     near the conductor surfaces.
   * No Coulomb gauge term: the curl-curl matrix is singular in the non-conducting regions (gradient fields). Gauge
     it with a tree-cotree decomposition, :func:`pytec3d.tree_gauge.tree_gauge` (:ref:`theory-gauge`).
   * Nodal load of A: ``CSG`` (:doc:`loads`), on the midside points: the load of the edge DOFs (``AMPS`` and
     ``CURR`` as for the nodal type).
   * Dirichlet values of A: component 2 of ``_A`` on the midside nodes, the tangential circulation along the edge
     (``_A = 0`` on the midside nodes of a face fixes :math:`\mathbf{n}\times\mathbf{A}=0` there). No nodal frame
     rotation is needed (the values are scalars). A non-zero value is :math:`\int\mathbf{A}\cdot d\mathbf{l}` from the
     lexicographically smaller end point of the edge to the larger one (compare x, then y, then z: libMesh's edge
     orientation), e.g. a uniform field :math:`\mathbf{B}_0` through :math:`\mathbf{A}=\tfrac12\mathbf{B}_0\times
     \mathbf{x}` (the same holds for every Nedelec field). :func:`pytec3d.tree_gauge.edges` lists the edges of
     the mesh (or of some cells) with their midside point and their end points in this orientation, and
     :func:`pytec3d.tree_gauge.edge_circulations` integrates a known field along them.
   * Cyclic couplings pair the midside nodes as well; the edge orientation is matched automatically.

Examples: ``examples/coil_3D_avoltemf_edge*``,
``examples/ITER_3D_evol_coils_edge``.

.. _el-timint:

Time-integrated variants (``_timint``, indices 3, 4, 5)
-------------------------------------------------------

``AVA2D_nodal_timint``, ``AVA3D_nodal_timint`` and ``AVA3D_edge_timint`` inherit everything from their base type, except
that **every electric unknown is a time integral**:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Field
     - Base type
     - ``_timint`` type
   * - V (2D circuit-coupled, 3D solid)
     - voltage / potential [V]
     - :math:`\int_0^tV\,dt` [V s]
   * - coil VOLT
     - terminal voltage [V]
     - :math:`\int_0^t\mathrm{VOLT}\,dt` [V s]
   * - coil EMF
     - induced voltage [V]
     - :math:`\int_0^t\mathrm{EMF}\,dt=\Psi` (linked flux) [Wb]

All the coupling terms then multiply a time derivative and move to **C**, which makes the coupled block symmetric
(:ref:`theory-timint`). Consequences:

* the meaning of ``_V`` / ``_VOLT`` / ``_EMF`` values changes (time integrals); ``AMPS`` and ``CURR`` keep their
  meaning (currents);
* the electric DOFs have no row in K: on a **static** step they are fixed to zero, so no conduction or coil current
  flows on a static step;
* ``Jt`` uses the time derivatives: e.g. for a coil :math:`n\mathbf{T}(\partial_t\mathrm{VOLT}-\partial_t\mathrm{EMF})`.

Examples: ``examples/coil_2D_avolt_timint``, ``examples/coil_3D_avoltemf_timint``,
``examples/ITER_3D_evol_fine2``.

.. _el-tphi3d:

``TPHI3D_edge`` (index 6): T-phi formulation
--------------------------------------------

The magnetic field :math:`\mathbf{H}=\mathbf{T}-\nabla\phi`, with :math:`\phi` (``PHI``) a Lagrange field on the whole
mesh and :math:`\mathbf{T}` (``T``) a Nedelec field of the conductors and the cut layers only, so that the unknowns are
few where the field diffuses (conductors, superconductors) and one per node elsewhere. The classic T-Omega form of the
literature (also called H-phi) [Carpenter1977]_, [Biro1999]_, [Kuczmann2009]_: Faraday's law in the conductors,
:math:`\nabla\cdot\mathbf{B}=0` everywhere at every step, so a static step solves magnetostatics (with
:math:`\mathbf{T}=0`) and a transient step eddy currents. There is **no source current**: the field is driven by the
boundary values of :math:`\phi` (an applied field), by magnets (a constant ``M``), by the current or the voltage of a
conductor through its cut, or by a coil given as **T fixed on its edges** (Setup). Domains: NON_CONDUCTING,
SOFT_MAGNETIC, HARD_MAGNETIC (``M``), SOLID_CONDUCTOR, CUT_LAYER (T without conduction). Nonlinear laws:
``power_law`` (:math:`\rho(J)`), ``bh_curve_h`` (:math:`\mu(H)`; in a conductor the rate of B by the time scheme).
Theory: :ref:`theory-potentials`.

Weak form
   .. list-table::
      :header-rows: 1
      :widths: 30 10 60

      * - Domain
        - Matrix
        - Term
      * - SOLID_CONDUCTOR
        - K
        - :math:`\int\nabla\times\mathbf{t}\cdot\rho\nabla\times\mathbf{T}`
      * - SOLID_CONDUCTOR, CUT_LAYER
        - C
        - :math:`\int\mathbf{t}\cdot\mu\partial_t(\mathbf{T}-\nabla\phi)`
      * - all
        - K
        - :math:`\int\nabla q\cdot\mu(\nabla\phi-\mathbf{T})`, :math:`\mu=\mu_0\,\mathrm{mur}`
      * - HARD_MAGNETIC
        - F
        - :math:`\int\nabla q\cdot\mathbf{B}_r`, :math:`\mathbf{B}_r=\mu_0\mu_r\mathbf{M}` (not the law's :math:`\mu`)

   The matrix is not symmetric (the coupling of T and :math:`\phi` is in C one way, in K the other).

Outputs
   ``H`` :math:`=\mathbf{T}-\nabla\phi`, ``B`` :math:`=\mu\mathbf{H}` (+ :math:`\mathbf{B}_r`) (nodes and elements);
   ``Jt`` :math:`=\nabla\times\mathbf{T}` (elements, SOLID_CONDUCTOR: the eddy currents, and the current of a fixed-T
   coil in its winding); ``E`` :math:`=\rho\nabla\times\mathbf{T}` (elements, SOLID_CONDUCTOR; with a power law
   :math:`E(J)`): the loss density :math:`\mathbf{E}\cdot\mathbf{J}` and the critical state (:math:`|E|/E_c`).

Setup
   * ``_T`` = 0 (component 2) on the midside nodes of the conductor surfaces, and ``_PHI`` = 0 on the conductor nodes
     that are not on a surface: :math:`\phi` only extends the outer potential into the conductors.
   * Or a tree gauge of T: ``_PHI`` free everywhere, ``_T`` = 0 on the surface and on the tree of the conductor edges
     rooted on the surface points, ``tree_gauge.tree_gauge(grid, conductor_points, surface_points, array_name="_T")``.
     The same solution and the same number of unknowns (each removes one unknown per inner node, a :math:`\phi` or a
     tree edge of T). Benchmark: ``eddy_sphere`` (``TPHI3D_edge_tree``).
   * First order only (``"second_order": false`` in the library, option2 = 2 refused before the solve): T is a
     first-order Nedelec field, and the gradients of first-order Lagrange functions lie in its edge space (the Whitney
     pair). The accuracy is first order.
   * ``_PHI`` fixes the tangential field: :math:`\phi=-\mathbf{H}_0(t)\cdot\mathbf{x}` on the outer boundary applies the
     uniform field :math:`\mathbf{H}_0`; a constant :math:`\phi` is a flux-normal plane. A free boundary is a flux wall
     (:math:`\mathbf{n}\cdot\mathbf{B}=0`). Fix :math:`\phi` on at least one node.
   * Static steps (the model's job, the solver does not do it): at rest :math:`\mathbf{H}=-\nabla\phi` in the
     conductors too, so set ``_T`` = 0 on every conductor edge (:math:`\rho\nabla\times\nabla\times` alone leaves the
     gradient part of T free, a singular system) and free ``_PHI`` (NaN) at the inner conductor nodes. A permanent
     magnet is a constant ``M``: it acts at every step, static or transient, so a transient may start from the
     magnetostatic state of a static step: before the transient steps put the transient gauges back (``_T`` NaN but 0
     on the surface, ``_PHI`` = 0 at the inner nodes). Benchmark: ``eddy_sphere_static``. In a static step a T DOF
     without stiffness (a CUT_LAYER edge) keeps a value given by ``_T`` (a fixed coil), else it is set to zero.
     A static step has no currents in the conductors: after a transient it wipes the induced currents, so trapped
     fields, persistent currents and their relaxation need transient steps only (a long hold: a transient with a
     long time step).
   * **No cancellation error** in high-permeability iron: there is no source field to cancel.
   * A conductor with a **net current** (a ring with an induced current, a wire or a tape with a transport current, a
     winding modelled as a conductor) needs a **cut**: :ref:`el-cuts`.
   * Without a cut: fill the hole with a weak conductor (SOLID_CONDUCTOR, resistivity :math:`10^4`-:math:`10^6` times
     the ring's) [Kuczmann2009]_: the conductor is then simply connected and T carries the net current across the hole.
     In ``ring`` a :math:`10^6` fill gives the current of the cut to :math:`10^{-5}` (:math:`10^2`: 3 %; :math:`10^8`:
     round-off). Simpler to mesh than a cut; the hole cells add T unknowns and cannot carry an imposed transport
     current or voltage.
   * Transport current or voltage through the cut, several conductors, conductors in series or in parallel:
     :ref:`el-cuts`.
   * **A coil as T fixed on its edges** (a stranded coil with current drive: uniform current, no eddy current in the
     winding, no circuit). Take any field :math:`\mathbf{T}_0` with :math:`\nabla\times\mathbf{T}_0=\mathbf{J}` (the
     current density of the winding at unit current) that is **zero outside the winding and the hole it encloses**,
     and impose its circulations: ``_T`` = :math:`I(t)\int_e\mathbf{T}_0\cdot d\mathbf{l}` on every edge of those cells
     (component 2 of the midside points, in libMesh's orientation: :func:`pytec3d.tree_gauge.edges` and
     :func:`pytec3d.tree_gauge.edge_circulations`), at each step. The winding cells are
     SOLID_CONDUCTOR (their ``rho`` is not used: T is fixed; ``Jt`` shows the coil current), the hole cells CUT_LAYER
     (T without conduction, with their ``mur``: air or iron). :math:`\phi` does the rest. Several coils add their
     :math:`\mathbf{T}_0`, each with its own current, also where they overlap. There is no voltage drive: the induced
     voltage is computed afterwards, from the flux :math:`\Psi=\int\mathbf{T}_0\cdot\mathbf{B}` per ampere.

     *The hole.* A winding with turns around an axis needs :math:`\mathbf{T}_0\ne0` across the area each turn
     encloses: its circulation around the winding's cross-section is the current. For an axisymmetric winding of
     rectangular cross-section :math:`r_1<r<r_2`, :math:`z_1<z<z_2` with the current density :math:`J`
     (:math:`J=NI/((r_2-r_1)(z_2-z_1))`):

     .. math::

        \mathbf{T}_0=T_z\,\mathbf{e}_z,\qquad
        T_z=\begin{cases}J\,(r_2-r) & r_1<r<r_2\ \text{(the winding)}\\
        J\,(r_2-r_1) & r<r_1\ \text{(the hole)}\\ 0 & \text{elsewhere}\end{cases}
        \qquad z_1<z<z_2,

     zero for :math:`z` outside :math:`(z_1,z_2)`. "The hole" is the cylinder :math:`r<r_1`, :math:`z_1<z<z_2`: the
     part of the region the winding encloses over the height of the winding, a disk as thick as the coil (for a thin
     turn, a thin disk). There :math:`\mathbf{T}_0` is uniform, curl-free: no current. :math:`\mathbf{n}\times
     \mathbf{T}_0=0` on the whole boundary of winding + hole (:math:`T_z\to0` at :math:`r_2`, and
     :math:`\mathbf{T}_0` is normal to the planes :math:`z_1`, :math:`z_2`), so H stays tangentially continuous. A
     square winding uses :math:`\max(|x|,|y|)` in place of :math:`r` (exact on a Cartesian grid).

     *A conductor in the hole* (a core, a vessel) has its own T, which cannot be fixed. Move :math:`\mathbf{T}_0` out
     of it by a gradient: :math:`\mathbf{T}_0'=\mathbf{T}_0-\nabla g` with :math:`g=T_h\,z\,\chi`, :math:`T_h` the
     uniform value in the hole and :math:`\chi` = 1 on the nodes of the conductor's cells, 0 elsewhere. On the edges
     :math:`\nabla g` is exact (:math:`g(\text{end})-g(\text{start})`), so the edge values of :math:`\mathbf{T}_0'` are
     zero on every edge of the conductor (its surface included: ``_T`` = 0 there as the gauge wants) and carry the
     same current; the support becomes the winding and the rest of the hole. The conductor keeps its usual gauge.

     Benchmarks: ``coil_core`` (a square solenoid around a conducting core block, ``pytec3d_driver.solenoid_t0``:
     against A, V-A within 3 % (B) and 6 % (core current)), ``closed_core`` (``TPHI3D_edge_coil``: a coil on an iron
     leg, the leg its hole, CUT_LAYER iron cells: the same solution as ``PHI3D_nodal`` with Hs given on the edges,
     to 10⁻¹⁴, linear and with a B-H curve).

     The alternative, the winding as a SOLID_CONDUCTOR of high resistivity driven through a cut, suppresses its eddy
     currents (:math:`\mu\sigma d^2\ll\Delta t`; not beyond, the conditioning degrades) but gives the current the
     distribution of conduction (about :math:`1/r` across a round winding, crowding at corners), not the uniform one.
   * Benchmarks of the ``formulations`` suite: ``phi3d_sphere``, ``phi3d_magnet``, ``phi_bh_sphere``, ``magnet_bh``
     (static, no conductor); ``eddy_sphere``, ``eddy_sphere_static``, ``ring``, ``ring_transport``, ``powerlaw_slab``,
     ``powerlaw_hts``, ``coil_core``, ``closed_core``.
   * Examples: ``examples/ring_3D_tphi`` (a cut, induced current), ``examples/coil_core_3D_tphi`` (a fixed-T coil
     around a conducting core); tutorial :doc:`../tutorials/tutorial_tphi_cuts` (the three drives of a cut).

.. _el-cuts:

Cuts: conductors with a net current (T-phi)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Why.* :math:`\phi` is one single-valued nodal field, so outside the conductors, where :math:`\mathbf{H}=-\nabla\phi`,
the circulation of H around any closed loop is zero. A loop around a conductor carrying a net current I must give I
(Ampère): the region outside the conductors must be made **simply connected** by **cuts**, surfaces :math:`\Sigma` on
which T carries the jump, :math:`\mathbf{T}=I\nabla c` with c going from 0 to 1 across a layer of cells. Every loop
around the conductor then crosses the cut once and picks up I from T. One cut per independent loop (per hole of the
non-conducting region): a ring needs one disk spanning its hole; a wire or a tape carrying a transport current needs
a cut from its surface to the outer boundary (or, in 3D, between its terminals on the boundary); a 2D cross-section of
a conductor carrying :math:`I_z`, a strip from the conductor to the boundary.

*Step by step* (the same in 3D, ``TPHI3D_edge``, and in 2D, ``TPHI2D_edge``):

1. **Choose** :math:`\Sigma`: a set of mesh faces (3D) or mesh edges (2D) in the non-conducting region, from the
   conductor surface to the outer boundary or across the hole, not crossing other conductors or other cuts.
2. **The layer**: the cells with a face (an edge in 2D) on :math:`\Sigma`, on one side of it, one cell thick, conforming.
   Give them a TYPE row with ``option1`` = ``CUT_LAYER`` (T without conduction) and their real material (``mur``: air
   or iron).
3. **Classify the edges of the layer**. *Planar* edges lie on :math:`\Sigma` (the side c = 0) or on the opposite face
   of the layer (c = 1); *crossing* edges have one end on each side. ``_T`` = 0 (component 2 of their midside nodes) on
   the planar edges. The crossing edges carry the current: the circulation of :math:`I\nabla c` along each one is
   :math:`\pm I`, + when its libMesh orientation (from the lexicographically smaller end point, x then y then z, to
   the larger) goes from the :math:`\Sigma` side to the other side, - otherwise (:func:`pytec3d.tree_gauge.edges` of
   the layer cells gives each edge with its start and end point in that orientation).
4. **The junction with the conductor**: the conductor has ``_T`` = 0 on its surface edges, except the crossing edges
   of the layer that lie on its surface: those belong to the cut (its value), not to the ``_T`` = 0 set.
5. **The drive**, one of:

   * *induced* (the current is an unknown, e.g. a ring in a field): couple all the crossing edges on one DOF, ``CP``
     id 1 on their midside nodes and ids 2 and 3 on the corner nodes of the two sides (they give the solver the
     orientation of each edge), ``add_dof_couplings({1: ["T"], 2: ["T"], 3: ["T"]})``; that DOF is the net
     current (benchmark ``ring``);
   * *imposed current*: ``_T`` = :math:`\pm I(t)` on every crossing edge (signed as in 3), no coupling (``ring_transport``
     ``current``, ``norris_wire``);
   * *imposed voltage* around the loop: the coupling of *induced* and ``VLTG`` = V(t) on one midside node of a crossing
     edge (``ring_transport`` ``voltage``).
6. **The outer boundary.** A free boundary (flux wall) needs nothing. Where the cut meets a boundary with a given
   :math:`\phi`, :math:`\phi` jumps by I across the cut, but it is one nodal field: the jump is across the layer
   (carried there by T), so the boundary nodes on :math:`\Sigma` take the value of the side away from the layer and the
   nodes on the far face of the layer that of the layer side (``norris_wire``: :math:`\phi=-I\theta/2\pi`, the nodes
   on the cut ray at :math:`\theta=2\pi`, those above it near 0).

*Several conductors.* Each conductor with a net current has its own cut, with its own coupling ids (1, 2, 3 for the
first, 4, 5, 6 for the second, ...); the cuts must not cross each other. Cuts are exactly what makes the air simply
connected, so making more of them is the way, with one rule: **one free cut DOF per independent hole**. A second free
cut on the same loop is singular (a constant shift of :math:`\phi` between the two cuts trades one DOF for the other);
an extra cut whose crossing edges are all given (0, or a current) is harmless. Connections:

* *series* (the same current): one coupling id for the crossing edges of all their cuts, or the same imposed I on each;
* *parallel* (the same voltage, the currents sharing by themselves): each its own coupled DOF, each loaded with the
  same ``VLTG``;
* a *total* current shared between conductors in parallel, a resistor in series, an external network: not available
  (a circuit of cut DOFs would be needed; the A, V-A types have coils with circuits).

Without a cut a ring can have its hole filled with a weak conductor (below). A coil given as T fixed on its edges is
a cut of the same kind, thick, with its current spread over the winding (:ref:`el-tphi3d`).

Examples: ``examples/ring_3D_tphi`` (3D, induced), ``examples/wire_2D_tphi`` (2D, imposed current, the jump of
:math:`\phi` on the box); :doc:`../tutorials/tutorial_tphi_cuts` runs the induced, current and voltage drives and a
ring without cut.

.. _el-tphi2d:

``TPHI2D_edge`` (index 8): 2D T-phi
-----------------------------------

The T-phi formulation on 2D cells: :math:`\mathbf{H}=\mathbf{T}-\nabla\phi` in the plane, T an in-plane Nedelec field
(NEDELEC_ONE on QUAD8 / TRI6 cells: quadratic cells, the DOFs on the midside nodes), the current density and the
electric field out of the plane (component 2 of ``Jt`` and ``E``): the 2D H-formulation of superconductor
cross-sections (tapes, wires, cables, bulks in a transverse field). The same weak form, domains, laws (``power_law``,
``power_law_b``, ``bh_curve_h``), gauges, static steps, cuts and fixed-T coils as ``TPHI3D_edge``.

* **Planar**: an infinite length along z, the integrals per ``thickness`` (``rc.thickness``, default 1 m: results per
  metre). A conductor carrying a net :math:`I_z` (a tape or a wire with a transport current) needs a cut: a one-cell
  strip of ``CUT_LAYER`` cells from its surface to the boundary (:ref:`el-cuts`). A conductor in a transverse field
  with no net current needs none.
* **Axisymmetric** (``option3`` = 1): (x, y) = (r, z), the integrals times :math:`2\pi r`, component 2 the
  azimuthal one: :math:`J_\theta=(\nabla\times\mathbf{T})_\theta=\partial_zT_r-\partial_rT_z`. The axis is a natural
  boundary of :math:`\phi` (:math:`B_r=0`). A conductor whose section touches the axis (a ball, a cylinder) needs no
  cut: the air of the half-plane is simply connected. A ring, a section away from the axis, with a net azimuthal
  current (a coil turn, a vessel) is encircled by air in the half-plane: it needs a cut, a strip from its section to
  the axis or to the boundary.
* Benchmarks: ``phi2d_cylinder``, ``phi2d_axisym`` (magnetostatics, 0.3 % and 0.9 % from the exact H),
  ``eddy_sphere`` (``TPHI2D_edge_axisym``: the conducting ball as a half model, its moment within 2.3 %),
  ``norris_wire`` (the AC transport loss of a round wire against Norris's formula).
* Examples: ``examples/wire_2D_tphi`` (planar, a cut, power law), ``examples/bulk_2D_tphi_axisym`` (axisymmetric,
  ``power_law_b``, a static step then transient steps); tutorial :doc:`../tutorials/tutorial_superconductors`.

.. _el-phi3d:

``PHI3D_nodal`` (index 7): magnetostatics with the scalar potential
-------------------------------------------------------------------

:math:`\mathbf{H}=\mathbf{H}_s-\nabla\phi`, with :math:`\phi` (``PHI``, Lagrange) and the source field
:math:`\mathbf{H}_s` (``Hs``, Nedelec) both unknowns of **one** static system: :math:`\mathbf{H}_s` has the curl of the
source current, :math:`\phi` makes B solenoidal. Domains: NON_CONDUCTING, SOFT_MAGNETIC (``bh_curve_h``),
HARD_MAGNETIC (``M``), STRANDED_CONDUCTOR (``me.Js``, the source current density, a cell array). Static
(``time_order = 0``), first order only (as T-phi). Theory: :ref:`theory-potentials`.

Weak form
   .. list-table::
      :header-rows: 1
      :widths: 30 10 60

      * - Domain
        - Matrix
        - Term
      * - all
        - K
        - :math:`\int\nabla\times\mathbf{h}\cdot\nabla\times\mathbf{H}_s` (the rows of Hs)
      * - STRANDED_CONDUCTOR
        - F
        - :math:`\int\nabla\times\mathbf{h}\cdot\mathbf{J}_s`: :math:`\nabla\times\mathbf{H}_s` is the solenoidal part
          of :math:`\mathbf{J}_s` (any ``Js`` is consistent)
      * - all
        - K
        - :math:`\int\nabla q\cdot\mu(\nabla\phi-\mathbf{H}_s)` (the rows of :math:`\phi`)
      * - HARD_MAGNETIC
        - F
        - :math:`\int\nabla q\cdot\mathbf{B}_r`, :math:`\mathbf{B}_r=\mu_0\mu_r\mathbf{M}`

   The coupling is one way (the Hs rows do not see :math:`\phi`): the matrix is block triangular, not symmetric. With
   the B-H law Newton refactors it all at each iteration, the Hs block included.

Outputs
   ``H`` :math:`=\mathbf{H}_s-\nabla\phi`, ``B`` :math:`=\mu\mathbf{H}` (+ :math:`\mathbf{B}_r`) (nodes and elements);
   ``Jt`` :math:`=\nabla\times\mathbf{H}_s` (elements: the current the source field carries); ``Hs`` (nodes).

Setup
   * Gauge of Hs: a tree of all the edges rooted on the box, ``tree_gauge.tree_gauge(grid, all_points, box_points,
     array_name="_Hs")``, and ``_Hs`` = 0 on the box edges (:math:`\mathbf{n}\times\mathbf{H}_s=0`). Needed also
     without a source (Hs = 0).
   * ``_PHI`` as in T-phi: the applied field on the boundary, or a free boundary (a flux wall) with one node fixed.
   * **No cancellation error** in iron: :math:`\mathbf{H}_s` and :math:`\nabla\phi` are in the same edge space (the
     gradients of first-order Lagrange functions are Nedelec functions), so where they cancel (in high-:math:`\mu`
     iron, :math:`\mathbf{H}\approx0`) they cancel exactly. A source field per cell (piecewise constant) is not in that
     space: it is the textbook cancellation case.
   * **Hs given instead of solved**: ``_Hs`` on every edge (the circulations of a known field with the right curl,
     e.g. of the Biot-Savart field of the coils, component 2 of the midside points in libMesh's orientation) removes
     the Hs rows: only :math:`\phi` is solved. Circulations, not point values: their discrete curl is then exactly
     the flux of the current through each face (Stokes), so the cancellation stays exact.
   * Cost: Hs adds one unknown per edge of the whole mesh (air included) to the system of :math:`\phi`; a separate
     source solve would need its own factorization anyway, and the edge values transferred.
   * Benchmarks: ``stranded_sphere``, ``closed_core`` (the iron frame around a coil: Js, Hs given, T0 as M; against
     A, V-A), ``phi3d_sphere``, ``phi3d_magnet``, ``phi_bh_sphere``, ``magnet_bh``.
   * Example: ``examples/core_3D_phi`` (a coil on a closed iron frame, B-H curve); tutorial
     :doc:`../tutorials/tutorial_scalar_potential` (Hs solved or given on the edges, linear and saturating iron).

.. _el-laws:

Nonlinear laws in the element types
-----------------------------------

The law ``bh_curve`` (source ``mat``: selected by MATERIAL id) replaces the coefficient ``nu`` (and ``nu_gauge`` of the
3D nodal types) by values computed at each quadrature point from :math:`|\mathbf{B}|`:

* secant :math:`\nu_s=H(|B|)/|B|` (times the identity) in the residual,
* tangent :math:`\nu_s\mathbf{I}+(H'(|B|)-\nu_s)\,\hat{\mathbf{b}}\otimes\hat{\mathbf{b}}` in the Jacobian of the
  curl-curl term; ``nu_gauge`` uses the secant only.

It applies to every element type and every domain (each has ``nu``). See :doc:`nonlinear` and
:doc:`../theory/nonlinear`.

How an element type is written
------------------------------

Each element type of ``element_types.json`` lists:

* ``dimension``, ``domains`` (the domains a TYPE row may give it), ``time_order``, ``out_of_plane_length`` (2D),
  ``constants`` (``mu0``);
* ``parameters``: ``mat.`` material properties, ``rc.`` real constants, ``me.`` element loads (cell arrays), ``mn.``
  nodal loads (point arrays), each with the domains where it exists and an optional default; the suffix ``[]`` or ``[][]`` makes it a vector or a tensor;
* ``expressions``: derived coefficients, C++ expressions on the libMesh types, optionally per domain;
* ``dof_fields``, ``nodal_fields``, ``element_fields``: the fields, their finite elements and domains; a Lagrange
  DOF field may name an ``mn.`` parameter of its own type as its ``nodal_loads`` (a scalar for a scalar field, a
  vector, in the nodal frames, for a vector field: :doc:`loads`); Lagrange fields have the order ``option2``,
  NEDELEC_ONE fields first order;
* ``second_order``: ``false`` when the type cannot run second order (``TPHI3D_edge``, ``PHI3D_nodal``: option2 = 2 is
  refused by the model check); absent by default;
* ``dof_kernels`` (``K`` and ``F`` lists), ``nodal_kernels``, ``element_kernels``: the terms, each
  ``{"trial": "grad(dt(V))", "coeff": "-sigma", "test": "A", "domain": [...]}``. A trial with ``dt()`` goes to C;
  a coefficient is a parameter or expression with an optional sign, or a sign alone (``"-"``: minus one);
* ``inherits``: start from another type and replace sections key by key (how the ``_edge`` and ``_timint`` types are
  defined).

The full grammar is in the docstring of ``solver/tools/generate_registries.py``, which checks the definitions (names,
types of every product, operators against the finite elements) and compiles them into C++. A new element type, domain
or law is added by editing the JSON files and rebuilding the solver module.

Writing expressions
-------------------

The ``constants``, the ``expressions`` and the C++ ``default`` values of an element type are C++ expressions. The
generator copies them into the solver, so they are written with the C++ operators and the libMesh vector and tensor
classes. Each name of a parameter or expression in the text is replaced by its value on the element. The key gives
the type of the result: ``nu`` a ``Real``, ``coil_j[]`` a ``RealVectorValue``, ``nu[][]`` a ``RealTensorValue``. An
expression that gives another type does not compile.

The types and their API
~~~~~~~~~~~~~~~~~~~~~~~

``Real`` is ``double``. ``RealVectorValue`` and ``RealTensorValue`` are libMesh's 3-component vector and 3 x 3 tensor
(``VectorValue<Real>`` and ``TensorValue<Real>``). Most of their operations come from the base classes ``TypeVector``
and ``TypeTensor``, documented in the libMesh Doxygen:

* `TypeVector <https://libmesh.github.io/doxygen/classlibMesh_1_1TypeVector.html>`_ and
  `VectorValue <https://libmesh.github.io/doxygen/classlibMesh_1_1VectorValue.html>`_;
* `TypeTensor <https://libmesh.github.io/doxygen/classlibMesh_1_1TypeTensor.html>`_ and
  `TensorValue <https://libmesh.github.io/doxygen/classlibMesh_1_1TensorValue.html>`_;
* the free functions (``outer_product``) are in the
  `libMesh namespace <https://libmesh.github.io/doxygen/namespacelibMesh.html>`_.

The operations most useful in a coefficient:

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Expression
     - Result
     - Meaning
   * - ``RealVectorValue(x, y, z)``
     - vector
     - a vector from its components
   * - ``RealTensorValue(xx, xy, xz, yx, yy, yz, zx, zy, zz)``
     - tensor
     - a tensor from its components, row by row
   * - ``v(i)``, ``T(i, j)``
     - Real
     - a component (0 = x, 1 = y, 2 = z)
   * - ``a * v``, ``a * T``, ``v + w``, ``T - S``, ``-v``
     - vector / tensor
     - scaling and sums
   * - ``v * w``
     - Real
     - dot product :math:`\mathbf{v}\cdot\mathbf{w}` (not component by component)
   * - ``v.cross(w)``
     - vector
     - cross product :math:`\mathbf{v}\times\mathbf{w}`
   * - ``v.norm()``, ``v.norm_sq()``, ``v.unit()``
     - Real / vector
     - length, squared length, unit vector
   * - ``outer_product(v, w)``
     - tensor
     - :math:`\mathbf{v}\otimes\mathbf{w}`, components :math:`v_i w_j`
   * - ``T * v``, ``v * T``
     - vector
     - :math:`\mathbf{T}\mathbf{v}` and :math:`\mathbf{v}^T\mathbf{T}`
   * - ``T * S``
     - tensor
     - matrix product
   * - ``T.transpose()``, ``T.inverse()``
     - tensor
     - transpose, inverse
   * - ``T.det()``, ``T.tr()``, ``T.norm()``
     - Real
     - determinant, trace, Frobenius norm
   * - ``T.contract(S)``
     - Real
     - double contraction :math:`\mathbf{T}:\mathbf{S} = \sum_{ij} T_{ij} S_{ij}`
   * - ``T.row(i)``
     - vector
     - row ``i``

Names available besides the parameters and the expressions:

* ``I``: the identity tensor; ``pi``; the constants of the element type (``mu0``);
* ``axisymmetric``: a ``bool``, true in an axisymmetric 2D model;
* ``abs``, ``sqrt``, ``cbrt``, ``pow``, ``exp``, ``log``, ``sin``, ``cos``, ``tan``, ``atan2``, ``hypot``, ``min``,
  ``max``: the ``std::`` functions of the same name, usable without the prefix;
* any other function of the C++ standard math library (``<cmath>``), written with its ``std::`` prefix (next
  section);
* the conditional ``c ? a : b``, with both branches of the same type.

Standard library functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The whole of ``<cmath>`` is available with the ``std::`` prefix, e.g. ``std::tanh``, ``std::sinh``, ``std::cosh``,
``std::asin``, ``std::acos``, ``std::atan``, ``std::erf``, ``std::erfc``, ``std::exp2``, ``std::expm1``,
``std::log10``, ``std::log2``, ``std::log1p``, ``std::floor``, ``std::ceil``, ``std::round``, ``std::trunc``,
``std::fmod``, ``std::copysign``, ``std::fmin``, ``std::fmax``, ``std::isfinite`` (see the
`cppreference list <https://en.cppreference.com/cpp/numeric/math>`_). The short names of the previous list may be
written with the prefix too (``std::sqrt`` is ``sqrt``). They take and give ``Real`` values; apply them to a
component or a norm of a vector or a tensor, not to the vector itself.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Expression
   * - ``rho_T``
     - ``rho0 * (1 + alpha * (T - T0))``: a resistivity linear in a temperature (scalar real constants)
   * - ``b_sat``
     - ``Bs * std::tanh(b / Bs)``: a smooth saturation at ``Bs``
   * - ``b_reg``
     - ``std::sqrt(B0.norm_sq() + h * h)``: a length that does not vanish
   * - ``angle``
     - ``std::atan2(T(1), T(0))``: the angle of a direction in the xy plane
   * - ``bounded``
     - ``std::fmin(std::fmax(x, lo), hi)``: ``x`` limited to ``[lo, hi]``
   * - ``finite``
     - ``std::isfinite(x) ? x : 0.0``

Only ``<cmath>`` and the libMesh types are included, and the solver is C++17. So ``std::clamp`` (``<algorithm>``)
does not compile; use ``std::fmin(std::fmax(x, lo), hi)``. ``std::numbers::pi`` (C++20) does not compile either;
use ``pi``. The name after ``std::`` is replaced like any other name: do not give a parameter or an expression the
name of a function you call (a parameter ``tanh`` breaks ``std::tanh``).

Rules
~~~~~

* **Order and domains**: an expression may only use the parameters and the expressions declared before it, and only
  those that exist in all of its domains (a ``domain`` list limits where it exists). The generator checks both.
* **Global axes**: vectors and tensors from the model (``mat.``, ``rc.``, ``me.``) are given in the element coordinate
  system of the cell (:doc:`coordinate_systems`). They are rotated to the global axes before any expression sees
  them, and the result of an expression is a global value. ``mur(0, 0)`` is therefore the global xx component.
  Write expressions so that they do not depend on the axes (``det``, ``tr``, ``norm``, products of the parameters),
  or take a direction from a parameter (``T`` below) instead of a fixed component.
* **Once per element**: an expression is evaluated once per element and kept. It is constant over the element, not
  evaluated at each quadrature point. A coefficient that varies in space varies through a cell array (``me.``).
  Coefficients that depend on the solution are written as nonlinear laws (:doc:`nonlinear`).
* **Point data** (``mn.``) cannot appear in an expression: it is only a nodal load.

Examples
~~~~~~~~

From the library (``element_types.json``):

.. code-block:: json

   {
       "constants": {"mu0": "4 * pi * 1e-7"},
       "expressions": {
           "nu[][]": "(mu0 * mur).inverse()",
           "nu_gauge": "cbrt(abs(nu.det()))",
           "sigma[][]": {"value": "rho.inverse()", "domain": ["SOLID_CONDUCTOR", "CIRCUIT_COUPLED_CONDUCTOR"]},
           "turns_density": {"value": "turns / (cross_section * resistance)", "domain": ["STRANDED_COIL"]},
           "coil_volume": {"value": "(axisymmetric ? 2 * pi * radius : 1.0) * cross_section * thickness",
                           "domain": ["STRANDED_COIL"]},
           "coil_j[]": {"value": "turns_density * T", "domain": ["STRANDED_COIL"]}
       }
   }

* ``nu`` is a tensor from a scalar constant and a tensor property: the reluctivity :math:`(\mu_0\mu_r)^{-1}`.
* ``nu_gauge`` makes a scalar out of a tensor, the cube root of its determinant (the geometric mean of the
  principal values). It does not depend on the axes.
* ``coil_volume`` uses ``axisymmetric``: the volume of the coil per unit length, or per radian times
  :math:`2\pi r`.
* ``coil_j`` scales a vector property, the winding direction ``T`` of the coil.

Combining parameters: a resistivity :math:`\rho_\parallel` along a direction ``T`` and :math:`\rho_\perp` across it,
from two scalar real constants and a vector one:

.. code-block:: json

   {
       "parameters": {"rc.rho_par": {}, "rc.rho_perp": {}, "rc.T[]": {}},
       "expressions": {
           "t[]": "T.unit()",
           "Pt[][]": "outer_product(t, t)",
           "rho[][]": "rho_par * Pt + rho_perp * (I - Pt)",
           "sigma[][]": "rho.inverse()"
       }
   }

``Pt`` projects a vector onto ``t`` and ``I - Pt`` onto the plane normal to it, so ``rho`` has the principal value
:math:`\rho_\parallel` along ``t`` and :math:`\rho_\perp` across it, in any orientation of ``T``.

More one-liners:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Expression
   * - ``Jm[]``
     - ``nu * M``: a tensor applied to a vector
   * - ``w``
     - ``M * B0`` (``B0`` a vector parameter): a dot product, a scalar
   * - ``n[]``
     - ``T.cross(B0).unit()``: the unit normal to two directions
   * - ``energy``
     - ``B0.norm_sq() / (2 * mu0)``
   * - ``mur_t[][]``
     - ``RealTensorValue(mur(0, 0), 0, 0, 0, mur(0, 0), 0, 0, 0, mur(2, 2))``: a tensor from chosen global
       components
   * - ``Q[][]``
     - ``RealTensorValue(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, 0, 0, 1)``: the rotation by
       ``alpha`` about z (``alpha`` a scalar parameter); ``Q * mur * Q.transpose()`` turns a tensor by it
   * - ``eps``
     - ``max(T.norm(), 1e-12)``: a guard before dividing by a length

C++ pitfalls
~~~~~~~~~~~~

* ``1 / 2`` is an integer division and gives 0. Write ``1.0 / 2`` or ``0.5``.
* ``^`` is not a power (it is a bitwise XOR and does not compile on ``Real``). Use ``pow(x, 2)`` or ``x * x``.
* ``v * w`` of two vectors is their dot product. There is no component-by-component product: build it with
  ``RealVectorValue(v(0) * w(0), v(1) * w(1), v(2) * w(2))``.
* A name that is neither a parameter nor one of the names above is passed to C++ unchanged, so a typo is reported by
  the C++ compiler, in the generated ``registries.h``, when the solver module is rebuilt.
* A default value of a ``mat.`` or ``rc.`` parameter may be a C++ expression of the constants, e.g.
  ``"default": "RealVectorValue(0, 0, 1)"``. It cannot use other parameters.
