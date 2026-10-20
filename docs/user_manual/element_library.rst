.. _element-library:

Element Library
===============

An **element type** is a complete finite-element formulation: its unknowns (DOF fields) and their finite elements,
the parameters it reads from the model data and the mesh, the terms of its weak form, and the fields it computes after
each solve. A **domain** (``option1``) selects, inside an element type, the physical region a cell represents: the
same element type gives air, a magnet, a coil or a massive conductor depending on the domain.

The library is defined in ``src/pytec3d/resources/element_types.json`` (with ``fields.json``, ``domains.json`` and
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
     - ``AV2D_nodal``
     - 2
     - A_z (Lagrange); V, VOLT, EMF (Lagrange, out of plane)
     - 2D planar or axisymmetric eddy currents, conductors with a terminal, stranded coils.
   * - 1
     - ``AV3D_nodal``
     - 3
     - A (Lagrange vector); V, VOLT, EMF (Lagrange)
     - 3D eddy currents, nodal vector potential with a Coulomb gauge (penalty) term.
   * - 2
     - ``AV3D_edge``
     - 3
     - A (Nedelec, first kind); V, VOLT, EMF (Lagrange)
     - 3D eddy currents with edge elements: tangential continuity only, tree-cotree gauge.
   * - 3
     - ``AV2D_nodal_timint``
     - 2
     - as 0, V, VOLT and EMF time-integrated
     - Symmetric coupling in 2D.
   * - 4
     - ``AV3D_nodal_timint``
     - 3
     - as 1, V, VOLT and EMF time-integrated
     - Symmetric coupling, nodal.
   * - 5
     - ``AV3D_edge_timint``
     - 3
     - as 2, V, VOLT and EMF time-integrated
     - Symmetric coupling, edge.

All the element types are first order in time (``time_order = 1``: the theta method, see
:ref:`theory-time-integration`). A model uses one element type for all its cells, with any of its domains.

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
       :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{M}`.
     - A. Material: ``mur``, ``M``.
   * - ``STRANDED_CONDUCTOR`` (3)
     - Imposed current density, no eddy currents: :math:`\mathbf{J}=\mathbf{J}_s`.
     - A. Material: ``mur`` (``M``). Element load ``Js``.
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
       ``resistance``, ``symmetry``, ``T``, and ``radius`` (2D) or ``volume`` (3D). Nodal load ``CURR``.

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

``AV2D_nodal`` (index 0)
------------------------

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
   * Examples: ``tests/apdl_benchmarks/coil_2D_avolt*`` (circuit-coupled conductor and a solid ring), ``tests/apdl_benchmarks/coil_2D_avoltemf*``
     (stranded coil and a solid ring), ``tests/apdl_benchmarks/two_wires_2D``, ``tests/examples/ITER_2D*`` (solid vessel structures).

.. _el-av3d-nodal:

``AV3D_nodal`` (index 1)
------------------------

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
   As ``AV2D_nodal``, with these differences: no ``thickness``; ``mn.AMPS`` is the nodal load of V in the solid
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
        - as in ``AV2D_nodal``, with 3D vectors.

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
   * Examples: ``tests/apdl_benchmarks/coil_3D_avoltemf*``, ``tests/examples/3D_loop_cond_base``, ``tests/apdl_benchmarks/sphere_coil_*``,
     ``tests/examples/ITER_3D_evol_coils``.

.. _el-av3d-edge:

``AV3D_edge`` (index 2)
-----------------------

``AV3D_nodal`` with the vector potential discretized by first-order Nedelec (edge) elements of the first kind,
``NEDELEC_ONE``: one DOF per edge, the circulation of A along it; A is tangentially continuous only, which represents
fields at material interfaces correctly. Theory: :ref:`theory-3d`.

Differences from ``AV3D_nodal``
   * The mesh must be **quadratic** (TET10, HEX20/27, PRISM15/18, PYRAMID13): the edge DOFs are stored on the midside
     nodes. ``option2`` is the order of the Lagrange fields V, VOLT and EMF.
   * No Coulomb gauge term: the curl-curl matrix is singular in the non-conducting regions (gradient fields). Gauge
     it with a tree-cotree decomposition, :func:`pytec3d.tree_gauge.tree_gauge` (:ref:`theory-gauge`).
   * Dirichlet values of A: component 2 of ``_A`` on the midside nodes, the tangential circulation along the edge
     (``_A = 0`` on the midside nodes of a face fixes :math:`\mathbf{n}\times\mathbf{A}=0` there). No nodal frame
     rotation is needed (the values are scalars).
   * Cyclic couplings pair the midside nodes as well; the edge orientation is matched automatically.

Examples: ``tests/apdl_benchmarks/coil_3D_avoltemf_edge*``, ``tests/examples/two_wires_3D_edge*``, ``tests/examples/CREATE_*``,
``tests/examples/ITER_3D_evol_*edge``.

.. _el-timint:

Time-integrated variants (``_timint``, indices 3, 4, 5)
-------------------------------------------------------

``AV2D_nodal_timint``, ``AV3D_nodal_timint`` and ``AV3D_edge_timint`` inherit everything from their base type, except
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

Examples: ``tests/apdl_benchmarks/coil_2D_avolt_timint``, ``tests/apdl_benchmarks/coil_2D_avoltemf_timint``, ``tests/apdl_benchmarks/coil_3D_avoltemf_timint``,
``tests/examples/ITER_3D_evol``.

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
  nodal loads (point arrays), each with the domains where it exists and an optional default; the suffix ``[]`` or
  ``[][]`` makes it a vector or a tensor;
* ``expressions``: derived coefficients, C++ expressions on the libMesh types, optionally per domain;
* ``dof_fields``, ``nodal_fields``, ``element_fields``: the fields, their finite elements and domains;
* ``dof_kernels`` (``K`` and ``F`` lists), ``nodal_kernels``, ``element_kernels``: the terms, each
  ``{"trial": "grad(dt(V))", "coeff": "-sigma", "test": "A", "domain": [...]}``. A trial with ``dt()`` goes to C;
* ``inherits``: start from another type and replace sections key by key (how the ``_edge`` and ``_timint`` types are
  defined).

The full grammar is in the docstring of ``src/tools/generate_registries.py``, which checks the definitions (names,
types of every product, operators against the finite elements) and compiles them into C++. A new element type, domain
or law is added by editing the JSON files and rebuilding the solver module.
