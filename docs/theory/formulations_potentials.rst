.. _theory-potentials:

Scalar potential formulations: PHI and T-phi
=============================================

Besides the A, V-A formulation (:doc:`formulation_2d`, :doc:`formulation_3d`), the element library has a
the magnetic scalar potential: for magnetostatics with source currents (``PHI3D_nodal``), and T-phi for eddy currents,
common in applied superconductivity (``TPHI3D_edge``). The names follow the usual convention of the literature: the comma groups the potentials of the
conductors and the dash separates them from those of the non-conducting regions. "A, V-A" is one formulation, with A and
V in the conductors and A outside them. The element type names spell them without the
punctuation (``AVA3D_edge`` is A, V-A). This page derives its weak form; the element types are listed in
:doc:`../user_manual/element_library`. The notation is that of :doc:`notation`: :math:`\nu = \mu^{-1}` the reluctivity,
:math:`\rho = \sigma^{-1}` the resistivity, and the material law :math:`\mathbf{B} = \mu(\mathbf{H} + \mathbf{M})`
(with a B-H curve :math:`\mathbf{B} = \mathbf{B}_\text{curve}(\mathbf{H}) + \mathbf{B}_r`, :math:`\mathbf{B}_r =
\mu_0\mu_r\mathbf{M}`: :doc:`constitutive`).

.. contents:: On this page
   :local:
   :depth: 1

Names of the formulations
-------------------------

A formulation is named by what it solves in the conductors and in the rest of the domain. The comma groups the potentials
of the conductors, the dash separates them from those outside [Kuczmann2009]_:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - In the conductors
     - Outside: A (vector potential)
     - Outside: :math:`\phi` (scalar potential)
   * - A and V
     - **A, V-A** (``AVA*``)
     - A, V-:math:`\phi` (not implemented)
   * - T (and :math:`\phi`), or H
     - T-A family: T-A, h-a, A, T-A (not implemented)
     - **T,** :math:`\phi` **-** :math:`\phi` = **T-phi** (``TPHI3D_edge``)

* **T-phi, T-Omega and H-phi are the same formulation**, with :math:`\mathbf{H}=\mathbf{T}-\nabla\phi` in the conductors and
  :math:`\mathbf{H}=-\nabla\phi` outside. The eddy-current literature calls it T-Omega or T, :math:`\phi`-:math:`\phi`,
  the superconductor literature H-phi, and the variants differ only in where the unknowns are supported. pyTEC3D
  calls it T-phi.
* **The T-A family replaces the scalar potential outside with A**, the conductors still being described by a current
  potential or by H: the T-A formulation is "basically an h(-:math:`\phi`)-a-formulation with a slab approximation"
  [Geuzaine2024]_.

Two element types use the scalar potential: ``PHI3D_nodal`` for magnetostatics with source currents, ``TPHI3D_edge``
for eddy currents in conductors and superconductors, without source currents.

The source field: ``PHI3D_nodal``
---------------------------------

**Field and potentials.** Where there is no current the field is a gradient. With source currents
:math:`\mathbf{J}_s` (stranded conductors), a source field :math:`\mathbf{H}_s` with :math:`\nabla\times\mathbf{H}_s =
\mathbf{J}_s` carries their curl and the scalar potential the rest [SimkinTrowbridge1979]_:

.. math::

   \mathbf{H} = \mathbf{H}_s - \nabla\phi, \qquad \mathbf{B} = \mu(\mathbf{H} + \mathbf{M}).

With :math:`\mathbf{H}_s = 0`, :math:`\phi` is the **total** potential; with :math:`\mathbf{H}_s \ne 0`, the
**reduced** potential. :math:`\mathbf{H}_s` is not the physical field of the currents (it need not vanish at infinity,
nor be their Biot-Savart field): only its curl matters, :math:`\phi` absorbs the rest.

**Equations.** Both are unknowns of one static system, :math:`\mathbf{H}_s` a Nedelec field and :math:`\phi` a Lagrange
field of the whole mesh, with the test functions :math:`\mathbf{h}` and :math:`q`:

.. math::

   \int\nabla\times\mathbf{h}\cdot\nabla\times\mathbf{H}_s &= \int\nabla\times\mathbf{h}\cdot\mathbf{J}_s, \\
   \int\nabla q\cdot\mu(\nabla\phi - \mathbf{H}_s) &= \int\nabla q\cdot\mathbf{B}_r - \oint q\,\mathbf{n}\cdot\mathbf{B}.

The first makes :math:`\nabla\times\mathbf{H}_s` the solenoidal part of :math:`\mathbf{J}_s` in the discrete space (any
:math:`\mathbf{J}_s`, e.g. cell values, is consistent); a tree gauge rooted on the boundary and
:math:`\mathbf{n}\times\mathbf{H}_s = 0` there make it regular. The coupling is one way: the matrix is block triangular,
not symmetric, and the Hs rows do not depend on :math:`\mu`. The boundary term is the natural condition: a free
boundary is a flux wall (:math:`\mathbf{n}\cdot\mathbf{B} = 0`); a fixed :math:`\phi` imposes the tangential field
(:math:`\phi = -\mathbf{H}_0\cdot\mathbf{x}`: a uniform applied field :math:`\mathbf{H}_0`).

* **No cancellation error.** The reduced potential is known for its cancellation error in iron
  (:math:`\mathbf{H}_s` and :math:`\nabla\phi` large and nearly opposite, [Mayergoyz1987]_). It comes from a source
  field outside the space of :math:`\nabla\phi` (e.g. Biot-Savart values at the nodes, or per cell). With
  :math:`\mathbf{H}_s` and :math:`\nabla\phi` in the same first-order edge space (the "compatible approximation" of
  [Kuczmann2009]_, §2.3.1 D) there is none, even in a closed iron loop around the coil, the case of MAPDL's
  generalized scalar potential: in the ``closed_core`` benchmark the scalar potential and A, V-A bracket the flux of
  the frame and the inductance of the coil from both sides, within 2-4 % on the standard mesh (1.4 % at half the cell
  size), independent of :math:`\mu_r`.
* **Hs given.** The circulations of a known field with the right curl on every edge (e.g. a Biot-Savart field of the
  coils) remove the Hs rows. Their discrete curl is the flux of the current through each face (Stokes), so the source
  stays exactly in the edge space: no cancellation error either (``closed_core``, ``PHI3D_nodal_fixed``).
* **Cost.** Hs adds one unknown per edge of the whole mesh, air included. Static, the solve is one factorization; a
  source computed in a separate session would need its own factorization of a system as large, and the transfer of its
  edge values.

T-phi: ``TPHI3D_edge``
----------------------

**Field and potentials.** In the conductors a current vector potential :math:`\mathbf{T}` carries the eddy currents,
:math:`\mathbf{J}=\nabla\times\mathbf{T}`:

.. math::

   \mathbf{H} = \mathbf{T} - \nabla\phi, \qquad \mathbf{B} = \mu(\mathbf{H} + \mathbf{M}),

T only in the conductors (and the cut layers), :math:`\phi` everywhere: few unknowns, the formulation of applied
superconductivity. There is no source current: the field comes from the boundary values of :math:`\phi`, from magnets
(a constant :math:`\mathbf{M}`), from currents and voltages imposed through the cuts of the conductors, and from coils
given as T fixed on their edges (below).

**Equations.** The classic T-Omega form of the eddy-current literature [Carpenter1977]_, [Biro1999]_: Faraday's law in the conductors
:math:`\Omega_c`, with :math:`\mathbf{E} = \rho\nabla\times\mathbf{T}`, and :math:`\nabla\cdot\mathbf{B} = 0`
everywhere. Galerkin, with the test functions :math:`\mathbf{t}` of T and :math:`q` of :math:`\phi`:

.. math::

   \int_{\Omega_c}\rho\nabla\times\mathbf{T}\cdot\nabla\times\mathbf{t} + \int_{\Omega_c}\partial_t\big(\mu(\mathbf{T}
   - \nabla\phi)\big)\cdot\mathbf{t} &= 0, \\
   \int\nabla q\cdot\mu(\nabla\phi - \mathbf{T}) &= \int\nabla q\cdot\mathbf{B}_r - \oint q\,\mathbf{n}\cdot\mathbf{B}.

**Rate form and divergence form.** The H-phi form of the superconductor literature [Bossavit1998]_,
[Dular2020]_, [Arsenault2021]_ tests Faraday's law with the gradients too: :math:`\nabla\cdot\partial_t\mathbf{B} = 0`
instead of :math:`\nabla\cdot\mathbf{B} = 0`. It is symmetric, but it keeps :math:`\nabla\cdot\mathbf{B}` at its initial
value: a constant magnetization never enters (it is the initial state), a static step has no equation for
:math:`\phi`, and a B-H law needs the rate of B even in non-conducting iron. The two forms have the same unknowns and the
same discrete space; with backward Euler they differ by the history term :math:`\int\nabla q\cdot\mathbf{B}^n`, which
the divergence form drops by imposing :math:`\nabla\cdot\mathbf{B}=0` at every step. The library uses the divergence
form:

* a **static step** solves magnetostatics (no eddy currents at rest: the model sets ``_T`` = 0 and frees
  :math:`\phi` inside the conductors, :math:`\mathbf{H}=-\nabla\phi` there too), and a transient may start from it;
* **permanent magnets** are a constant ``M``, a load on :math:`\phi` at every step;
* **B-H law** (``bh_curve_h``, :math:`\mu(H)` with the tangent :math:`dB/dH`): an ordinary secant and tangent in the
  :math:`\phi` equation; the rate of B, by the time scheme on B (backward Euler :math:`(B(H)-B(H^n))/\Delta t`), only
  where the iron conducts (any theta, :ref:`theory-nonlinear-elements`);
* the matrix is **not symmetric**: the coupling of T and :math:`\phi` is in C in the conductor equation and in K in
  the divergence equation (a direct solver does not mind; an iterative one needs GMRES).

**Discretization.** :math:`\phi` is a Lagrange field over the whole mesh and :math:`\mathbf{T}` (``T``) a Nedelec field
of the conductors with zero tangential trace on their surface (``_T`` = 0 on the surface edges). In the conductors
:math:`\phi` is fixed to zero at the inner nodes (``_PHI`` = 0), so that it only extends the outer potential into the
first layer of cells; :math:`\mathbf{T}` adds the rest. This decomposition is unique and spans the H(curl) space of the
conductors (the h-phi space of [Dular2020]_). T needs no gauge of its own: :math:`\mu\,\partial_t\mathbf{T}` is a
mass term; the only freedom is a gradient carried by T or by :math:`\nabla\phi`, removed either by :math:`\phi` = 0 at
the inner nodes or, equivalently, by a tree gauge of T in the conductors (T = 0 on a tree rooted on the surface,
:math:`\phi` free): the same solution and number of unknowns.

* :math:`\phi` is **first order** (``option2`` = 1, checked before the solve; ``PHI3D_nodal`` too): the gradients of
  first-order Lagrange functions lie in the first-order Nedelec space, which makes :math:`\mathbf{T} - \nabla\phi`
  (:math:`\mathbf{H}_s - \nabla\phi`) consistent (the Whitney complex). With second-order :math:`\phi` the uniform field
  is not even represented inside a conductor. There is no second-order Nedelec element: the formulation is first order.
* **A coil as T fixed on its edges.** A stranded coil with current drive needs no source field: any
  :math:`\mathbf{T}_0` with :math:`\nabla\times\mathbf{T}_0 = \mathbf{J}` (the current density at unit current) that
  vanishes outside the winding and the hole it encloses is a T of those cells, imposed on their edges as
  :math:`I(t)` times its circulations; :math:`\phi` does the rest. The winding has no eddy currents (its T is fixed)
  and the current is exactly uniform; several coils add their :math:`\mathbf{T}_0`, each with its current. In the hole
  :math:`\mathbf{T}_0` is curl-free (e.g. uniform, :math:`J(r_2-r_1)\,\mathbf{e}_z` inside a winding
  :math:`r_1<r<r_2`, :math:`z_1<z<z_2`, over its height). A conductor in the hole keeps its own T: a discrete gradient
  moves :math:`\mathbf{T}_0` out of it (:math:`\mathbf{T}_0-\nabla(T_h z\chi)`, :math:`\chi` = 1 on its nodes). It is
  the cut-layer construction of a ring current, thick, with the current distributed over the winding. The induced
  voltage, not an unknown, is :math:`d\Psi/dt` with :math:`\Psi = \int\mathbf{T}_0\cdot\mathbf{B}` per ampere
  (reciprocity, :math:`\mathbf{n}\times\mathbf{T}_0 = 0` on the boundary of its support). Setup and an example:
  :ref:`el-tphi3d`; benchmarks ``coil_core``, ``closed_core``.
* **Topology.** A conductor carrying a net current makes :math:`\phi` multivalued around it: cuts (thick or thin,
  [Arsenault2021b]_) or cohomology basis functions [Lahtinen2015]_ are needed. The library has the cut functions of
  [Geuzaine2024]_ through the ``CUT_LAYER`` domain and a DOF coupling (:ref:`el-tphi3d`). Without a cut, the hole can
  be filled with a weak conductor (:math:`\rho` 10⁴-10⁶ times the ring's), which makes the conductor simply connected
  ([Kuczmann2009]_, §2.3.3; benchmark ``ring``); A, V-A needs neither.
* **Nonlinear resistivity.** ``power_law`` replaces :math:`\rho` by :math:`E(J)/J` of a tabulated E-J curve (e.g. the
  power law :math:`E = \rho_0 J + E_c (J/J_c)^n` with a normal-state floor), with the tangent from :math:`dE/dJ`.
  ``power_law_b`` adds a critical current depending on the field, :math:`E=E_0(J/s(|B|))`, :math:`s=J_c(B)/J_c(0)`
  a second curve; its tangent is at fixed B (no cross-coupling block). The settings that make Newton reliable at
  high n (an absolute residual tolerance, a decreasing residual, a normal-state end of the E-J curve):
  :ref:`superconductor-settings`.
* **2D** (``TPHI2D_edge``). The same equations on 2D cells: T and H in the plane (in-plane Nedelec), J and E out of
  it, :math:`J_z=\partial_xT_y-\partial_yT_x` (axisymmetric: :math:`J_\theta=\partial_zT_r-\partial_rT_z`, the
  integrals times :math:`2\pi r`). It is the 2D H-formulation of superconductor cross-sections; a net transport
  current needs a cut, a strip from the conductor to the boundary.
* **High-permeability conductors** on a coarse mesh: the scalar and the vector potential bracket the field of an
  iron block from both sides (``coil_core``, :math:`\mu_r` = 100: the core field of the magnetostatics 24 % apart at
  5 cells across the core, 10 % at 10, converging at first order): a discretization error, not a cancellation.

2D magnetostatics uses ``AVA2D_nodal`` (planar or axisymmetric): the out-of-plane vector potential is the 2D scalar
unknown, without a source field.

Validation
----------

The ``formulations`` suite of ``pyTEC3D_APDL_benchmarks`` checks every type against analytical solutions, against
the A, V-A types on the same mesh, and, where MAPDL has the formulation, against MAPDL
(:doc:`../user_manual/validation`):

* magnetostatics: a ball of permeability :math:`\mu_r` in a uniform field (3D, ``PHI3D_nodal`` and ``TPHI3D_edge``;
  2D planar and axisymmetric with ``AVA2D_nodal`` and ``TPHI2D_edge``), a uniformly magnetized ball, both with a B-H
  curve, a rotating source current in a ball (``PHI3D_nodal``), a coil on a closed iron frame (``PHI3D_nodal`` with Hs
  solved or given, ``TPHI3D_edge`` with a fixed-T coil, against A);
* eddy currents: a conducting ball in a field switched on (3D and axisymmetric), with a static step first; a ring
  with a cut (induced, imposed current, imposed voltage) or a filled hole; a fixed-T solenoid around a conducting
  core, linear and iron (against A, V-A);
* superconductors: power-law slabs against a 1D reference (n = 20, Kim model Jc(B) with ``power_law_b``), the AC loss
  of a slab and of a round wire with a transport current through a cut against Bean's and Norris's critical state.

The examples ``core_3D_phi``, ``ring_3D_tphi``, ``coil_core_3D_tphi``, ``wire_2D_tphi`` and ``bulk_2D_tphi_axisym``
(:doc:`../examples/formulations`) and their tutorials show the same constructions in complete analyses.

The conducting ball in a field switched on at :math:`t = 0` has the induced moment

.. math::

   m(t) = -2\pi a^3 H_0\,\frac{6}{\pi^2}\sum_{k\ge1}\frac{e^{-k^2\pi^2t/\tau}}{k^2}, \qquad \tau = \mu_0\sigma a^2,

is reproduced by ``AVA3D_edge`` and ``TPHI3D_edge``, the box carrying the exact exterior
field. The edge DOF of a Dirichlet array is the circulation of the field along the edge **from its lexicographically
smaller end point to the larger one** (x, then y, then z: libMesh's edge orientation), which is how the benchmarks
impose a non-zero tangential A (:func:`pytec3d.tree_gauge.edges`, :func:`pytec3d.tree_gauge.edge_circulations`).
