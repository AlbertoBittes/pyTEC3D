Circuits and Conductors
=======================

The domains that carry a net current through terminals: stranded coils with their circuit equation, conductors with a
terminal in 2D, currents injected into 3D conductors, and the time-integrated variants of their unknowns. The treatment
of stranded and massive conductors follows the standard potential formulations of the literature [Meunier2008]_.

.. _theory-stranded-coil:

Stranded coils
--------------

A **stranded coil** (``STRANDED_COIL``) is a winding of :math:`N` thin insulated turns of total cross-section
:math:`S`, wound along the unit direction :math:`\mathbf{T}` and fed through its two terminals. No eddy currents flow
inside the winding and the current density is uniform over its cross-section,

.. math::

   \mathbf{J} = \frac{N}{S}\, i\, \mathbf{T},

with :math:`i` the current of one turn. The current follows from the **circuit equation** of the coil,

.. math::

   \mathrm{VOLT} = R\, i + \mathrm{EMF},
   \qquad
   \mathrm{EMF} = \frac{d\Psi}{dt},
   \qquad
   \Psi = \frac{N}{S}\int_{\Omega_\text{coil}} \mathbf{T}\cdot\mathbf{A}\; d\Omega,

where VOLT is the terminal voltage, :math:`R` the resistance and :math:`\Psi` the flux linked with the winding (the
circulation of :math:`\mathbf{A}` along each turn, averaged over the cross-section, times the number of turns).

Unknowns and weak form
~~~~~~~~~~~~~~~~~~~~~~

The circuit unknowns ``VOLT`` and ``EMF`` are scalar Lagrange fields on the coil cells, tied to one value per coil by a
DOF coupling (:doc:`../user_manual/couplings`). With the coefficients of the element types

.. math::

   n = \frac{N}{S\,R},
   \qquad
   k = \frac{1}{R\, V_c},
   \qquad
   s = \text{symmetry},

where :math:`V_c` is the volume of the **whole** coil (``volume`` in 3D, :math:`2\pi\,\text{radius}\,S\,t` in 2D
axisymmetry, :math:`S\,t` planar) and :math:`s` the ratio between the whole coil and the modelled part (e.g. 12 for a
30 degree sector), the weak form of the coil is (the A-VOLT-EMF formulation of [MAPDLTheory2024]_)

.. math::

   \int_\Omega \nabla\times\mathbf{w}\cdot\nu\nabla\times\mathbf{A}
   - \int_{\Omega_\text{coil}} n\,\mathbf{w}\cdot\mathbf{T}\,(\mathrm{VOLT}-\mathrm{EMF}) = \dots
   \qquad\text{(A)}

.. math::

   \int_{\Omega_\text{coil}} q\, s\,\big(k\,\mathrm{VOLT} - n\,\mathbf{T}\cdot\partial_t\mathbf{A}\big)
   = q_m\,\mathrm{CURR}
   \qquad\text{(VOLT)}

.. math::

   \int_{\Omega_\text{coil}} q\,\big(n\,\mathbf{T}\cdot\partial_t\mathbf{A} - k\,\mathrm{EMF}\big) = 0
   \qquad\text{(EMF)}

with :math:`q_m` the test function of the master node, where the nodal load ``CURR`` is applied. Because VOLT and EMF
are coupled, :math:`q=1` on the coil, and with the modelled volume :math:`V_c/s`:

* the EMF row gives :math:`\mathrm{EMF} = s\,\frac{N}{S}\int\mathbf{T}\cdot\partial_t\mathbf{A} = d\Psi/dt` of the
  whole coil;
* the VOLT row gives :math:`(\mathrm{VOLT}-\mathrm{EMF})/R = \mathrm{CURR}`: the nodal load is the coil current
  :math:`i`, and the circuit equation holds;
* in the A row the coupling term is the source current density :math:`n(\mathrm{VOLT}-\mathrm{EMF})\,\mathbf{T} =
  (N/S)\,i\,\mathbf{T}`.

The coil is driven by its **current** (``CURR`` on the master node, VOLT free) or by its **terminal voltage**
(``_VOLT`` fixed, ``CURR = 0``); the excitation can change from one step to the next.

Block form
~~~~~~~~~~

With the coil unknowns :math:`(\mathbf{a}, \mathrm{VOLT}, \mathrm{EMF})`:

.. math::

   \begin{bmatrix} K_{AA} & K_{AU} & -K_{AU} \\ 0 & s\,K_{UU} & 0 \\ 0 & 0 & -K_{UU} \end{bmatrix}
   \begin{bmatrix} \mathbf{a} \\ \mathrm{VOLT} \\ \mathrm{EMF} \end{bmatrix}
   +
   \begin{bmatrix} 0 & 0 & 0 \\ s\,K_{AU}^T & 0 & 0 \\ -K_{AU}^T & 0 & 0 \end{bmatrix}
   \begin{bmatrix} \dot{\mathbf{a}} \\ \dot{\mathrm{VOLT}} \\ \dot{\mathrm{EMF}} \end{bmatrix}
   =
   \begin{bmatrix} \mathbf{F}_A \\ \mathrm{CURR} \\ 0 \end{bmatrix},
   \qquad
   (K_{AU})_{i} = -\int n\,\boldsymbol\phi_i\cdot\mathbf{T}\,q,
   \quad
   K_{UU} = \int k\, q\, q.

The coupling is in K in the A rows and in C in the circuit rows: the system is not symmetric. With the time-integrated
unknowns of the ``_timint`` element types every coupling term multiplies a time derivative, the whole coil block moves
to C and becomes symmetric for :math:`s = 1` (:ref:`theory-timint`).

.. _theory-solid-2d:

Conductors with a terminal in 2D
--------------------------------

In 2D the current of a massive conductor flows out of plane. A ``SOLID_CONDUCTOR`` has no applied voltage along that
direction: :math:`J_z=-\sigma\,\partial_tA_z`, the current induced in a closed ring. A
``CIRCUIT_COUPLED_CONDUCTOR`` adds a **terminal**: the voltage :math:`V` across its out-of-plane length, with
:math:`J_z=\sigma(V/L-\partial_tA_z)` (:ref:`theory-2d`). Coupled over the cross-section (one :math:`V` per
conductor), its equation is

.. math::

   \int_S J_z\, dS = \mathrm{AMPS}:

the total current through the cross-section equals the nodal load of the master node. Three ways of driving it:

* **current driven**: ``AMPS`` on the master node, :math:`V` free; :math:`V` is the voltage the source must supply;
* **voltage driven**: ``_V`` fixed, ``AMPS = 0``; the current is computed;
* **shorted**: ``_V = 0``; only induced currents flow, as in a solid conductor.

Several circuit-coupled conductors can be connected in series by coupling their :math:`V` with the same CP id (one
current, the sum of the voltages), and an external circuit can be closed by updating ``_V`` or ``AMPS`` from step to
step.

.. _theory-circuit-conductor:

Injected currents in 3D conductors
----------------------------------

In 3D every ``SOLID_CONDUCTOR`` has the potential :math:`V` and carries :math:`\mathbf{J}=-\sigma(\partial_t\mathbf{A}+
\nabla V)` with :math:`\nabla\cdot\mathbf{J}=0` (:ref:`theory-3d`). Its weak form, with the point currents :math:`I_k`
of the nodal load ``AMPS``,

.. math::

   \int_{\Omega_c}\nabla q\cdot\sigma\,(\partial_t\mathbf{A}+\nabla V)\,d\Omega = \sum_k I_k\, q(\mathbf{x}_k),

makes ``AMPS`` a current source at the node: :math:`I_k>0` enters the conductor there. A net current is imposed
through a conductor by injecting :math:`+I` at the nodes of one terminal and either :math:`-I` at the other terminal
or :math:`V` fixed there; a terminal with all its nodes coupled (one :math:`V`) is an equipotential electrode. The
potential of a conductor without fixed :math:`V` is defined up to a constant: fix it on one node (``_V``).

The halo currents of the plasma algorithm are injected this way into the wall (:ref:`theory-plasma-halo`). A load on a
node coupled to another one (a slave node) is carried to its master by the elimination of the coupling: the two loads
add up.

.. _theory-timint:

Time-integrated electric unknowns
---------------------------------

The ``_timint`` element types replace every electric unknown :math:`U` (the potential :math:`V` of the conductors, the
VOLT and EMF of the coils) by its time integral :math:`u=\int_0^tU\,dt`, so that :math:`U=\partial_tu`. Every term
that couples :math:`\mathbf{A}` with an electric unknown, and every term of the electric rows, then multiplies a time
derivative and moves to C:

.. math::

   \begin{bmatrix} C_{AA} & C_{Au} \\ C_{Au}^T & C_{uu} \end{bmatrix}
   \begin{bmatrix} \dot{\mathbf{a}} \\ \dot{\mathbf{u}} \end{bmatrix}
   + \begin{bmatrix} K_{AA} & 0 \\ 0 & 0 \end{bmatrix}
   \begin{bmatrix} \mathbf{a} \\ \mathbf{u} \end{bmatrix} = \mathbf{F}.

The pair (K, C) is symmetric (for the coils when :math:`s=1`), and so is the matrix :math:`K + C/(\theta\Delta t)` of
every step. For the user:

* Dirichlet values of the electric unknowns are time integrals ([V s], or [Wb] for the EMF of a coil); the nodal loads
  ``AMPS`` and ``CURR`` are still currents;
* on a **static** step (no C) the electric unknowns have no equation and are fixed to zero: no conduction or coil
  current flows on a static step;
* the solution is the same as with the base type for :math:`\theta=1` (backward Euler); for :math:`\theta<1` the two
  discretizations of the electric terms differ slightly.

``_timint`` domains: ``CIRCUIT_COUPLED_CONDUCTOR`` and ``STRANDED_COIL`` in ``AV2D_nodal_timint``,
``SOLID_CONDUCTOR`` and ``STRANDED_COIL`` in ``AV3D_nodal_timint`` and ``AV3D_edge_timint`` (:ref:`el-timint`).
