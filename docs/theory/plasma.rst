.. _theory-plasma:

Plasma Equilibrium and its Coupling with the Field
==================================================

pyTEC3D computes the electromagnetic transients of a tokamak in which the plasma is treated as a **sequence of
axisymmetric MHD equilibria**, coupled with the eddy currents of the surrounding conductors, which can be 2D
(axisymmetric) or 3D. This chapter derives the equations of the plasma model, the free-boundary equilibrium problem
and its iterative solution, and explains how the plasma is coupled, at every time step, with the finite-element
solution of :doc:`formulation_2d` and :doc:`formulation_3d`. The approach follows the APEC code developed at Fusion
for Energy [APEC2018]_ [APEC2019]_ [APEC2020]_, validated against the CarMa0NL code [Villone2013]_; the implementation in
pyTEC3D is the class
:class:`pytec3d.algorithm.Algorithm` (APEC-M), described for the user in :doc:`../user_manual/plasma_algorithm`.

.. contents:: On this page
   :local:
   :depth: 2

Axisymmetric equilibrium
------------------------

Ideal MHD equilibrium
~~~~~~~~~~~~~~~~~~~~~

On the time scales of the plasma evolution and of the eddy currents (milliseconds) the plasma is in force balance at
every instant: inertia is negligible and the ideal MHD equilibrium holds,

.. math::

   \mathbf{J}\times\mathbf{B} = \nabla p, \qquad \nabla\times\mathbf{B} = \mu_0\mathbf{J}, \qquad \nabla\cdot\mathbf{B} = 0,

with :math:`p` the plasma pressure. Two consequences follow directly: :math:`\mathbf{B}\cdot\nabla p = 0` and
:math:`\mathbf{J}\cdot\nabla p = 0`; the magnetic field lines and the current lines lie on surfaces of constant
pressure, the **flux surfaces**.

Poloidal flux and field
~~~~~~~~~~~~~~~~~~~~~~~

In a cylindrical frame :math:`(R,\varphi,Z)` with :math:`Z` the axis of the machine (the global y axis of the
pyTEC3D models), an axisymmetric divergence-free field can be written as

.. math::

   \mathbf{B} = \nabla\psi\times\nabla\varphi + F\,\nabla\varphi
   = \underbrace{-\frac{1}{R}\frac{\partial\psi}{\partial Z}\,\mathbf{e}_R
     + \frac{1}{R}\frac{\partial\psi}{\partial R}\,\mathbf{e}_Z}_{\mathbf{B}_p}
     + \underbrace{\frac{F}{R}\,\mathbf{e}_\varphi}_{B_\varphi\mathbf{e}_\varphi},

where :math:`\psi(R,Z)` is the **poloidal flux function** (the poloidal flux per radian through a circle of radius
:math:`R` at height :math:`Z`) and :math:`F(R,Z) = R B_\varphi` the **poloidal current function**
(:math:`2\pi F/\mu_0` is the poloidal current flowing through that circle). In terms of the toroidal component of the
vector potential,

.. math::

   \psi = R\,A_\varphi ,

which is how pyTEC3D computes the flux from the solution (``model_edit.poloidal_flux``: point array ``PSI`` from ``A``
and cell array ``PSIe``; in a 3D model :math:`A_\varphi` is first averaged along the torus, see
:ref:`theory-plasma-3d`).

The Grad-Shafranov equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Force balance and axisymmetry give :math:`p = p(\psi)` and :math:`F = F(\psi)`: pressure and poloidal current function
are flux functions. Ampere's law then splits into a toroidal and a poloidal part. The poloidal current density is
parallel to the poloidal field,

.. math::

   \mathbf{J}_p = \frac{F'(\psi)}{\mu_0}\,\mathbf{B}_p ,

and the toroidal current density is

.. math::

   J_\varphi = R\,p'(\psi) + \frac{F(\psi)F'(\psi)}{\mu_0 R} .

With :math:`\mu_0 J_\varphi\mathbf{e}_\varphi = \nabla\times\mathbf{B}_p` this gives the **Grad-Shafranov equation**
[Shafranov1958]_

.. math::

   \Delta^*\psi \equiv R\frac{\partial}{\partial R}\Big(\frac{1}{R}\frac{\partial\psi}{\partial R}\Big)
   + \frac{\partial^2\psi}{\partial Z^2} = -\mu_0 R\,J_\varphi
   = -\mu_0R^2p'(\psi) - F(\psi)F'(\psi)

inside the plasma, and :math:`\Delta^*\psi = -\mu_0 R J_\varphi` with the prescribed or induced current density
elsewhere (coils, passive structures; :math:`J_\varphi = 0` in vacuum).

Relation with the A formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 2D axisymmetric finite elements of the A, V-A formulation [BiroPreis1989]_ solve, for :math:`\mathbf{A}=A_\varphi\mathbf{e}_\varphi` and a uniform
reluctivity :math:`\nu=1/\mu_0`,

.. math::

   \nabla\times(\nu\nabla\times\mathbf{A}) = J_\varphi\mathbf{e}_\varphi
   \quad\Longleftrightarrow\quad
   -\frac{1}{\mu_0}\Big(\nabla^2A_\varphi - \frac{A_\varphi}{R^2}\Big) = J_\varphi ,

and since :math:`\Delta^*\psi = R(\nabla^2A_\varphi - A_\varphi/R^2)` for :math:`\psi = RA_\varphi`, this **is** the
Grad-Shafranov operator. The equilibrium is therefore obtained with the ordinary eddy-current machinery: the plasma
current density is written into the source current density ``Js`` of the plasma region (a ``STRANDED_CONDUCTOR``
domain), and the field of the plasma, of the coils and of the conductors is computed together, with any material
(e.g. saturating iron) and any boundary condition of the code. The poloidal plasma current :math:`\mathbf{J}_p` is
imposed as well (components R and Z of ``Js``): it produces the toroidal-field perturbation (diamagnetic or
paramagnetic) of the plasma and, in 3D, the poloidal halo currents that close through the structures.

Profiles of the plasma current
------------------------------

The functions :math:`p(\psi)` and :math:`F(\psi)` are not known: they are described by a few parameters chosen to
reproduce global quantities of the plasma. pyTEC3D uses the parametrization of APEC [APEC2018]_ [APEC2020]_,
generalized to two exponents:

.. math::

   J_\varphi = \lambda\Big[\beta\,\frac{R}{R_0} + (1-\beta)\,\frac{R_0}{R}\Big]\,g(\bar\psi),
   \qquad g(\bar\psi) = \big(1-\bar\psi^{\,\alpha_m}\big)^{\alpha_n},
   \qquad \bar\psi = \frac{\psi-\psi_\text{axis}}{\psi_b-\psi_\text{axis}} \in [0,1],

with :math:`\psi_\text{axis}` the flux at the magnetic axis, :math:`\psi_b` the flux at the plasma boundary (or at the
halo boundary, see :ref:`theory-plasma-halo`), :math:`R_0` the major radius (``R0``), and the parameters
``beta_p`` (:math:`\beta`), ``alpha_mp`` (:math:`\alpha_m`), ``alpha_np`` (:math:`\alpha_n`). Comparing with the
Grad-Shafranov current density identifies the two flux functions:

.. math::

   p'(\psi) = \frac{\lambda\beta}{R_0}\,g(\bar\psi), \qquad
   F F'(\psi) = \mu_0\lambda(1-\beta)R_0\,g(\bar\psi).

* :math:`\beta` weights the pressure-driven part of the current (the term in :math:`R/R_0`, larger on the outboard
  side) against the force-free part (:math:`R_0/R`); it is close to the **poloidal beta**
  :math:`\beta_p = \int p\,dV\big/\int \frac{B_p^2}{2\mu_0}dV` of the plasma.
* the exponents shape the radial profile: a peaked profile (low :math:`\alpha`) concentrates the current on the axis
  and increases the **internal inductance** :math:`l_i = 2L_i/(\mu_0R_0)`, with
  :math:`\tfrac12 L_i I_p^2 = \int_\text{plasma} \frac{B_p^2}{2\mu_0}dV`. The original APEC parametrization uses a single
  exponent, :math:`g = 1-\bar\psi^{\alpha}`, i.e. :math:`\alpha_m=\alpha`, :math:`\alpha_n = 1`, with :math:`\alpha`
  initialized to :math:`1/l_i` [APEC2018]_.
* :math:`\lambda` scales the total toroidal current to the prescribed plasma current (``plasma_current``):

  .. math::

     \lambda = I_p \Big/ \int_{\Omega_p}\Big[\beta\frac{R}{R_0} + (1-\beta)\frac{R_0}{R}\Big]g(\bar\psi)\,dS ,

  the integral running over the plasma cross-section (the cells of the plasma region, with their areas).

The pressure and the poloidal current function follow by integration from the boundary, where :math:`p=0` and
:math:`F = F_0 = R_0B_{\varphi 0}` (``AF0``, the vacuum toroidal field):

.. math::

   p(\bar\psi) = \frac{\lambda\beta}{R_0}(\psi_\text{axis}-\psi_b)\int_{\bar\psi}^1 g(s)\,ds,
   \qquad
   F^2(\bar\psi) = F_0^2 + 2\mu_0\lambda(1-\beta)R_0(\psi_\text{axis}-\psi_b)\int_{\bar\psi}^1 g(s)\,ds .

For :math:`\alpha_m=1` the integral is :math:`(1-\bar\psi)^{\alpha_n+1}/(\alpha_n+1)`, and pyTEC3D evaluates

.. math::

   F(\bar\psi) = \sqrt{F_0^2 + \frac{2\mu_0R_0(1-\beta)(\psi_\text{axis}-\psi_b)\lambda}{\alpha_n+1}
   \big(1-\bar\psi^{\,\alpha_m}\big)^{\alpha_n+1}},

exact for :math:`\alpha_m = 1` and a close approximation for :math:`\alpha_m\approx1` (the default ``alpha_mp`` is
1.0001). For the single-exponent profile the closed form is

.. math::

   F(\bar\psi) = \Big\{F^2(\psi_b) + \xi\big[\bar\psi^{\,\alpha+1} - (\alpha+1)\bar\psi + \alpha\big]\Big\}^{1/2},
   \qquad \xi = \frac{2\lambda\mu_0R_0(\psi_\text{axis}-\psi_b)(1-\beta)}{\alpha+1},

as in [APEC2020]_. The poloidal current density is then

.. math::

   \mathbf{J}_p = \frac{F'}{\mu_0}\mathbf{B}_p = \frac{\lambda(1-\beta)R_0\,g(\bar\psi)}{F(\bar\psi)}\,\mathbf{B}_p ,
   \qquad \mathbf{B}_p = \Big(-\frac{1}{R}\frac{\partial\psi}{\partial Z},\ \frac{1}{R}\frac{\partial\psi}{\partial R}\Big),

with the gradient of :math:`\psi` averaged on each cell. pyTEC3D writes :math:`(J_R, J_Z, J_\varphi)` of each plasma
cell into ``Js`` (cell components 0, 1, 2 in 2D: R, Z and toroidal).

The free-boundary problem
-------------------------

The plasma boundary is not known in advance: it is the **last closed flux surface** (LCFS) of the total field, which
depends on the plasma current, which depends on the boundary. The problem is:

   find :math:`\psi` and the plasma region :math:`\Omega_p(\psi)` such that :math:`\Delta^*\psi = -\mu_0 R J_\varphi`
   everywhere, with :math:`J_\varphi` given by the profile in :math:`\Omega_p`, by the coil currents in the coils, by
   the induced currents in the conductors, and :math:`0` elsewhere, for prescribed :math:`I_p`, :math:`\beta`,
   :math:`\alpha_m`, :math:`\alpha_n` and coil currents.

Magnetic axis, X point, limiter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The search is done on the vacuum-chamber mesh (the 2D cells where the plasma can be, ``VACUUM_CHAMBER = 1``), with
the gradient and the Hessian of :math:`\psi` computed on its nodes:

* the **magnetic axis** is the extremum of :math:`\psi` (its maximum with the sign convention of the ITER examples);
* the **X point** is the saddle point of :math:`\psi` (:math:`\det\nabla\nabla\psi<0`) where :math:`|\nabla\psi|` is
  smallest;
* the **limiter point** is the node of the chamber wall with the largest :math:`\psi`, excluding the private-flux
  region of the divertor (where :math:`(\mathbf{r}-\mathbf{r}_\text{axis})\cdot\nabla\psi` has the sign opposite to
  that of the core).

The boundary flux is the one reached first when moving out from the axis: :math:`\psi_b = \psi_\text{limiter}` if it
is larger than :math:`\psi_X` (a **limiter** configuration, the plasma touches the wall), :math:`\psi_b = \psi_X`
otherwise (a **diverted**, X-point configuration). The plasma region is found in two steps: the cells inside the flux
surface of the X point that are connected to the axis, grown from the axis through cells sharing an edge (so that the
region does not leak through the X point into the private-flux region), then, among them, the cells inside the
boundary of the current, :math:`\psi_2` (:math:`\psi_b`, or the halo boundary when a halo is active): in a limiter
configuration the region ends at the limiter flux surface. These cells, marked in the cell array ``PLASMA``, are
exactly the cells where the profile, hence the plasma current, is not zero.

Picard iteration
~~~~~~~~~~~~~~~~

The non-linearity of the free-boundary problem (the region and :math:`\bar\psi` depend on :math:`\psi`) is resolved by
fixed-point (Picard) iterations, as in APEC [APEC2018]_:

#. from :math:`\psi^k`: axis, X point, limiter, boundary flux, plasma region :math:`\Omega_p^k`;
#. the current density :math:`\mathbf{J}^k(\psi^k)` in :math:`\Omega_p^k`, scaled to :math:`I_p` (:math:`\lambda^k`);
#. (equilibrium only) position control: corrections of the control-coil currents, see below;
#. a finite-element solve with the new source, giving :math:`\mathbf{A}^{k+1}` and :math:`\psi^{k+1}`;
#. convergence when :math:`\max_{\text{chamber nodes}}|\psi^{k+1}-\psi^k| <` ``tolerance_value`` (an absolute flux
   tolerance, e.g. 1e-4 Wb/rad), or failure after ``max_picard_iterations``.

The first iteration starts from an initial guess of the current: typically a uniform current density over the vacuum
chamber carrying :math:`I_p`. Each iteration appends a row to ``__pytec3d__/pytec3d_iteration_diagnostics.csv``
(axis, X point, limiter, fluxes, configuration, correction currents, :math:`\lambda`, plasma cells).

Position control
~~~~~~~~~~~~~~~~

A free-boundary equilibrium of an elongated plasma is vertically unstable, and a plain Picard iteration drifts: the
plasma moves up or down from one iteration to the next. As in a real machine, a feedback on the coil currents keeps
it in place. pyTEC3D measures the poloidal field :math:`\mathbf{B}_t = (B_{R,t}, B_{Z,t})` at the **target position of
the magnetic axis** (``target_central_axis_x``, ``_y``) and computes a radial and a vertical correction current from
the field that a unit correction current produces there,

.. math::

   \begin{bmatrix} b^\text{rad}_R & b^\text{ver}_R \\ b^\text{rad}_Z & b^\text{ver}_Z \end{bmatrix}
   \begin{bmatrix} I_\text{rad} \\ I_\text{ver} \end{bmatrix} = -\begin{bmatrix} B_{R,t} \\ B_{Z,t}\end{bmatrix},

with :math:`b^\text{rad}` = (``bradx_value``, ``brady_value``) and :math:`b^\text{ver}` = (``bverx_value``,
``bvery_value``). Each control coil :math:`c` of ``control_coils`` then carries its initial current plus
:math:`g^\text{rad}_c I_\text{rad} + g^\text{ver}_c I_\text{ver}` (the gains of the coil), written as a current density
added to its ``Js``. The correction is computed at every iteration from the field of the current state, which
already contains the previous correction, and it *replaces* the previous one (the initial coil currents are kept as
the base): the control is **proportional**, not integral. Writing :math:`M` for the matrix of the ``b*_value``
parameters and :math:`M_c` for the actual response of the control coils with their gains, the converged correction
satisfies :math:`M\mathbf I = -(\mathbf B_0 + M_c\mathbf I)`, with :math:`\mathbf B_0` the field at the target without
correction, so a residual field

.. math::

   \mathbf B_t = \big(\mathbb 1 + M_c M^{-1}\big)^{-1}\mathbf B_0

remains at the target (neglecting the reaction of the plasma, which moves with the field): the axis converges
*near* the target, closer for a larger loop gain :math:`M_cM^{-1}`. In ``tests/examples/ITER_2D`` the ``b*_value`` are the
fields of unit gains while the gains are 0.5 (:math:`M_c = M/2`), and the axis converges about 5 cm from the target
with a residual field of a few mT (see :doc:`../tutorials/tutorial_iter_equilibrium`). The control currents are a
small fraction of the coil currents (a few percent in [APEC2018]_). Control is used only for the equilibrium; during
the evolution the plasma moves freely.

.. _theory-plasma-evolution:

Plasma evolution coupled with the eddy currents
-----------------------------------------------

Quasi-static evolution
~~~~~~~~~~~~~~~~~~~~~~

During a disruption, a vertical displacement or a current quench the plasma passes through a sequence of equilibria driven by prescribed
waveforms of :math:`I_p(t)` and :math:`\beta_p(t)` (and possibly of the profile exponents), while the changing flux
induces currents in the vacuum vessel, the passive structures and the coils. Those currents change the field, hence
the equilibrium: the plasma and the conductors are coupled. Each time step :math:`t^n\to t^{n+1}` must find the
equilibrium at :math:`t^{n+1}` consistent with the eddy currents at :math:`t^{n+1}`:

.. math::

   \mathsf{K}\,\mathbf{u}^{n+1} + \mathsf{C}\,\dot{\mathbf{u}}^{n+1} = \mathbf{F}\big(\mathbf{J}_\text{plasma}(\psi^{n+1})\big),
   \qquad \psi^{n+1} = R\,A_\varphi(\mathbf{u}^{n+1}),

where the conduction terms of :math:`\mathsf{C}` carry the eddy currents (see :ref:`theory-time-integration`) and the
plasma enters only through the source :math:`\mathbf{F}`, the plasma being a ``STRANDED_CONDUCTOR`` region whose
current density is imposed by the profile.

Two-step scheme of APEC, implicit iteration of pyTEC3D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

APEC [APEC2019]_ solves each step with a two-step loop: a transient solve with the current plasma current density gives
the eddy currents of the structures; the plasma region and the profile are updated; a static solve (transient effects
off) with the plasma current and the eddy currents applied as imposed currents gives the new flux; the loop is repeated
until the flux converges.

pyTEC3D folds the two solves into one: at every Picard iteration it **solves the same time step again**,
:math:`t^{n+1}` with the history :math:`(\mathbf{u}^n,\dot{\mathbf{u}}^n)` unchanged, with the updated plasma current
density (the solver accepts a repeated solve at the same time, see :doc:`../user_manual/analysis`). The eddy currents
are therefore recomputed with the theta scheme from the latest plasma current at every iteration, and at convergence
:math:`(\mathbf{u}^{n+1},\psi^{n+1})` satisfies the coupled implicit equations above: the plasma and the structures
are advanced together, with no lag. Since the time step, the scheme and the constraints do not change between
iterations, the factorized matrix is reused and each iteration costs one forward-backward substitution.

The first step is a static solve (``transient_effects = False``, at a small time) that starts the evolution from the
equilibrium (its current densities are copied into ``Js`` of the evolution model); the following steps are transient.
In the evolution there is no position control: the coils keep their equilibrium currents (stranded conductors), or
follow their circuit equations if they are modelled as voltage-driven stranded coils.

Voltage-driven coils
~~~~~~~~~~~~~~~~~~~~

The poloidal field coils and the central solenoid react to a plasma transient: currents are induced in them, which
have a stabilizing effect on the plasma [APEC2020]_. Modelled as **stranded coils** (:ref:`theory-stranded-coil`), each
coil has the A-VOLT-EMF formulation,

.. math::

   \mathrm{VOLT} = R\,i + \frac{d\Psi}{dt}, \qquad \mathbf{J} = \frac{N}{S}\,i\,\mathbf{T},
   \qquad \Psi = \frac{N}{S}\int\mathbf{T}\cdot\mathbf{A}\,dV,

with its terminal voltage fixed to the value that sustains its equilibrium current, :math:`\mathrm{VOLT} = R\,i_0`.
During the transient the induced EMF changes the coil currents; a larger resistance reduces the stabilizing effect (the
limit :math:`R\to\infty` gives fixed currents). See ``tests/examples/ITER_2D_evol_coils`` and ``tests/examples/ITER_3D_evol_coils``.

.. _theory-plasma-3d:

Coupling the axisymmetric plasma with 3D structures
---------------------------------------------------

The main feature of the approach [APEC2019]_ is the coupling of the **axisymmetric** plasma equilibrium with **3D**
eddy currents in non-axisymmetric structures. The plasma problem stays on the 2D vacuum-chamber mesh; the field is
solved on a 3D mesh obtained by revolving the 2D mesh about the axis over a sector (``extrusion_angle``, with
``n_toroidal_divisions`` cells along the torus) with cyclic symmetry on its two faces. Two maps link the meshes:

* **projection**: the current density of each 2D plasma cell, :math:`(J_R, J_Z, J_\varphi)`, is written on the 3D cells
  of every toroidal station that correspond to it, converted to the global components of each station (the toroidal
  component along :math:`\mathbf{e}_\varphi` of the station). In 3D the cells of the plasma region use the global
  coordinate system (their ``Js`` holds global components);
* **toroidal average**: the flux used by the plasma is :math:`\psi = R\,\langle A_\varphi\rangle`, the toroidal component
  of the 3D vector potential averaged over the stations, on the nodes (``PSI``) and the cells (``PSIe``).

The plasma thus sees the axisymmetric part of the field of the 3D currents, while the structures carry fully 3D eddy
currents (with the periodicity of the sector). The vessel and the passive structures are ``SOLID_CONDUCTOR`` cells
with the electric potential :math:`V` (time-integrated in the ``_timint`` types, the nodal formulation used by APEC,
or the potential itself); :math:`V` must be fixed on one node of each floating conductor, and :math:`\mathbf{A}` and
:math:`V` are coupled between the sector faces. In the 2D evolutions the structures are ``SOLID_CONDUCTOR`` cells too,
without electric unknown: closed toroidal shells carrying :math:`J_\varphi=-\sigma\partial_tA_\varphi`.

.. _theory-plasma-halo:

Halo currents
-------------

In vertical displacement events the plasma touches the wall and part of its current flows in the **halo**, the region
of open field lines outside the last closed flux surface, closing through the wall [APEC2020]_. The model extends the
current-carrying region beyond the LCFS to a halo boundary flux

.. math::

   \psi_2 = \psi_b - \Delta_1\big(1+(t-t_0)\Delta_2\big)\big(\psi_\text{axis}-\psi_b\big), \qquad t \ge t_0,

with ``halo_delta1`` (:math:`\Delta_1`, the width of the halo as a fraction of the core flux), ``halo_delta2``
(:math:`\Delta_2`, its growth rate) and ``t0_halo`` (:math:`t_0`, typically the limiter / X-point transition). The
normalized flux is redefined with :math:`\psi_2`, :math:`\bar\psi = (\psi-\psi_\text{axis})/(\psi_2-\psi_\text{axis})`, so
that the same profile, toroidal and poloidal, applies in the core and in the halo, and :math:`\lambda` scales the total
(core + halo) current to :math:`I_p`.

The poloidal halo current that reaches the wall is **injected into the conductors**. Along the wall nodes wetted by the
halo (between the flux surfaces :math:`\psi_b` and :math:`\psi_2`), the poloidal current crossing the wall between two
consecutive nodes :math:`h` and :math:`h+1` follows from the poloidal current function,

.. math::

   I_h = \frac{2\pi}{\mu_0}\Big(F\big(\bar\psi_h\big) - F\big(\bar\psi_{h+1}\big)\Big)\,\frac{\Theta}{360^\circ},

with :math:`\Theta` the toroidal extent of the model; half of it is added to each of the two nodes as the nodal load
``AMPS`` of the potential :math:`V` of the conductor (a current source, :ref:`theory-circuit-conductor`). On the 3D
mesh this current of the sector is shared among the :math:`n+1` node stations of the sector with the trapezoid weights
:math:`1/(2n)` on the two end stations and :math:`1/n` on the inner ones, so that the sector receives it once: when the
two end stations are a periodic pair, the load of the slave node is carried to its master and the two halves add up.
The halo acts de facto as a current generator, driven by the equilibrium itself. Halo injection needs a 3D model with
conducting structures (the ``AMPS`` load of their potential).

Summary of the coupled model
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Ingredient
     - Model
   * - Plasma
     - Sequence of axisymmetric ideal MHD equilibria; current density from the Grad-Shafranov profile with prescribed
       :math:`I_p(t)`, :math:`\beta_p(t)`, exponents; imposed as a source current density (toroidal and poloidal).
   * - Plasma boundary
     - Last closed flux surface from the axis, X point and limiter of the toroidally averaged flux; halo extension from
       :math:`t_0`.
   * - Equilibrium
     - Picard iterations with position control on the field at the target axis.
   * - Evolution
     - Implicit coupling: Picard iterations repeat the transient solve of the step until the flux converges.
   * - Structures
     - 2D or 3D eddy currents (A-V formulation), any material, cyclic symmetry in 3D.
   * - Coils
     - Fixed currents (stranded conductors) or voltage-driven multi-turn coils (stranded coils, A-VOLT-EMF).
   * - Halo
     - Core + halo current profile; poloidal halo current injected into the wall as nodal currents.

Limits: the plasma has no inertia and no resistivity of its own (its current follows the prescribed waveform); the
profile is described by a few parameters; the plasma is axisymmetric even when the structures are not.
