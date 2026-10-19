Implementation in pyTEC3D
=========================

How the formulations of the previous chapters become the discrete system the code solves: the kernels of the element
library, the time scheme, the constraints and the output fields.

.. _theory-kernels:

From the element library to the discrete system
-----------------------------------------------

The element types of :doc:`../user_manual/element_library` are lists of **kernels**. A kernel of the DOF system is

.. math::

   \int_{\Omega_e} \mathcal{T}(\boldsymbol{\phi}_i)^T\; P\; \mathcal{L}\big(\partial_t^d\,\mathbf{u}\big)\; d\Omega,

with :math:`\mathcal{T}` and :math:`\mathcal{L}` the test and trial operators (identity, ``grad``, ``curl`` or
``div``), :math:`P` a coefficient (a parameter, an expression or a law output: scalar, vector or tensor), and :math:`d`
the number of ``dt()`` on the trial. In ``element_types.json`` it reads, for example,
``{"trial": "grad(V)", "coeff": "sigma", "test": "A", "domain": ["SOLID_CONDUCTOR"]}``. A kernel goes into the stiffness
K for :math:`d=0`, the damping C for :math:`d=1` and the mass M for :math:`d=2`; a kernel without trial is a
right-hand side term :math:`\int\mathcal{T}(\boldsymbol\phi_i)^T P`. The semi-discrete system of a model is

.. math::

   \mathbf{M}\,\ddot{\mathbf{u}} + \mathbf{C}\,\dot{\mathbf{u}} + \mathbf{K}\,\mathbf{u} = \mathbf{F}(t),

with :math:`\mathbf{u}` all the DOFs (A, V, VOLT, EMF). Each section of the theory gives the kernels of a domain: the
2D ones in :ref:`theory-2d`, the 3D ones in :ref:`theory-3d`, the coils in :ref:`theory-stranded-coil`.

* **Quadrature**: Gauss rule of the default order of the finite elements of the element type.
* **2D weights**: the weight of each quadrature point is multiplied by the out-of-plane length :math:`L`
  (``thickness``, times :math:`2\pi r` in axisymmetry), so every 2D integral is an integral over the 3D body.
* **Coefficients**: material properties and real constants are constant over a cell (or computed per quadrature point
  by a law). Vector and tensor parameters are rotated from the element coordinate system of the cell to the global
  frame (:ref:`theory-element-cs`).
* **Element matrices** K, C (and M) and the operators of the output fields are computed at the first solve of a
  session and stored in the folder ``emat`` of the working directory; the following steps reuse them (the
  coefficients do not change in time). The right-hand side is assembled at every step from the current loads; nodal
  loads (``CURR``, ``AMPS``) are added to the row of the DOF they load.

.. _theory-time-integration:

Time integration
----------------

Each step solves the semi-discrete system at the end time :math:`t^{n+1}` of the step, writing the velocity
:math:`\mathbf{v}=\dot{\mathbf{u}}` and the acceleration :math:`\mathbf{a}=\ddot{\mathbf{u}}` as linear functions of the
unknown :math:`\mathbf{u}^{n+1}` and of the state :math:`(\mathbf{u}^n,\mathbf{v}^n,\mathbf{a}^n)` at the start of the
step: :math:`\mathbf{v}^{n+1} = c_v\,\mathbf{u}^{n+1} + \mathbf{h}_v`,
:math:`\mathbf{a}^{n+1} = c_a\,\mathbf{u}^{n+1} + \mathbf{h}_a`.

**Theta method** (first-order formulations, all the element types of the library):

.. math::

   \mathbf{v}^{n+1} = \frac{\mathbf{u}^{n+1}-\mathbf{u}^n}{\theta\,\Delta t} - \frac{1-\theta}{\theta}\,\mathbf{v}^n,
   \qquad c_v = \frac{1}{\theta\,\Delta t},\quad c_a = 0.

:math:`\theta=1` is backward Euler (first order, strongly damped), :math:`\theta=1/2` Crank-Nicolson (second order,
no numerical damping: stiff components may oscillate from step to step); values in between trade accuracy for
damping. ``time_step["theta"]`` must be in :math:`(0, 1]`.

**Newmark** (second-order formulations, ``time_order = 2``; parameters ``beta`` and ``gamma``):

.. math::

   \mathbf{a}^{n+1} = \frac{\mathbf{u}^{n+1}-\mathbf{u}^n-\Delta t\,\mathbf{v}^n}{\beta\,\Delta t^2}
                      - \Big(\frac{1}{2\beta}-1\Big)\mathbf{a}^n,
   \qquad
   \mathbf{v}^{n+1} = \mathbf{v}^n + \Delta t\big((1-\gamma)\,\mathbf{a}^n+\gamma\,\mathbf{a}^{n+1}\big);

:math:`\gamma=1/2`, :math:`\beta=1/4` (the defaults) is the unconditionally stable average acceleration scheme.

The **effective system** of a step is

.. math::

   \big(\mathbf{K} + c_v\,\mathbf{C} + c_a\,\mathbf{M}\big)\,\mathbf{u}^{n+1}
   = \mathbf{F}(t^{n+1}) - \mathbf{C}\,\mathbf{h}_v - \mathbf{M}\,\mathbf{h}_a,

and :math:`\mathbf{v}^{n+1}`, :math:`\mathbf{a}^{n+1}` are updated from :math:`\mathbf{u}^{n+1}` after the solve (with a
nonlinear law, :doc:`nonlinear`).

* A **static** step (``transient_effects = False``) uses :math:`c_v=c_a=0` and :math:`\mathbf{v}=\mathbf{a}=0`: C and
  M drop out. DOFs with no stiffness row on any cell around them (the time-integrated electric unknowns) are fixed to
  zero, since their equation would be empty.
* The first step starts from :math:`\mathbf{u}=\mathbf{v}=\mathbf{a}=0`: a transient first step with a non-zero load
  is a step change of the load at :math:`t=0`. Start with a static step to begin from a magnetostatic state.
* The effective matrix and its factorization are **reused** while :math:`\Delta t`, :math:`\theta`, :math:`\beta`,
  :math:`\gamma`, ``transient_effects`` and the set of constrained DOFs do not change; a change of any of them rebuilds
  and refactorizes it.

.. _theory-constraints:

Constraints: Dirichlet values, nodal frames and couplings
---------------------------------------------------------

**Dirichlet values** are read at every step from the point arrays ``_F`` of the DOF fields: a number fixes the
component, NaN leaves it free. The constrained DOFs are eliminated from the system (heterogeneous constraints), so
their equation is replaced by the imposed value.

**Nodal frames.** The values of a vector field at a node are given in the frame of the node, defined by the angles
``node_angles`` (:math:`\theta_{xy},\theta_{yz},\theta_{zx}`, rotations about z, x and y applied in this order) through
the rotation :math:`\mathsf{R}_k` of node :math:`k`: :math:`\mathbf{A}_\text{global}=\mathsf{R}_k\mathbf{A}_\text{local}`.
When all three local components are fixed the global values follow directly; when only some are, the solver writes
the corresponding linear constraints between the global components. The solution arrays are always in global
components. The rotations of the nodal frames, like those of the element coordinate systems, are computed **once, when
the model is loaded** into a solver session: ``node_angles`` (and ``COORDINATE_SYSTEM``) must be set before
``load_model_into_solver()``; the solver stops with an error if they are changed afterwards.

**DOF couplings.** The nodes with the same CP id in a point array ``CP_k`` share the values of the fields listed in
``dof_couplings``: the node with the lowest id is the master and every other node a slave,
:math:`\mathbf{u}_\text{slave}=\mathbf{u}_\text{master}` (for vector fields, component by component in the nodal frames
of the two nodes: how cyclic symmetry is written). Edge DOFs are coupled through their midside nodes, with the relative
orientation of the two edges. The couplings are set up when the model is loaded. A load on a slave node is carried to
its master by the elimination: the loads of the coupled nodes add up.

.. _theory-postprocessing:

Output fields
-------------

After each solve the derived fields of the element types are computed from the solution:

* **Nodal fields** (``B``, ``H`` on the points) are :math:`L^2` projections on discontinuous Lagrange elements of the
  order of the model, cell by cell, :math:`\mathbf{f}_e = M_e^{-1}\int_e\boldsymbol\phi\cdot P\,\mathcal{L}(\mathbf{u})`
  with :math:`M_e` the mass matrix of the cell; the values of the cells sharing a node are averaged at the node.
* **Element fields** (``B``, ``H``, ``Jt`` on the cells) are the plain mean of the values
  :math:`P\,\mathcal{L}(\mathbf{u})` at the quadrature points of the cell (not weighted by the quadrature weights).
  Where a field varies strongly inside a cell (the current density of a coarse 2D conductor) the mean is a poor
  estimate of the integral: compute integral quantities (currents, losses, forces) from the nodal fields or on a finer
  mesh.
* The terms of a field that depend on the solution are linear operators computed with the element matrices; the terms
  without trial (a load or a parameter alone, e.g. ``Jt = Js`` of a stranded conductor) are evaluated at every output
  from the current loads.
* Where a law applies, its secant is used (e.g. :math:`\mathbf{H}=\nu_s(B)\,\mathbf{B}`).
* ``H`` is :math:`\nu\nabla\times\mathbf{A}`: in a magnetized cell the magnetization is not subtracted.
* The time derivatives of the output (e.g. :math:`\partial_t\mathbf{A}` in the eddy current density) are the **mean
  rates over the step**, :math:`(\mathbf{u}^{n+1}-\mathbf{u}^n)/\Delta t`, for first-order formulations, so that
  :math:`\int\mathbf{J}\,dt` is consistent over the steps. This equals the rate of the discrete equations only for
  :math:`\theta=1`; for :math:`\theta<1` the eddy current densities of the output differ from those of the
  equations, by a large amount after sudden changes of the excitation (:doc:`../user_manual/outputs`). Second-order
  formulations use :math:`\mathbf{v}^{n+1}`.
