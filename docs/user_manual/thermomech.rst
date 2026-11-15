Thermal and Elastic Analyses
============================

Element types
-------------

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Type
     - Domain
     - Formulation
   * - ``TH2D_nodal``, ``TH3D_nodal``
     - ``THERMAL``
     - Heat conduction: the temperature ``TEMP`` (Lagrange, order 1 or 2; 2D planar or axisymmetric),
       :math:`\int k\nabla T\cdot\nabla v + \rho c_p\,\dot T\,v = \int Q v`. Material ``kt`` (conductivity tensor,
       W/(m K)), ``rho_cp`` (volumetric heat capacity, J/(m³ K), default 0: steady); element data ``QV`` (W/m³);
       nodal heat ``HEAT`` (W). Outputs: ``TC`` (cell temperature), ``GRADT``, ``HFLUX`` = :math:`-k\nabla T` (cells;
       ``HFLUX`` also at the points).
   * - ``EL2D_nodal``, ``EL3D_nodal``
     - ``ELASTIC``
     - Linear isotropic elasticity: the displacement ``U`` (Lagrange, order 1 or 2; 2D plane strain, or plane stress
       with the real constant ``plane_stress`` = 1), :math:`\int 2\mu\,\varepsilon(u):\varepsilon(v) + \lambda\,
       \nabla\cdot u\,\nabla\cdot v = \int f\cdot v`. Material ``Ey`` (Young's modulus), ``nu_p`` (Poisson's ratio);
       element data ``FB`` (body force, N/m³); nodal forces ``FORCE`` (N). Outputs: the stress and strain rows ``SX``,
       ``SY``, ``SZ``, ``EX``, ``EY``, ``EZ`` (cells and points; 2D: the z rows hold only :math:`\sigma_{zz}`,
       :math:`\varepsilon_{zz}`) and the **tensors** ``SIGMA`` and ``EPS``.

The nonlinear law ``conductivity_t`` makes ``kt`` a function of the temperature, ``k(T)``: points ``T``, ``k``, or an
analytic law ``"k": "<expression of T>"`` with ``"dk"`` (:doc:`nonlinear`). Its Jacobian has the secant ``k(T)`` only
(the term :math:`k'(T)\,\delta T\,\nabla T` is left out), so Newton converges linearly there (12 iterations for k
varying by a factor 3 in ``thermal_nonlinear_slab``).

Tensor fields
-------------

The solver computes a tensor as its three rows (vector fields). After every ``solve()`` the API writes ``SIGMA``
and ``EPS`` as 9-component cell and point arrays, row-major (xx, xy, xz, yx, yy, yz, zx, zy, zz), the layout of VTK
tensors (ParaView shows their components, invariants and eigenvalues). ``pytec3d.tensor_fields.TENSOR_FIELDS`` lists
them.

Contact
-------

``pytec3d.contact`` couples bodies meshed apart with **coincident nodes** on their interface, through the DOF
couplings (:doc:`couplings`):

.. code-block:: python

   from pytec3d import contact

   pairs = contact.interface_pairs(mesh, nodes_of_body_1, nodes_of_body_2)   # (n, 2): coincident nodes
   contact.bonded(solver.model, pairs)              # linear: all the components of U tied
   contact.sliding(solver.model, pairs, "z")        # linear, frictionless: only the normal component tied
   result = contact.solve_unilateral(solver, pairs, "z")   # nonlinear: frictionless, no tension

``sliding`` ties ``"U[z]"``, the z component **in the nodal frame** of the first node: with ``node_angles`` the
normal may have any direction. A coupling ``"<vector field>[x|y|z]"`` couples one component in general.

``solve_unilateral`` is an active-set iteration of static solves: the active pairs tied in the normal direction; a
pair whose normal stress is tensile is released, an open pair whose nodes interpenetrate is closed; until the set does
not change (each iteration reloads the model, couplings being fixed at the load). It returns the active pairs, the
contact pressure :math:`-n\cdot\sigma\cdot n` at them and the gaps. Benchmarks: ``contact_blocks`` (bonded, sliding,
unilateral in compression, tension and under a tilted load: the contact conditions hold to round-off).

Verification
------------

``patch_thermal`` and ``patch_elastic`` (the patch tests of every shape and order), ``thermal_nonlinear_slab``,
``thermal_transient_slab``, ``elastic_lame`` and ``contact_blocks`` in the formulations benchmarks.
