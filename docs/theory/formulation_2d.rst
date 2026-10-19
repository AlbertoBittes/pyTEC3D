.. _theory-2d:

2D formulation: ``AV2D_nodal``
------------------------------

The 2D element type ``AV2D_nodal`` (:ref:`el-av2d`) models bodies that are invariant along a direction: **planar**
(prismatic, invariant along :math:`z`, ``option3 = 0``) or **axisymmetric** (invariant around the :math:`y` axis,
``option3 = 1``, with :math:`x = r`). The currents flow out of plane and the field lies in the plane.

Geometry: the out-of-plane length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 2D cell stands for a 3D body of out-of-plane length

.. math::

   L = t \ \ (\text{planar}),
   \qquad
   L = 2\pi r\, t \ \ (\text{axisymmetric}),

with :math:`t` the real constant ``thickness`` (default 1: a 1 m slice, or the full :math:`360^\circ` ring; another
value models a fraction of the ring). Every volume integral of the weak form is computed as
:math:`\int_{\Omega^\text{3D}}f\,dV=\int_\Omega f\,L\,dS`: the quadrature weights of the 2D cells are multiplied by
:math:`L`, so loads and results (currents, energies, forces) refer to the 3D body.

Unknowns
~~~~~~~~

**Vector potential.** Only the out-of-plane component, a scalar Lagrange field of order ``option2``:

.. math::

   \mathbf{A} = A_z(x,y,t)\,\mathbf{e}_z\ \ (\text{planar}),
   \qquad
   \mathbf{A} = A_\theta(r,z,t)\,\mathbf{e}_\theta\ \ (\text{axisymmetric}),

stored as component 2 of the point array ``A``. The flux density is in plane:

.. math::

   \mathbf{B} = \Big(\frac{\partial A_z}{\partial y},\; -\frac{\partial A_z}{\partial x}\Big)\ \ (\text{planar}),
   \qquad
   \mathbf{B} = \Big(-\frac{\partial A_\theta}{\partial z},\; \frac{1}{r}\frac{\partial (rA_\theta)}{\partial r}\Big)
   \ \ (\text{axisymmetric}).

In axisymmetry :math:`rA_\theta` is the poloidal flux per radian, :math:`\psi`: the lines :math:`\psi=` const are the
field lines (``PSI`` of the plasma algorithm, :doc:`plasma`).

**Terminal voltage** (``CIRCUIT_COUPLED_CONDUCTOR`` only). The electric field of an out-of-plane current is
:math:`E_z=-\partial_tA_z + V/L`, with :math:`V` the voltage applied across the out-of-plane length of the conductor
(between the two ends of a planar bar, across the gap of a ring). pyTEC3D writes :math:`V` as a scalar field whose
gradient is out of plane,

.. math::

   \nabla V = -\frac{V}{L}\,\mathbf{e}_z,

so that the 3D expression :math:`\mathbf{E}=-\partial_t\mathbf{A}-\nabla V` holds unchanged. Its nodal values are
coupled over the cross-section of each conductor (one value per conductor, :doc:`../user_manual/couplings`).

**Coil unknowns** (``STRANDED_COIL``): ``VOLT`` and ``EMF``, out-of-plane scalars coupled over each coil
(:ref:`theory-stranded-coil`).

Current density
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Domain
     - :math:`J_z` (output ``Jt``, component 2)
   * - ``NON_CONDUCTING``, ``SOFT_MAGNETIC``, ``HARD_MAGNETIC``
     - :math:`0`
   * - ``STRANDED_CONDUCTOR``
     - :math:`J_{s,z}` (the cell load ``Js``)
   * - ``SOLID_CONDUCTOR``
     - :math:`-\sigma\,\partial_tA_z`: the induced current of a closed, shorted conductor
   * - ``CIRCUIT_COUPLED_CONDUCTOR``
     - :math:`\sigma\,(V/L-\partial_tA_z)`
   * - ``STRANDED_COIL``
     - :math:`n\,T_z\,(\mathrm{VOLT}-\mathrm{EMF})=(N/S)\,i\,T_z`

with :math:`\sigma` the out-of-plane conductivity :math:`\sigma_{zz}`. A 2D solid conductor needs no unknown: its
net current is whatever the induction drives (in a planar model, the current of an infinitely long bar whose ends are
connected at infinity; the net current of all the solid conductors of a planar model need not be zero).

Weak form and element matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`A_h=\sum_jN_j\,a_j`, the test functions :math:`w_i=N_i\mathbf{e}_z` and the weight :math:`L` in every
integral (:math:`\int_e f = \sum_g \omega_g\,|J_g|\,L_g\,f(\mathbf{x}_g)`), the terms of ``AV2D_nodal`` are, domain by
domain:

.. list-table::
   :header-rows: 1
   :widths: 24 10 66

   * - Domain
     - Matrix
     - Term
   * - all
     - K
     - :math:`(K_{AA})_{ij}=\int\nabla\times w_i\cdot\nu\,\nabla\times w_j`
   * - STRANDED_CONDUCTOR
     - F
     - :math:`(F_A)_i=\int w_i\cdot\mathbf{J}_s`
   * - all but NON_CONDUCTING, SOFT_MAGNETIC
     - F
     - :math:`(F_A)_i \mathrel{+}= \int\nabla\times w_i\cdot\mathbf{M}` (zero when ``M`` is not given)
   * - SOLID, CIRCUIT_COUPLED
     - C
     - :math:`(C_{AA})_{ij}=\int w_i\cdot\sigma\,w_j`
   * - CIRCUIT_COUPLED
     - K
     - :math:`(K_{AV})_{ij}=\int w_i\cdot\sigma\,\nabla N_j`, :math:`(K_{VV})_{ij}=\int\nabla N_i\cdot\sigma\,\nabla N_j`
   * - CIRCUIT_COUPLED
     - C
     - :math:`(C_{VA})_{ij}=\int\nabla N_i\cdot\sigma\,w_j`
   * - STRANDED_COIL
     - K, C
     - the coil blocks of :ref:`theory-stranded-coil`

with :math:`\nabla N_j=-N_j/L\,\mathbf{e}_z` for the out-of-plane scalar :math:`V`. In the matrices of a
circuit-coupled conductor, :math:`(K_{AV})_{ij}=-\int\sigma_{zz}N_iN_j/L` and :math:`(K_{VV})_{ij}=\int\sigma_{zz}
N_iN_j/L^2`: the voltage :math:`V` acts as a uniform out-of-plane electric field :math:`V/L` over the conductor.

For the cells of a circuit-coupled conductor the element system reads

.. math::

   \begin{bmatrix} C_{AA} & 0 \\ C_{VA} & 0 \end{bmatrix}
   \begin{bmatrix} \dot{\mathbf{a}} \\ \dot{\mathbf{V}} \end{bmatrix}
   +
   \begin{bmatrix} K_{AA} & K_{AV} \\ 0 & K_{VV} \end{bmatrix}
   \begin{bmatrix} \mathbf{a} \\ \mathbf{V} \end{bmatrix}
   =
   \begin{bmatrix} \mathbf{F}_A \\ \mathbf{I} \end{bmatrix},

with :math:`\mathbf{I}` the nodal loads ``AMPS``. Coupled over the conductor, the test function of :math:`V` is 1 on
its cross-section and the :math:`V` row becomes

.. math::

   \int_S \sigma\Big(\frac{V}{L}-\partial_tA_z\Big)\,dS = \int_S J_z\,dS = I,

the total current through the cross-section equals the load ``AMPS`` of the master node. The conductor is therefore
**current driven** (``AMPS`` on the master node, :math:`V` computed), **voltage driven** (``_V`` fixed, the current is
the reaction) or **shorted** (``_V = 0``: the same as a solid conductor). See :ref:`theory-solid-2d`.

The axisymmetric operator
~~~~~~~~~~~~~~~~~~~~~~~~~

With isotropic :math:`\nu` and the weight :math:`L=2\pi r\,t`, the curl-curl term is the weak form of

.. math::

   -\frac{\partial}{\partial r}\Big(\frac{\nu}{r}\frac{\partial(rA_\theta)}{\partial r}\Big)
   - \frac{\partial}{\partial z}\Big(\nu\frac{\partial A_\theta}{\partial z}\Big) = J_\theta,

which, written for :math:`\psi=rA_\theta` with :math:`\nu=1/\mu_0`, is the operator :math:`\Delta^*\psi=-\mu_0rJ_\theta`
of the Grad-Shafranov equation (:ref:`theory-plasma`). The term :math:`A_\theta/r` of the flux density is integrated
by the quadrature like the rest; on the axis :math:`A_\theta=0` is imposed.

Notes
~~~~~

* No gauge is needed: the out-of-plane potential has no gradient freedom (:ref:`theory-gauge`).
* Element fields: ``B``, ``H`` (nodal and element) and ``Jt`` (element, the table above). The constant part of
  ``Jt`` (``Js`` of a stranded conductor) is evaluated at every output from the current loads.
* The time-integrated variant ``AV2D_nodal_timint`` uses :math:`v=\int V dt` for the circuit-coupled conductors and
  the time integrals of VOLT and EMF for the coils (:ref:`theory-timint`).
