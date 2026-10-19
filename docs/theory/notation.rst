Notation
--------

The symbols of the theory, with the names they have in the code: the DOF fields and output fields (``fields.json``),
the parameters of the element types (``mat.``, ``rc.``, ``me.``, ``mn.``) and the mesh arrays.

Regions
~~~~~~~

* :math:`\Omega` the model (a 3D region, or the 2D cross-section of a planar or axisymmetric body), :math:`\Omega_e` a
  cell, :math:`\partial\Omega` the external boundary.
* :math:`\Omega_0` the regions without current (``NON_CONDUCTING``, ``SOFT_MAGNETIC``, ``HARD_MAGNETIC``),
  :math:`\Omega_s` the regions with an imposed current density (``STRANDED_CONDUCTOR``), :math:`\Omega_c` the massive
  conductors (``SOLID_CONDUCTOR``, ``CIRCUIT_COUPLED_CONDUCTOR``), :math:`\Omega_\text{coil}` a stranded coil
  (``STRANDED_COIL``). See :doc:`domains`.

Fields
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 18 14 50

   * - Symbol
     - Code
     - Unit
     - Meaning
   * - :math:`\mathbf{A}`
     - ``A``
     - Wb/m
     - Magnetic vector potential, :math:`\mathbf{B}=\nabla\times\mathbf{A}`. In 2D only its out-of-plane component
       :math:`A_z` (or :math:`A_\theta`), stored as component 2.
   * - :math:`V`
     - ``V``
     - V
     - Electric scalar potential of a 3D solid conductor; in 2D, the terminal voltage of a circuit-coupled conductor
       across its out-of-plane length.
   * - :math:`v=\int_0^tV\,dt`
     - ``V`` (``_timint``)
     - V s
     - Time-integrated potential of the ``_timint`` element types.
   * - :math:`\mathrm{VOLT}`, :math:`\mathrm{EMF}`
     - ``VOLT``, ``EMF``
     - V
     - Terminal voltage and induced voltage of a stranded coil.
   * - :math:`\mathbf{B}`, :math:`\mathbf{H}`
     - ``B``, ``H``
     - T, A/m
     - Flux density and magnetic field (output fields, nodal and element).
   * - :math:`\mathbf{J}`
     - ``Jt``
     - A/m\ :sup:`2`
     - Total current density (output, element field).
   * - :math:`\mathbf{J}_s`
     - ``Js``
     - A/m\ :sup:`2`
     - Imposed current density of a stranded conductor (cell load).
   * - :math:`I_k`
     - ``AMPS``
     - A
     - Current injected at node :math:`k` into the equation of :math:`V` (nodal load).
   * - :math:`i`
     - ``CURR``
     - A
     - Current of a stranded coil, load of its VOLT equation (nodal load on the master node).
   * - :math:`\psi = rA_\theta`
     - ``PSI``
     - Wb/rad
     - Poloidal flux per radian (plasma algorithm).

The name of a Dirichlet value is the name of the field with a leading underscore (``_A``, ``_V``, ``_VOLT``,
``_EMF``); NaN means free.

Material and element data
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Symbol
     - Code
     - Meaning
   * - :math:`\mu_r`, :math:`\nu=(\mu_0\mu_r)^{-1}`
     - ``mat.mur``
     - Relative permeability (3x3 tensor) and reluctivity tensor; :math:`\mu_0=4\pi\cdot10^{-7}` H/m.
   * - :math:`\rho`, :math:`\sigma=\rho^{-1}`
     - ``mat.rho``
     - Resistivity [Ohm m] (3x3 tensor) and conductivity tensor.
   * - :math:`\mathbf{M}`
     - ``mat.M``
     - Magnetization [A/m]: :math:`\mathbf{B}=\mu_0\mu_r(\mathbf{H}+\mathbf{M})`.
   * - :math:`t`
     - ``rc.thickness``
     - Out-of-plane thickness factor of 2D cells (default 1).
   * - :math:`L`
     - (derived)
     - Out-of-plane length of a 2D cell: :math:`L=t` (planar), :math:`L=2\pi r\,t` (axisymmetric).
   * - :math:`N, S, R, \mathbf{T}, s`
     - ``rc.turns``, ``rc.cross_section``, ``rc.resistance``, ``rc.T``, ``rc.symmetry``
     - Turns, cross-section, resistance and winding direction of a coil; :math:`s` the ratio between the whole coil and
       the modelled part.
   * - :math:`V_c`
     - ``rc.volume`` (3D), :math:`2\pi\,\mathrm{radius}\,S\,t` (2D)
     - Volume of the whole coil.
   * - :math:`\mathsf{R}`
     - ``COORDINATE_SYSTEM``
     - Rotation of the element coordinate system of a cell, in which its vector and tensor parameters are given.
   * - :math:`\mathsf{R}_k`
     - ``node_angles``
     - Rotation of the nodal frame of node :math:`k`, in which its Dirichlet values and couplings are given.

Discretization
~~~~~~~~~~~~~~

* :math:`N_i` scalar Lagrange shape functions, :math:`\mathbf{W}_k` Nedelec (edge) shape functions; :math:`\mathbf{w}`
  the test function of :math:`\mathbf{A}`, :math:`q` that of a scalar field.
* :math:`\mathbf{a}`, :math:`\mathbf{V}` (or :math:`\mathbf{v}`), ... the vectors of DOFs; :math:`\mathbf{u}` all the
  DOFs of the model.
* :math:`\mathbf{K}`, :math:`\mathbf{C}`, :math:`\mathbf{M}` the stiffness, damping and mass matrices, :math:`\mathbf{F}`
  the right-hand side; :math:`K_{AA}`, :math:`C_{AV}`, ... their blocks between the fields.
* :math:`t^n`, :math:`\Delta t`, :math:`\theta` time levels, time step and parameter of the theta method; :math:`k` the
  Newton iteration, :math:`\mathbf{R}` the residual, :math:`\mathbf{J}` the Jacobian (from the context, not the
  current density).
