.. _theory-3d:

3D formulations: ``AV3D_nodal`` and ``AV3D_edge``
-------------------------------------------------

The 3D element types solve the A-V formulation of :doc:`maxwell` on volume meshes, with two discretizations of the
vector potential: Lagrange (nodal) components in ``AV3D_nodal`` (:ref:`el-av3d-nodal`) and Nedelec edge elements in
``AV3D_edge`` (:ref:`el-av3d-edge`). The scalar potential and the coil unknowns are Lagrange fields in both.

Unknowns and currents
~~~~~~~~~~~~~~~~~~~~~

* :math:`\mathbf{A}` on the whole model;
* :math:`V` in the solid conductors (``SOLID_CONDUCTOR``), nodal load ``AMPS``;
* ``VOLT`` and ``EMF`` in the stranded coils, nodal load ``CURR`` (:ref:`theory-stranded-coil`).

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Domain
     - :math:`\mathbf{J}` (output ``Jt``)
   * - ``NON_CONDUCTING``, ``SOFT_MAGNETIC``, ``HARD_MAGNETIC``
     - :math:`\mathbf{0}`
   * - ``STRANDED_CONDUCTOR``
     - :math:`\mathbf{J}_s`
   * - ``SOLID_CONDUCTOR``
     - :math:`-\sigma\,(\partial_t\mathbf{A}+\nabla V)`
   * - ``STRANDED_COIL``
     - :math:`n\,\mathbf{T}\,(\mathrm{VOLT}-\mathrm{EMF})=(N/S)\,i\,\mathbf{T}`

In 3D a massive conductor always has the potential :math:`V`: charge conservation is part of the problem, and a
conductor with terminals is a solid conductor where currents are injected (``AMPS``) or potentials fixed (``_V``).
``CIRCUIT_COUPLED_CONDUCTOR``, which in 2D adds the terminal voltage, is therefore not a domain of the 3D element types
(:doc:`domains`).

Weak form
~~~~~~~~~

Find :math:`(\mathbf{A},V)` such that, for every test pair :math:`(\mathbf{w},q)`,

.. math::

   \int_\Omega\nabla\times\mathbf{w}\cdot\nu\,\nabla\times\mathbf{A}
   + g(\mathbf{w},\mathbf{A})
   + \int_{\Omega_c}\mathbf{w}\cdot\sigma\,(\partial_t\mathbf{A}+\nabla V)
   + \text{(coil terms)}
   = \int_{\Omega_s}\mathbf{w}\cdot\mathbf{J}_s + \int_\Omega\nabla\times\mathbf{w}\cdot\mathbf{M},

.. math::

   \int_{\Omega_c}\nabla q\cdot\sigma\,(\partial_t\mathbf{A}+\nabla V) = \sum_k I_k\,q(\mathbf{x}_k),

with :math:`g` the gauge term of the nodal type (below) and :math:`I_k` the loads ``AMPS``. The second equation is the
weak form of :math:`\nabla\cdot\mathbf{J}=0` in the conductors with the currents :math:`I_k` entering at the nodes.

Discretization
~~~~~~~~~~~~~~

**Nodal** (``AV3D_nodal``). Each component of :math:`\mathbf{A}` is a Lagrange field of order ``option2``:
:math:`\mathbf{A}_h=\sum_jN_j\,\mathbf{a}_j`, three DOFs per node. The components are continuous everywhere, including
the normal component at material interfaces, where the exact :math:`\mathbf{A}` may jump; the gauge term
:math:`g(\mathbf{w},\mathbf{A})=\int\nu_g\,(\nabla\cdot\mathbf{w})(\nabla\cdot\mathbf{A})`, with the scalar reluctivity
:math:`\nu_g=\sqrt[3]{|\det\nu|}` of the cell, makes the curl-curl operator definite (Coulomb gauge in penalty form)
[BiroPreis1989]_. Nodal elements accept any first- or second-order mesh, and their Dirichlet values and couplings are
vectors in the nodal frames (``node_angles``).

**Edge** (``AV3D_edge``). :math:`\mathbf{A}_h=\sum_kW_k\,a_k` with the first-order Nedelec shape functions of the first
kind (``NEDELEC_ONE``): one DOF per edge, the circulation :math:`a_k=\int_{\text{edge }k}\mathbf{A}\cdot d\mathbf{l}`.
Only the tangential component is continuous across faces, which is what the physics requires at material interfaces
[Biro1999]_. There is no gauge term: the gradient part of :math:`\mathbf{A}` is removed with a tree-cotree gauge
(:ref:`theory-gauge`). The edge DOFs are stored on the midside nodes of a quadratic mesh: an edge model needs TET10,
HEX20/27, PRISM15/18 or PYRAMID13 cells, ``option2`` being the order of :math:`V`, VOLT and EMF.

The shape functions on a cell are those of the reference cell mapped by the geometric Jacobian :math:`J_e` (the
covariant map :math:`\mathbf{W}=J_e^{-T}\hat{\mathbf{W}}` for edge functions, with
:math:`\nabla\times\mathbf{W}=J_e\hat\nabla\times\hat{\mathbf{W}}/\det J_e`); this geometric mapping is independent of
the rotation of the element coordinate system, which only rotates the coefficients (:ref:`theory-element-cs`).

Element matrices
~~~~~~~~~~~~~~~~

With :math:`\boldsymbol\phi_i` the vector shape functions of :math:`\mathbf{A}` (:math:`N_i\mathbf{e}_c` nodal,
:math:`\mathbf{W}_i` edge) and :math:`N_l` those of :math:`V`:

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Domain
     - Matrix
     - Term
   * - all
     - K
     - :math:`(K_{AA})_{ij}=\int\nabla\times\boldsymbol\phi_i\cdot\nu\nabla\times\boldsymbol\phi_j`
       (+ :math:`\int\nu_g\nabla\cdot\boldsymbol\phi_i\,\nabla\cdot\boldsymbol\phi_j`, nodal)
   * - STRANDED_CONDUCTOR
     - F
     - :math:`\int\boldsymbol\phi_i\cdot\mathbf{J}_s`
   * - all but NON_CONDUCTING, SOFT_MAGNETIC
     - F
     - :math:`\int\nabla\times\boldsymbol\phi_i\cdot\mathbf{M}`
   * - SOLID
     - C
     - :math:`(C_{AA})_{ij}=\int\boldsymbol\phi_i\cdot\sigma\boldsymbol\phi_j`,
       :math:`(C_{VA})_{lj}=\int\nabla N_l\cdot\sigma\boldsymbol\phi_j`
   * - SOLID
     - K
     - :math:`(K_{AV})_{il}=\int\boldsymbol\phi_i\cdot\sigma\nabla N_l`,
       :math:`(K_{VV})_{lm}=\int\nabla N_l\cdot\sigma\nabla N_m`
   * - STRANDED_COIL
     - K, C
     - the coil blocks of :ref:`theory-stranded-coil`

and the element system of a solid conductor

.. math::

   \begin{bmatrix} C_{AA} & 0 \\ C_{VA} & 0 \end{bmatrix}
   \begin{bmatrix} \dot{\mathbf{a}} \\ \dot{\mathbf{V}} \end{bmatrix}
   +
   \begin{bmatrix} K_{AA} & K_{AV} \\ 0 & K_{VV} \end{bmatrix}
   \begin{bmatrix} \mathbf{a} \\ \mathbf{V} \end{bmatrix}
   =
   \begin{bmatrix} \mathbf{F}_A \\ \mathbf{I} \end{bmatrix},
   \qquad C_{VA}=K_{AV}^T.

The pair (K, C) is not symmetric: the coupling of :math:`\mathbf{A}` with :math:`V` is in K in the A rows and in C in
the V rows. The time-integrated types make it symmetric (:ref:`theory-timint`).

Notes
~~~~~

* On a **static** step (no C) the V rows reduce to :math:`K_{VV}\mathbf{V}=\mathbf{I}` (a conduction problem with
  the injected currents) and the A rows to magnetostatics with the conduction current :math:`-\sigma\nabla V`.
* Each conductor needs :math:`V` fixed somewhere (``_V``), or a symmetry plane with :math:`V=0`; otherwise the V block
  is singular.
* Cyclic symmetry couples A and V between the faces of the sector, in rotated nodal frames for the nodal type;
  edge DOFs are coupled through their midside nodes with the orientation of the edges
  (:doc:`../user_manual/couplings`).
* ``Js`` is a load of the cell, given in global components or in the element coordinate system of the cell; the
  output ``Jt`` evaluates it at every output, so it follows the current loads.
