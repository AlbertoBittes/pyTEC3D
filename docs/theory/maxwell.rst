The eddy-current problem
------------------------

Magneto-quasi-static approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pyTEC3D solves the low-frequency (eddy-current) regime of Maxwell's equations: the fields change slowly enough that
the displacement current :math:`\partial_t\mathbf{D}` is negligible against the conduction current, i.e. the size of
the model is small against the wavelength and capacitive effects are negligible [Bossavit1998]_ [Meunier2008]_. The
equations are then

.. math::

   \nabla\times\mathbf{H} = \mathbf{J},
   \qquad
   \nabla\cdot\mathbf{B} = 0,
   \qquad
   \nabla\times\mathbf{E} = -\partial_t\mathbf{B},
   \qquad
   \nabla\cdot\mathbf{J} = 0,

with the constitutive laws of :doc:`constitutive`, :math:`\mathbf{H}=\mathbf{H}(\mathbf{B})` and, in the massive
conductors, Ohm's law :math:`\mathbf{J}=\sigma\mathbf{E}`. Charge conservation :math:`\nabla\cdot\mathbf{J}=0` follows
from Ampere's law; it is an independent equation only where the electric potential is an unknown.

Potentials
~~~~~~~~~~

Since :math:`\nabla\cdot\mathbf{B}=0`, :math:`\mathbf{B}=\nabla\times\mathbf{A}`; Faraday's law then gives
:math:`\mathbf{E}=-\partial_t\mathbf{A}-\nabla V`. The current density of each domain (:doc:`domains`) is

.. math::

   \mathbf{J} =
   \begin{cases}
   \mathbf{0} & \Omega_0 \ (\text{non-conducting, soft and hard magnetic}),\\
   \mathbf{J}_s & \Omega_s \ (\text{stranded conductor: imposed}),\\
   -\sigma\,(\partial_t\mathbf{A} + \nabla V) & \Omega_c \ (\text{solid and circuit-coupled conductors}),\\
   \dfrac{N}{S}\,i\,\mathbf{T} & \Omega_\text{coil} \ (\text{stranded coil}),
   \end{cases}

and Ampere's law becomes the equation of :math:`\mathbf{A}`,

.. math::

   \nabla\times\big(\nu\,\nabla\times\mathbf{A} - \mathbf{M}\big) = \mathbf{J},

completed in :math:`\Omega_c` by charge conservation,
:math:`\nabla\cdot\big(\sigma(\partial_t\mathbf{A}+\nabla V)\big)=0`, and in each coil by its circuit equation. The
term :math:`-\partial_t\mathbf{A}` is the electric field induced by the change of flux (eddy currents),
:math:`-\nabla V` the field due to charges and terminals (applied voltages, injected currents, the redistribution that
keeps the currents inside the conductors).

Weak form
~~~~~~~~~

Multiplying by a test function :math:`\mathbf{w}` and integrating by parts over :math:`\Omega` gives the equation the
element types discretize:

.. math::

   \int_\Omega \nabla\times\mathbf{w}\cdot\nu\,\nabla\times\mathbf{A}
   + \int_{\Omega_c}\mathbf{w}\cdot\sigma\,(\partial_t\mathbf{A}+\nabla V)
   - \int_{\Omega_\text{coil}}\mathbf{w}\cdot\frac{N}{S}\,i\,\mathbf{T}
   = \int_{\Omega_s}\mathbf{w}\cdot\mathbf{J}_s + \int_\Omega\nabla\times\mathbf{w}\cdot\mathbf{M}
   + \oint_{\partial\Omega}\mathbf{w}\cdot(\mathbf{n}\times\mathbf{H}),

with, in the conductors, the weak charge conservation (test function :math:`q`)

.. math::

   \int_{\Omega_c}\nabla q\cdot\sigma\,(\partial_t\mathbf{A}+\nabla V) = \sum_k I_k\,q(\mathbf{x}_k),

where :math:`I_k` are currents injected at nodes (the nodal load ``AMPS``): the boundary term of the integration by
parts, :math:`-\oint q\,\mathbf{J}\cdot\mathbf{n}`, is zero on the insulated surface of a conductor and is the current
entering through a terminal. Each formulation (:doc:`formulation_2d`, :doc:`formulation_3d`) reduces these integrals
to the fields it has.

Boundary and interface conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weak form imposes the continuity of the tangential :math:`\mathbf{H}` and of the normal :math:`\mathbf{B}` between
cells. On the external boundary:

* **Dirichlet** :math:`\mathbf{n}\times\mathbf{A}=0` (``_A`` = 0 on the tangential components): :math:`\mathbf{B}\cdot
  \mathbf{n}=0`, the flux is parallel to the boundary. Used on far boundaries (the field of the sources decays) and on
  symmetry planes to which the currents are normal.
* **Natural** (no condition): :math:`\mathbf{n}\times\mathbf{H}=0`, the flux is normal to the boundary. A symmetry
  plane parallel to the currents.
* **Axis** of an axisymmetric model or of a sector: :math:`A_\theta=0` (2D), the components of :math:`\mathbf{A}`
  normal to the axis vanish (3D).
* **Cyclic symmetry** between the faces of a sector: the values on one face equal those on the other, in rotated
  frames (:doc:`../user_manual/couplings`).

The scalar potential :math:`V` of a conductor is defined up to a constant: it must be fixed on one node of every
conductor (``_V``), unless the conductor touches a symmetry plane where :math:`V` is fixed, or has a terminal at a
fixed potential.

.. _theory-gauge:

Gauge
~~~~~

The vector potential is defined up to a gradient: :math:`\mathbf{A}+\nabla\phi` gives the same :math:`\mathbf{B}`. In
the conductors the term :math:`\sigma\partial_t\mathbf{A}` fixes this freedom (together with :math:`V`); in the regions
without conductivity the discrete curl-curl operator is singular and a gauge is needed:

* the **nodal** 3D element types add the Coulomb gauge term
  :math:`\int\nu_g\,(\nabla\cdot\mathbf{w})(\nabla\cdot\mathbf{A})` with :math:`\nu_g=\sqrt[3]{|\det\nu|}`, which makes
  the operator definite (penalty form of :math:`\nabla\cdot\mathbf{A}=0`; the Lagrange interpolation imposes the
  continuity of all the components of :math:`\mathbf{A}`) [BiroPreis1989]_;
* the **edge** 3D element types use a **tree-cotree gauge**: the circulations of :math:`\mathbf{A}` on the edges of a
  spanning tree of the mesh are fixed to zero (they are the gradient part), the cotree edges carry the solution
  [AlbaneseRubinacci1988]_ [MangesCendes1995]_ (:func:`pytec3d.tree_gauge.tree_gauge`);
* in **2D** the out-of-plane potential :math:`A_z` has no gradient freedom: no gauge is needed.
