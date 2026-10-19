.. _theory-constitutive:

Materials
---------

The material of a cell (its ``MATERIAL`` row) gives the reluctivity, the conductivity and the magnetization; a
nonlinear law can replace the reluctivity by a B-H curve. The parameters are tensors and vectors in the element
coordinate system of the cell.

Magnetic law
~~~~~~~~~~~~

pyTEC3D uses the law of a linear magnetic material with a magnetization,

.. math::

   \mathbf{B} = \mu_0\,\mu_r\,(\mathbf{H}+\mathbf{M})
   \qquad\Longleftrightarrow\qquad
   \mathbf{H} = \nu\,\mathbf{B} - \mathbf{M},
   \qquad \nu = (\mu_0\,\mu_r)^{-1},

with :math:`\mu_r` the relative permeability tensor (``mat.mur``, 3x3) and :math:`\mathbf{M}` the magnetization
[A/m] (``mat.M``, default 0). In the weak form the magnetization is a source,
:math:`\int\nabla\times\mathbf{w}\cdot\mathbf{M}` (a magnetized region is equivalent to the magnetization currents
:math:`\nabla\times\mathbf{M}` in its volume and :math:`\mathbf{M}\times\mathbf{n}` on its surface).

* **Where.** The magnetization is a parameter of every domain except ``NON_CONDUCTING`` and ``SOFT_MAGNETIC``:
  ``HARD_MAGNETIC`` (permanent magnets), and the conducting domains, so that a magnetized region can also carry
  currents. In ``NON_CONDUCTING`` and ``SOFT_MAGNETIC`` cells it is not read.
* **Permanent magnets.** For a magnet with remanence :math:`B_r` and recoil permeability :math:`\mu_\text{rec}` along
  its easy axis :math:`\mathbf{e}`, :math:`\mu_r=\mu_\text{rec}` and :math:`\mathbf{M}=B_r/(\mu_0\mu_\text{rec})\,
  \mathbf{e}`, i.e. the coercive field. The easy axis is conveniently given in an element coordinate system.
* **Output.** The field ``H`` of the output is :math:`\nu\,\mathbf{B}`: in a magnetized cell :math:`\mathbf{M}` is
  not subtracted.

Conductivity
~~~~~~~~~~~~

In the massive conductors (``SOLID_CONDUCTOR``, ``CIRCUIT_COUPLED_CONDUCTOR``) :math:`\mathbf{J}=\sigma\mathbf{E}` with
:math:`\sigma=\rho^{-1}`, :math:`\rho` the resistivity tensor (``mat.rho``, 3x3, Ohm m). An anisotropic resistivity
models laminated or layered conductors (e.g. a much larger resistivity across the laminations). The other domains do
not read ``rho``.

B-H curves
~~~~~~~~~~

A saturating material is described by the law ``bh_curve`` (``nonlinear_laws.json``), selected for a MATERIAL id in
the section ``nonlinear_laws`` of the model data (:doc:`../user_manual/nonlinear`). It replaces the isotropic
reluctivity by a function of :math:`B=|\mathbf{B}|`, evaluated at every quadrature point:

.. math::

   \mathbf{H} = h(B)\,\hat{\mathbf{b}} = \nu_s(B)\,\mathbf{B},
   \qquad
   \nu_s(B) = \frac{h(B)}{B},
   \qquad
   \hat{\mathbf{b}} = \frac{\mathbf{B}}{B},

with :math:`h` the curve :math:`H(B)` of the material. The derivative of :math:`\mathbf{H}` with respect to
:math:`\mathbf{B}`, needed by Newton's method, is the **differential reluctivity**

.. math::

   \nu_t(\mathbf{B}) = \frac{\partial\mathbf{H}}{\partial\mathbf{B}}
   = \nu_s(B)\,\big(\mathbf{I}-\hat{\mathbf{b}}\otimes\hat{\mathbf{b}}\big) + h'(B)\,\hat{\mathbf{b}}\otimes\hat{\mathbf{b}},

symmetric and transversely isotropic about the field: the secant :math:`\nu_s` governs a rotation of :math:`\mathbf{B}`
at constant magnitude, the slope :math:`h'` a change of magnitude. In saturation :math:`h'\to\nu_0` while
:math:`\nu_s` is still small, so :math:`\nu_t` is strongly anisotropic, which slows Newton's method down
(:doc:`nonlinear`).

**What the law does with the table.** The table of the model data (rows :math:`[H, B]`, :math:`B` strictly
increasing) is used as given: :math:`h(B)` is the piecewise-linear interpolant of the points, extended beyond the
first and the last point along the first and the last segment; :math:`h'(B)` is the slope of the segment;
:math:`\nu_s=h(B)/B` for :math:`B>10^{-14}` T and :math:`h'(0)` at :math:`B=0`. Two consequences:

* the behaviour in deep saturation is the slope of the last segment: end the table with a point far in saturation on
  the slope wanted there (close to :math:`\mu_0`, i.e. :math:`h'\approx1/\mu_0`);
* a monotone table (:math:`H` increasing with :math:`B`) keeps :math:`\nu_t` positive definite and the Newton
  iteration well posed.

The solver reports at every step the largest :math:`|B|` reached in each nonlinear material, with a warning when it is
beyond the last point of the table.

.. _theory-element-cs:

Element coordinate systems
~~~~~~~~~~~~~~~~~~~~~~~~~~

The vector and tensor parameters of a cell (``mur``, ``rho``, ``M``, the coil direction ``T``, the load ``Js``) are
given in its element coordinate system, the entry of ``coordinate_systems`` named by its ``COORDINATE_SYSTEM`` cell
value (0: the global axes). With :math:`\mathsf{R}` the rotation from the local to the global axes at the centre of
the cell, the element type uses

.. math::

   \mathbf{p}_\text{glob} = \mathsf{R}\,\mathbf{p}_\text{loc},
   \qquad
   P_\text{glob} = \mathsf{R}\,P_\text{loc}\,\mathsf{R}^T,

once per cell: the shape functions stay global and only the coefficients are rotated. A cylindrical system makes a
radial or toroidal direction constant in the local components (e.g. ``T = (0, 1, 0)`` is the toroidal direction of
every cell of a 3D coil, :doc:`../user_manual/coordinate_systems`). The rotations are computed when the model is
loaded into a solver session.
