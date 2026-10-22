.. _theory-nonlinear-elements:

Nonlinear materials and Newton-Raphson
--------------------------------------

With a saturating material (the law ``bh_curve``, :doc:`constitutive`) the reluctivity depends on the solution and
each time step is a nonlinear problem, solved with Newton's method. This page gives the residual and the Jacobian of
the A-V formulations in block form; the conduction and circuit terms stay linear.

Residual
~~~~~~~~

Write the step of the time scheme (:ref:`theory-time-integration`) as
:math:`\dot{\mathbf{u}}^{n+1}=c_v\,\mathbf{u}^{n+1}+\mathbf{h}_v`, with :math:`c_v=1/(\theta\Delta t)` and
:math:`\mathbf{h}_v` known from the start of the step. The residual of the step is

.. math::

   \mathbf{R}(\mathbf{u}) = \mathbf{g}(\mathbf{u}) + (\mathbf{K}_\text{lin} + c_v\mathbf{C})\,\mathbf{u}
   - \big(\mathbf{F}(t^{n+1}) - \mathbf{C}\,\mathbf{h}_v\big),
   \qquad
   g_i(\mathbf{u}) = \int_{\Omega_\text{nl}}\nabla\times\boldsymbol\phi_i\cdot\nu_s(B)\,\mathbf{B}_h,

where :math:`\mathbf{g}` collects the curl-curl terms of the cells with a law (:math:`\Omega_\text{nl}`) and
:math:`\mathbf{K}_\text{lin}` every other stiffness term: the curl-curl term of the linear cells, the gauge term, the
conduction and coil terms. The magnetization and the source currents are in :math:`\mathbf{F}`: they do not depend on
the solution.

Jacobian
~~~~~~~~

.. math::

   \mathbf{J}(\mathbf{u}) = \frac{\partial\mathbf{R}}{\partial\mathbf{u}}
   = \mathbf{K}_\nu(\mathbf{u}) + \mathbf{K}_\text{lin} + c_v\,\mathbf{C},
   \qquad
   (K_\nu)_{ij} = \int_{\Omega_\text{nl}}\nabla\times\boldsymbol\phi_i\cdot\nu_t(\mathbf{B}_h)\,\nabla\times\boldsymbol\phi_j,

with the differential reluctivity :math:`\nu_t` of :doc:`constitutive`. In the blocks of the fields
:math:`(\mathbf{A},V,\mathrm{VOLT},\mathrm{EMF})` of a model with conductors and coils:

.. math::

   \mathbf{J} =
   \begin{bmatrix}
   K_\nu(\mathbf{B}) + K_{AA}^\text{lin} + c_v\, C_{AA} & K_{AV} & K_{AU} & -K_{AU} \\
   c_v\, C_{VA} & K_{VV} & 0 & 0 \\
   c_v\,s\,K_{AU}^T & 0 & s\,K_{UU} & 0 \\
   -c_v\, K_{AU}^T & 0 & 0 & -K_{UU}
   \end{bmatrix},

where only the first block depends on the state. With the time-integrated unknowns (``_timint``) the coupling blocks
move to the :math:`c_v` terms and the matrix is symmetric. In the nodal 3D types the gauge coefficient :math:`\nu_g` is
also replaced by the law (its secant value), and its term enters :math:`K_\nu` with the secant only.

Iteration
~~~~~~~~~

.. math::

   \mathbf{J}(\mathbf{u}^k)\,\delta\mathbf{u} = \mathbf{R}(\mathbf{u}^k),
   \qquad
   \mathbf{u}^{k+1} = \mathbf{u}^k - \lambda\,\delta\mathbf{u},

starting from the solution of the previous step. :math:`\lambda\in[\lambda_\text{min},1]` is chosen by a line search
when the full step does not reduce the residual [DennisSchnabel1996]_. The linear system of each iteration is solved
with the linear solver
settings; the iterations stop when the residual or the step is small enough
(``solver_settings["nonlinear"]``, :doc:`../user_manual/solver_settings`). The Dirichlet values and the couplings are
applied to every iterate.

On a static step :math:`c_v=0`: the conduction terms disappear and Newton solves nonlinear magnetostatics.

.. _theory-newton:

What the law contributes, in general
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The element types are written as kernels (:ref:`theory-kernels`), and a law can replace any of their coefficients.
For a kernel :math:`\int\mathcal{T}(\boldsymbol\phi_i)^T P\,\mathcal{L}(\partial_t^d\mathbf{u})` whose coefficient
:math:`P` is replaced by a law with input :math:`\mathbf{x}` (for ``bh_curve``, :math:`\mathbf{x}=\mathbf{B}`), the law
gives a secant :math:`S(\mathbf{x})` and, optionally, a tangent :math:`\mathbf{T}(\mathbf{x})=\partial(S\mathbf{x})
/\partial\mathbf{x}`, and

.. math::

   \mathbf{R} \mathrel{+}= \int\mathcal{T}(\boldsymbol\phi)^T S(\mathbf{x})\,\mathcal{L}(\partial_t^d\mathbf{u}),
   \qquad
   \mathbf{J} \mathrel{+}= c_d\int\mathcal{T}(\boldsymbol\phi)^T
   \begin{cases}\mathbf{T}(\mathbf{x}) & \mathcal{L}(\mathbf{u})=\mathbf{x},\ d=0\\ S(\mathbf{x}) & \text{otherwise}
   \end{cases}\;\mathcal{L},

with :math:`c_0=1`, :math:`c_1=c_v`, :math:`c_2=c_a`. When the law coefficient multiplies a trial that is not its own
input (a field-dependent conductivity on :math:`\partial_t\mathbf{A}`, the gauge term), the Jacobian uses the secant
only: Newton still converges to the solution, more slowly.

Convergence in practice
~~~~~~~~~~~~~~~~~~~~~~~

* **Saturation.** Deep in saturation the slope :math:`h'` tends to :math:`1/\mu_0` while the secant :math:`\nu_s` is
  still much smaller: the tangent is very anisotropic (the secant governs rotations of the field, the slope changes of
  its magnitude). Large load increments into saturation are the usual cause of slow convergence: use smaller time
  steps or ramp the load.
* **The table.** A monotone B-H table keeps :math:`\nu_t` positive definite; the slope of its last segment sets the
  behaviour beyond it. The report of each step gives the largest :math:`|B|` reached in every nonlinear material.
* **Cost.** Only the cells with a law are re-integrated at each iteration; the linear terms are assembled once.
* **Permanent magnets.** The magnetization is a constant source: it does not enter the Jacobian. The model has a
  fixed recoil line and does not represent irreversible demagnetization; check the operating point of the magnets
  against their demagnetization curve in post-processing.
