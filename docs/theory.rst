Theory
======

Electromagnetic field domains
-----------------------------

Overview
~~~~~~~~

Magnetic fields may exist as a result of an electric current, a permanent magnet, or an applied
external field. In electromagnetic analysis, components are typically categorized by their
electric and magnetic properties. These categories define which physics are active in each
subdomain and which constitutive laws apply.

.. figure:: _static/field_domains.svg
   :alt: Conceptual electromagnetic field domains: background air, conductors, and magnetic media.
   :align: center
   :width: 90%

   Conceptual domain types used to organize electric and magnetic material behavior.

Electric characteristics
~~~~~~~~~~~~~~~~~~~~~~~~

Possible electric characteristics include:

* **No current** in the body (insulating region)
* **Solid conductor**, current known or imposed, no eddy-current effects (static or circuit-fed)
* **Stranded conductor**, current density known, no eddy-current effects (e.g., source coils)
* **Solid conductor with eddy currents** (harmonic or transient analyses only)

Magnetic characteristics
~~~~~~~~~~~~~~~~~~~~~~~~

Magnetic characteristics include:

* **Non-magnetic** (typically air, copper, aluminum)
* **Soft magnetic** (typically iron or steel)
* **Hard magnetic** (permanent magnets such as samarium cobalt or Alnico)

The analysis model combines these electric and magnetic properties to define the material
domains and their governing equations.

Maxwell equations for the eddy-current regime
---------------------------------------------

Scope and assumptions
~~~~~~~~~~~~~~~~~~~~~

The **eddy-current regime** refers to electromagnetic phenomena where **time-varying magnetic
fields induce currents** in conducting media, while **wave propagation effects are negligible**.
This is typically modeled with the **magneto-quasi-static (MQS)** approximation, in which the
displacement current is neglected in Ampere's law.

Let:

* **E** be the electric field (:math:`\mathbf{E}`) ``[V/m]``
* **H** be the magnetic field (:math:`\mathbf{H}`) ``[A/m]``
* **B** be the magnetic flux density (:math:`\mathbf{B}`) ``[T]``
* **D** be the electric flux density (:math:`\mathbf{D}`) ``[C/m^2]``
* **rho** be the free charge density (:math:`\rho`) ``[C/m^3]``
* **J** be the (conduction) current density (:math:`\mathbf{J}`) ``[A/m^2]``
* **J_s** be an externally imposed source current density (:math:`\mathbf{J}_s`) ``[A/m^2]``

Material laws (generic form):

* :math:`\mathbf{B} = \mu \mathbf{H}`
* :math:`\mathbf{D} = \epsilon \mathbf{E}`

Governing Maxwell equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In differential form, Maxwell's equations are:

1. **Gauss's law (electric):**

   .. math::

      \nabla \cdot \mathbf{D} = \rho

2. **Gauss's law (magnetic):**

   .. math::

      \nabla \cdot \mathbf{B} = 0

3. **Faraday's law (electromagnetic induction):**

   .. math::

      \nabla \times \mathbf{E} = - \frac{\partial \mathbf{B}}{\partial t}

4. **Ampere-Maxwell law (full form):**

   .. math::

      \nabla \times \mathbf{H} = \mathbf{J} + \frac{\partial \mathbf{D}}{\partial t}

Eddy-current (MQS) approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the eddy-current regime, the displacement current term is neglected:

.. math::

   \frac{\partial \mathbf{D}}{\partial t} \approx \mathbf{0}

so Ampere's law reduces to:

.. math::

   \nabla \times \mathbf{H} = \mathbf{J}

(or, when explicitly separating impressed currents, see below).

This approximation is valid when the characteristic length scales are small compared to the
electromagnetic wavelength and when capacitive effects are negligible relative to conductive
effects.

Current density modeling: imposed vs. evaluated currents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Imposed (source) current density (:math:`\mathbf{J}_s`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, part (or all) of the current is **prescribed** by the user or excitation
model, e.g. coils driven by a current source. This is represented as an **impressed current
density**:

.. math::

   \mathbf{J}_s

It is treated as a known input to the formulation.

Evaluated current density (:math:`\mathbf{J}`) in conducting media
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In conducting regions, the current density may be **computed** from the electric field using a
constitutive law. A common decomposition is:

.. math::

   \mathbf{J} = \sigma \left(\underbrace{\mathbf{E}_i}_{\text{induction term}}
   + \underbrace{\mathbf{E}_c}_{\text{conduction term}}\right)
   = \underbrace{\mathbf{J}_i}_{\text{induced}}
   + \underbrace{\mathbf{J}_c}_{\text{conduction}}

where:

* :math:`\sigma` is the electrical conductivity ``[S/m]``
* :math:`\mathbf{E}_i` is the **induction-related** electric field contribution (associated with
  time-varying magnetic flux)
* :math:`\mathbf{E}_c` is the **conduction (potential-driven)** electric field contribution

Correspondingly:

.. math::

   \mathbf{J}_i = \sigma \mathbf{E}_i,
   \qquad
   \mathbf{J}_c = \sigma \mathbf{E}_c

**Interpretation**

* The **induction component** (:math:`\mathbf{J}_i`) captures eddy currents generated by
  Faraday's law due to :math:`\partial \mathbf{B}/\partial t`.
* The **conduction component** (:math:`\mathbf{J}_c`) captures currents driven by an electric
  potential difference (e.g., electrodes, imposed voltages).

Combined statement used in formulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the modeling choice, the current entering Ampere's law may be:

* **Imposed directly** as :math:`\mathbf{J}_s`, and/or
* **Evaluated** as :math:`\mathbf{J} = \sigma(\mathbf{E}_i + \mathbf{E}_c)
  = \mathbf{J}_i + \mathbf{J}_c`

In implementations that separate impressed and conductive currents, a common form is:

.. math::

   \nabla \times \mathbf{H} = \mathbf{J}_s + \mathbf{J}

with :math:`\mathbf{J}` computed from the conductivity model above.

PLANE233 2D FE formulation (A-only, A-V, time-integrated V)
-----------------------------------------------------------

This section documents a 2D finite-element formulation aligned with the PLANE233 element
capabilities for planar and axisymmetric cases, using local notation and the PLANE233 DOFs:
``AZ`` plus optional ``VOLT`` (electric potential) or time-integrated ``VOLT``. [1]_

Geometry reduction and the THK (out-of-plane measure)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let the 2D computational domain be :math:`\Omega \subset \mathbb{R}^2`.

Planar (2D "plane")
^^^^^^^^^^^^^^^^^^^

A 3D prismatic body is represented by extruding :math:`\Omega` in the out-of-plane direction
with thickness

.. math::

   \ell \equiv t = \mathrm{THK}.

All volume integrals become

.. math::

   \int_{\Omega^{3D}} (\cdot)\, dV
   = \int_{\Omega} (\cdot)\, t\, d\Omega.

PLANE233 explicitly supports THK as thickness for the plane option. [1]_

Axisymmetric (2D r-z)
^^^^^^^^^^^^^^^^^^^^^

A 3D body is formed by revolving the meridian :math:`\Omega \subset \{(r,z): r \ge 0\}`.
PLANE233 allows modeling a fraction :math:`f = \mathrm{THK}` of the full :math:`360^\circ`. [1]_
Hence the azimuthal span is :math:`\Delta \theta = 2\pi f`, and

.. math::

   \int_{\Omega^{3D}} (\cdot)\, dV
   = \int_{\Omega} (\cdot)\, (2\pi f\, r)\, d\Omega.

So, throughout:

* **planar weight:** :math:`w_{\text{geom}} = t`
* **axisym weight:** :math:`w_{\text{geom}} = 2\pi f\, r`

Unknowns and fields
~~~~~~~~~~~~~~~~~~~

Magnetic potential (2D reduction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PLANE233 uses a single scalar DOF representing the only nonzero component of
:math:`\mathbf{A}`:

* **planar:** :math:`\mathbf{A} = A(x,y,t)\, \mathbf{e}_z`
* **axisym:** :math:`\mathbf{A} = A(r,z,t)\, \mathbf{e}_\theta`

Flux density (:math:`\mathbf{B}=\nabla\times\mathbf{A}`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Planar (in-plane x-y)**

.. math::

   \mathbf{B} =
   \begin{bmatrix}
   \partial_y A \\
   -\partial_x A \\
   0
   \end{bmatrix}
   \equiv \mathbf{C}\, \nabla A,
   \qquad
   \mathbf{C} =
   \begin{bmatrix}
   0 & 1 \\
   -1 & 0
   \end{bmatrix}.

**Axisymmetric (meridian r-z)**

.. math::

   B_r = -\partial_z A,
   \qquad
   B_z = \frac{1}{r}\partial_r(rA) = \partial_r A + \frac{A}{r}.

Constitutive law
^^^^^^^^^^^^^^^^

Let :math:`\boldsymbol{\nu}` be the (possibly orthotropic in-plane) reluctivity tensor in the
2D plane:

.. math::

   \mathbf{H} = \boldsymbol{\nu}\, \mathbf{B}.

Electric field and currents (for A-V cases)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The out-of-plane electric field component driving conduction is modeled via a **lumped
potential drop** across the out-of-plane length :math:`\ell`:

* planar: :math:`\ell = t`
* axisym: :math:`\ell(r) = r, \; \Delta\theta = 2\pi f\, r`

Define the **electric potential drop DOF** :math:`\varphi` (units V). Then

.. math::

   E_{\parallel} = -\partial_t A - \frac{\varphi}{\ell}.

(":math:`\parallel`" means along :math:`z` for planar or along :math:`\theta` for axisym.)

Let :math:`\sigma` be the conductivity in that direction (PLANE233 uses out-of-plane
conductivity for this 2D formulation). The **total current density** in a solid conductor is

.. math::

   J = J_s + \sigma E_{\parallel}
   = J_s - \sigma\left(\partial_t A + \frac{\varphi}{\ell}\right).

Per PLANE233 documentation: in 2D electromagnetic regions, the VOLT DOF represents a
single voltage drop across the out-of-plane measure, enforced by coupling. [1]_

For the **time-integrated** potential option, define

.. math::

   \eta(t) = \int_0^t \varphi(\tau)\, d\tau,
   \qquad
   \varphi = \partial_t \eta,

which is exactly the documented "time-integrated electric potential" meaning of VOLT. [2]_

Governing strong form (magneto-quasistatic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`A` be the scalar magnetic potential DOF. The reduced Ampere law becomes, in both
planar and axisym,

.. math::

   \mathcal{L}(A) = J,

where :math:`\mathcal{L}(\cdot)` is the standard 2D vector-potential operator.

Planar
^^^^^^

.. math::

   -\nabla\cdot\left(\boldsymbol{\nu}\, \nabla A\right) = J

Axisymmetric
^^^^^^^^^^^^

Using :math:`B_r,B_z` above, one may write

.. math::

   -\frac{1}{r}\partial_r\left(r\, \nu\, \partial_r A\right)
   - \partial_z\left(\nu\, \partial_z A\right)
   + \frac{\nu}{r^2} A
   = J

for isotropic :math:`\nu`; the FE form below handles the operator robustly.

Weak forms
~~~~~~~~~~

Let :math:`w` be the test function for :math:`A`. For A-V, let :math:`q` be the test function
for the nodal potential DOF.

A-only (stranded / no eddy currents)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is PLANE233 with only :math:`A` DOF; in time-varying analyses it behaves as a stranded
conductor (eddy currents not modeled). [2]_ Thus :math:`J \equiv J_s` is **imposed**.

Find :math:`A` such that for all :math:`w`:

.. math::

   \int_\Omega (\nabla w)^T \boldsymbol{\nu} (\nabla A)\, w_{\text{geom}}\, d\Omega
   =
   \int_\Omega w\, J_s\, w_{\text{geom}}\, d\Omega
   + \text{(natural boundary terms)}.

Axisymmetric uses the same pattern but with the axisymmetric :math:`\mathbf{B}` operator
(next section shows the discrete operator).

A-V with electric potential (:math:`\varphi`) (non time-integrated VOLT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here :math:`J = J_s - \sigma(\partial_t A + \varphi/\ell)`.

**A-equation (weak Ampere law)**: find :math:`(A,\varphi)` such that for all :math:`w`,

.. math::

   \int_\Omega (\nabla w)^T \boldsymbol{\nu} (\nabla A)\, w_{\text{geom}}\, d\Omega
   + \int_\Omega w\, \sigma\, \partial_t A\, w_{\text{geom}}\, d\Omega
   + \int_\Omega w\, \sigma\, \frac{\varphi}{\ell}\, w_{\text{geom}}\, d\Omega
   =
   \int_\Omega w\, J_s\, w_{\text{geom}}\, d\Omega.

**V-equation (global current/terminal constraint, Galerkin form)**: PLANE233 allows
prescribing total current or potential in massive conductors (VOLT via constraint, AMPS via
nodal load). [1]_ A consistent weak statement is obtained by enforcing that
:math:`\varphi` is work-conjugate to current; in Galerkin form:

.. math::

   \int_\Omega q\, \sigma\, \partial_t A\, \frac{w_{\text{geom}}}{\ell}\, d\Omega
   + \int_\Omega q\, \sigma\, \frac{\varphi}{\ell^2}\, w_{\text{geom}}\, d\Omega
   =
   \langle q, I^{\text{ext}} \rangle,

where :math:`\langle q, I^{\text{ext}} \rangle` represents the external current excitation
(e.g., distributed to coupled VOLT DOFs so that the assembled equation corresponds to the
imposed total current). If instead :math:`\varphi` is prescribed (voltage-driven), this
equation yields the reaction current.

Practical note (matches PLANE233 behavior): in 2D electromagnetic regions all nodal VOLT are
coupled so :math:`\varphi` is constant over the conductor cross-section, i.e. one effective
DOF. [1]_

A-V with time-integrated potential (:math:`\eta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :math:`\varphi = \partial_t \eta`. Then
:math:`J = J_s - \sigma(\partial_t A + \partial_t\eta/\ell)`.

Find :math:`(A,\eta)` such that for all :math:`(w,q)`:

.. math::

   \int_\Omega (\nabla w)^T \boldsymbol{\nu} (\nabla A)\, w_{\text{geom}}\, d\Omega
   + \int_\Omega w\, \sigma\, \partial_t A\, w_{\text{geom}}\, d\Omega
   + \int_\Omega w\, \sigma\, \frac{\partial_t \eta}{\ell}\, w_{\text{geom}}\, d\Omega
   =
   \int_\Omega w\, J_s\, w_{\text{geom}}\, d\Omega.

.. math::

   \int_\Omega q\, \sigma\, \partial_t A\, \frac{w_{\text{geom}}}{\ell}\, d\Omega
   + \int_\Omega q\, \sigma\, \frac{\partial_t \eta}{\ell^2}\, w_{\text{geom}}\, d\Omega
   =
   \langle q, I^{\text{ext}} \rangle.

This is the **symmetric strong-coupling** time-integrated option documented for PLANE233. [1]_

FE discretization (shape functions + Gauss quadrature)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let the element have :math:`n` nodes and scalar shape functions :math:`N_i(\xi,\eta)`. Collect

.. math::

   \mathbf{N}(\xi,\eta) = [N_1,\dots,N_n],
   \quad
   A_h = \mathbf{N}\, \mathbf{a}_e,
   \quad
   \varphi_h = \mathbf{N}\, \boldsymbol{\varphi}_e,
   \quad
   \eta_h = \mathbf{N}\, \boldsymbol{\eta}_e.

At Gauss point :math:`g`: weight :math:`\omega_g`, Jacobian determinant :math:`|J_g|`, and
geometry weight :math:`w_{\text{geom}}(g)` (planar :math:`t`, axisym :math:`2\pi f r_g`).

Magnetic "curl" operator matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Planar
^^^^^^

Define the in-plane gradients :math:`\nabla N_i=[\partial_x N_i,\partial_y N_i]^T`.
Define the :math:`2\times n` operator that maps nodal :math:`\mathbf{a}_e` to
:math:`(B_x,B_y)`:

.. math::

   \mathbf{B}_A(g) =
   \begin{bmatrix}
   \partial_y \mathbf{N}(g) \\
   -\partial_x \mathbf{N}(g)
   \end{bmatrix},
   \quad
   \mathbf{b}(g) = \mathbf{B}_A(g)\, \mathbf{a}_e.

Axisymmetric
^^^^^^^^^^^^

Using :math:`(r,z)` and :math:`A = A_\theta`, define the :math:`2\times n` operator mapping
:math:`\mathbf{a}_e` to :math:`(B_r,B_z)`:

.. math::

   \mathbf{B}_A^{\text{axi}}(g) =
   \begin{bmatrix}
   -\partial_z \mathbf{N}(g) \\
   \partial_r \mathbf{N}(g) + \dfrac{1}{r_g}\mathbf{N}(g)
   \end{bmatrix},
   \quad
   \mathbf{b}(g) = \mathbf{B}_A^{\text{axi}}(g)\, \mathbf{a}_e.

Element matrices (A-only and A-V)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below, :math:`\nu_g` denotes reluctivity (scalar or :math:`2\times 2`) evaluated at Gauss
point :math:`g`. :math:`\sigma_g` is the conductivity in the out-of-plane current direction.
The out-of-plane length is :math:`\ell_g` (planar :math:`\ell_g=t`; axisym
:math:`\ell_g = 2\pi f r_g`).

A-only formulation (element stiffness and source vector)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Planar**

.. math::

   \mathbf{K}_{AA}^e =
   \sum_{g}
   \omega_g |J_g|\, t\,
   \mathbf{B}_A(g)^T\, \nu_g\, \mathbf{B}_A(g)

.. math::

   \mathbf{f}_A^e =
   \sum_{g}
   \omega_g |J_g|\, t\,
   \mathbf{N}(g)^T\, J_s(g)

**Axisymmetric**

.. math::

   \mathbf{K}_{AA}^e =
   \sum_{g}
   \omega_g |J_g|\, (2\pi f r_g)\,
   \mathbf{B}_A^{\text{axi}}(g)^T\, \nu_g\, \mathbf{B}_A^{\text{axi}}(g)

.. math::

   \mathbf{f}_A^e =
   \sum_{g}
   \omega_g |J_g|\, (2\pi f r_g)\,
   \mathbf{N}(g)^T\, J_s(g)

Integration rules: PLANE233 uses Gauss quadrature; e.g. 3x3 for the quad coefficient matrix
(and 2x2 for axisymmetric massive conductor option). [3]_

A-V with electric potential (:math:`\varphi`) (non time-integrated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Define the **conductive (damping-type) matrix** for the :math:`\partial_t A` term:

.. math::

   \mathbf{C}_{AA}^e =
   \sum_g \omega_g |J_g|\, w_{\text{geom}}(g)\,
   \mathbf{N}(g)^T\, \sigma_g\, \mathbf{N}(g)

Define the **magneto-electric coupling** (A-V) matrix:

.. math::

   \mathbf{K}_{A\varphi}^e =
   \sum_g \omega_g |J_g|\, w_{\text{geom}}(g)\,
   \mathbf{N}(g)^T\, \sigma_g\, \mathbf{N}(g)\, \frac{1}{\ell_g}

Define the **electric (V-V) matrix**:

.. math::

   \mathbf{K}_{\varphi\varphi}^e =
   \sum_g \omega_g |J_g|\, w_{\text{geom}}(g)\,
   \mathbf{N}(g)^T\, \sigma_g\, \mathbf{N}(g)\, \frac{1}{\ell_g^2}

Then the element equations are

.. math::

   \mathbf{K}_{AA}^e\, \mathbf{a}_e
   + \mathbf{C}_{AA}^e\, \dot{\mathbf{a}}_e
   + \mathbf{K}_{A\varphi}^e\, \boldsymbol{\varphi}_e
   =
   \mathbf{f}_A^e

.. math::

   (\mathbf{K}_{A\varphi}^e)^T\, \dot{\mathbf{a}}_e
   + \mathbf{K}_{\varphi\varphi}^e\, \boldsymbol{\varphi}_e
   =
   \mathbf{f}_\varphi^e

where :math:`\mathbf{f}_\varphi^e` represents the imposed total current excitation (or is
zero if voltage-driven and current is a reaction).

**Specialization of :math:`\ell_g` and :math:`w_{\text{geom}}(g)`**

* planar: :math:`w_{\text{geom}}=t`, :math:`\ell_g=t` (factors simplify)
* axisym: :math:`w_{\text{geom}}=2\pi f r_g`, :math:`\ell_g=2\pi f r_g`

This matches the PLANE233 statement that VOLT represents a single out-of-plane voltage drop,
and that the same THK real constant controls thickness (plane) vs fraction of :math:`360^\circ`
(axisym). [1]_

A-V with time-integrated potential (:math:`\eta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :math:`\boldsymbol{\varphi}_e = \dot{\boldsymbol{\eta}}_e`. The same matrices
:math:`\mathbf{C}_{AA}^e`, :math:`\mathbf{K}_{A\varphi}^e`,
:math:`\mathbf{K}_{\varphi\varphi}^e` appear, but they multiply time derivatives of
:math:`\eta`:

.. math::

   \mathbf{K}_{AA}^e\, \mathbf{a}_e
   + \mathbf{C}_{AA}^e\, \dot{\mathbf{a}}_e
   + \mathbf{K}_{A\varphi}^e\, \dot{\boldsymbol{\eta}}_e
   =
   \mathbf{f}_A^e

.. math::

   (\mathbf{K}_{A\varphi}^e)^T\, \dot{\mathbf{a}}_e
   + \mathbf{K}_{\varphi\varphi}^e\, \dot{\boldsymbol{\eta}}_e
   =
   \mathbf{f}_\eta^e

This is the documented **strong coupling with time-integrated VOLT**, yielding a symmetric
matrix for the electromagnetic option. [1]_

Block element matrix forms
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let the element unknown vector be:

* **non-integrated:** :math:`\mathbf{x}_e = [\mathbf{a}_e; \boldsymbol{\varphi}_e]`
* **time-integrated:** :math:`\mathbf{x}_e = [\mathbf{a}_e; \boldsymbol{\eta}_e]`

A-only
^^^^^^

.. math::

   \mathbf{K}_{AA}^e\, \mathbf{a}_e = \mathbf{f}_A^e
   \quad \text{(static)}

or with transient source (still no eddy currents in the "stranded" behavior case): same form
but time-dependent loads.

A-V with electric potential (:math:`\varphi`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A convenient first-order block form:

.. math::

   \underbrace{\begin{bmatrix}
   \mathbf{C}_{AA}^e & \mathbf{0} \\
   (\mathbf{K}_{A\varphi}^e)^T & \mathbf{0}
   \end{bmatrix}}_{\mathbf{D}^e}
   \begin{bmatrix}
   \dot{\mathbf{a}}_e \\
   \dot{\boldsymbol{\varphi}}_e
   \end{bmatrix}
   +
   \underbrace{\begin{bmatrix}
   \mathbf{K}_{AA}^e & \mathbf{K}_{A\varphi}^e \\
   \mathbf{0} & \mathbf{K}_{\varphi\varphi}^e
   \end{bmatrix}}_{\mathbf{K}^e}
   \begin{bmatrix}
   \mathbf{a}_e \\
   \boldsymbol{\varphi}_e
   \end{bmatrix}
   =
   \begin{bmatrix}
   \mathbf{f}_A^e \\
   \mathbf{f}_\varphi^e
   \end{bmatrix}.

(Here :math:`\dot{\boldsymbol{\varphi}}_e` does not appear; the lower-right block of
:math:`\mathbf{D}^e` is zero.)

A-V with time-integrated potential (:math:`\eta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now the coupling is fully in the "damping" (time-derivative) block:

.. math::

   \underbrace{\begin{bmatrix}
   \mathbf{C}_{AA}^e & \mathbf{K}_{A\varphi}^e \\
   (\mathbf{K}_{A\varphi}^e)^T & \mathbf{K}_{\varphi\varphi}^e
   \end{bmatrix}}_{\mathbf{D}^e \; (\text{symmetric})}
   \begin{bmatrix}
   \dot{\mathbf{a}}_e \\
   \dot{\boldsymbol{\eta}}_e
   \end{bmatrix}
   +
   \underbrace{\begin{bmatrix}
   \mathbf{K}_{AA}^e & \mathbf{0} \\
   \mathbf{0} & \mathbf{0}
   \end{bmatrix}}_{\mathbf{K}^e}
   \begin{bmatrix}
   \mathbf{a}_e \\
   \boldsymbol{\eta}_e
   \end{bmatrix}
   =
   \begin{bmatrix}
   \mathbf{f}_A^e \\
   \mathbf{f}_\eta^e
   \end{bmatrix}.

This is the form that corresponds to the documented "time-integrated VOLT" option used to
obtain symmetric coupling. [1]_

Current interpretation (as requested)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **A-only (stranded behavior):** :math:`J \equiv J_s` is **imposed** (body load), eddy currents
  are not solved self-consistently in time varying analyses without VOLT. [2]_
* **A-V (solid/massive conductor):**

  .. math::

     J = J_s - \sigma\left(\partial_t A + \frac{\varphi}{\ell}\right)
     \quad \text{or} \quad
     J = J_s - \sigma\left(\partial_t A + \frac{\partial_t \eta}{\ell}\right),

  so the conduction/induction part is computed **self-consistently**, while :math:`J_s`
  represents any impressed source current density you include.

If you want, I can also add the **harmonic** versions (replace
:math:`\partial_t(\cdot) \mapsto i\omega(\cdot)`) and explicitly show how the axisymmetric
coupling factors simplify because :math:`w_{\text{geom}}/\ell` becomes 1 when both are
:math:`2\pi f r`.

3D A-VOLT formulation (SOLID236 and SOLID97)
--------------------------------------------

This section documents the 3D (full-volume) counterpart of the :math:`\mathbf{A}`-VOLT
formulation for:

* **SOLID236** (edge-based electromagnetic solid, Nedelec edge-flux DOFs + nodal scalar
  potential). [4]_
* **SOLID97** (nodal electromagnetic solid, nodal vector potential DOFs + time-integrated
  scalar potential), legacy technology included for completeness. [5]_

Both non-time-integrated electric potential :math:`\varphi` and time-integrated potential
:math:`\eta(t)=\int_0^t \varphi(\tau)\, d\tau` (so :math:`\varphi=\dot{\eta}`) are shown.
SOLID236 supports both meanings of VOLT via options. [4]_ For SOLID97, comparisons report the
VOLT DOF is **always time-integrated**; the non-integrated variant is shown for theoretical
parity. [6]_

Continuous 3D MQS model (common to both element technologies)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\Omega \subset \mathbb{R}^3`. Unknowns:

* magnetic vector potential :math:`\mathbf{A}(\mathbf{x},t)`
* electric scalar potential :math:`\varphi(\mathbf{x},t)` (non-integrated) or
  :math:`\eta(\mathbf{x},t)` (time-integrated)

Fields and constitutive laws
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{B} = \nabla\times \mathbf{A},
   \qquad
   \mathbf{H} = \boldsymbol{\nu}\, \mathbf{B}

.. math::

   \mathbf{E} = -\partial_t \mathbf{A} - \nabla \varphi
   \quad \text{or} \quad
   \mathbf{E} = -\partial_t \mathbf{A} - \nabla \dot{\eta}

**Currents**

* **Stranded / imposed current density:** :math:`\mathbf{J}=\mathbf{J}_s` (applied as a
  source/body load).
* **Solid conductor (eddy currents):** the **total** current density is computed
  self-consistently as

  .. math::

     \mathbf{J} = \mathbf{J}_s + \mathbf{J}_c,
     \qquad
     \mathbf{J}_c = \boldsymbol{\sigma}\, \mathbf{E}
     = -\boldsymbol{\sigma}\left(\partial_t \mathbf{A} + \nabla \varphi\right)
     \quad \text{or} \quad
     -\boldsymbol{\sigma}\left(\partial_t \mathbf{A} + \nabla \dot{\eta}\right)

  This is the standard :math:`\mathbf{A}`-VOLT formulation used for massive conductors and
  terminal conditions. [7]_

Governing equations (magneto-quasistatic, no displacement current)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \nabla\times\mathbf{H} = \mathbf{J}
   \quad\Rightarrow\quad
   \nabla\times(\boldsymbol{\nu}\,\nabla\times\mathbf{A})
   = \mathbf{J}_s - \boldsymbol{\sigma}(\partial_t\mathbf{A}+\nabla\varphi)

and charge conservation in conducting regions:

.. math::

   \nabla\cdot\mathbf{J} = 0
   \quad\Rightarrow\quad
   \nabla\cdot\left(\boldsymbol{\sigma}(\partial_t\mathbf{A}+\nabla\varphi)\right)
   = \nabla\cdot\mathbf{J}_s

Gauge / uniqueness
^^^^^^^^^^^^^^^^^^

The :math:`\mathbf{A}` formulation needs a gauge constraint to eliminate non-uniqueness.

* **SOLID236 (edge-based)** uses **tree gauging** (a discrete gauge) to produce a unique
  solution. [4]_
* **SOLID97 (nodal)** typically relies on a **Coulomb-type gauge** (or equivalent
  stabilization/constraint). [5]_

To keep matrix expressions explicit, include a generic gauge bilinear form (zero for the
edge/tree-gauged case):

.. math::

   g(\mathbf{A},\mathbf{w}) =
   \int_\Omega \gamma\, (\nabla\cdot\mathbf{w})(\nabla\cdot\mathbf{A})\, d\Omega,

with :math:`\gamma \ge 0` a stabilization/penalty parameter.

Weak formulations
~~~~~~~~~~~~~~~~~

Let :math:`\mathbf{w}` be a vector test function for :math:`\mathbf{A}`; let :math:`q` be a
scalar test function for :math:`\varphi` (or :math:`\eta`).

Non-time-integrated electric potential (:math:`\varphi`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Find :math:`(\mathbf{A},\varphi)` such that for all :math:`(\mathbf{w},q)`:

**A-equation**

.. math::

   \int_\Omega (\nabla\times\mathbf{w})^T \boldsymbol{\nu}
   (\nabla\times\mathbf{A})\, d\Omega
   + g(\mathbf{A},\mathbf{w})
   + \int_\Omega \mathbf{w}^T \boldsymbol{\sigma}\, \partial_t\mathbf{A}\, d\Omega
   + \int_\Omega \mathbf{w}^T \boldsymbol{\sigma}\, \nabla\varphi\, d\Omega
   =
   \int_\Omega \mathbf{w}^T \mathbf{J}_s\, d\Omega

**VOLT-equation (conduction continuity)**

.. math::

   \int_\Omega (\nabla q)^T \boldsymbol{\sigma}\, \partial_t\mathbf{A}\, d\Omega
   + \int_\Omega (\nabla q)^T \boldsymbol{\sigma}\, \nabla\varphi\, d\Omega
   =
   \int_\Omega q\, (\nabla\cdot \mathbf{J}_s)\, d\Omega

In typical terminal-condition usage, :math:`\varphi` is constrained/coupled on conductor
terminals and the right-hand side is assembled as net injected current (or :math:`\mathbf{J}_s`
is used for impressed source regions). This matches the documented use of VOLT for current-fed
massive conductors and terminal enforcement. [8]_

Time-integrated electric potential (:math:`\eta`) (VOLT = :math:`\eta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Find :math:`(\mathbf{A},\eta)` such that for all :math:`(\mathbf{w},q)`:

**A-equation**

.. math::

   \int_\Omega (\nabla\times\mathbf{w})^T \boldsymbol{\nu}
   (\nabla\times\mathbf{A})\, d\Omega
   + g(\mathbf{A},\mathbf{w})
   + \int_\Omega \mathbf{w}^T \boldsymbol{\sigma}\, \partial_t\mathbf{A}\, d\Omega
   + \int_\Omega \mathbf{w}^T \boldsymbol{\sigma}\, \nabla\dot{\eta}\, d\Omega
   =
   \int_\Omega \mathbf{w}^T \mathbf{J}_s\, d\Omega

**VOLT-equation**

.. math::

   \int_\Omega (\nabla q)^T \boldsymbol{\sigma}\, \partial_t\mathbf{A}\, d\Omega
   + \int_\Omega (\nabla q)^T \boldsymbol{\sigma}\, \nabla\dot{\eta}\, d\Omega
   =
   \int_\Omega q\, (\nabla\cdot \mathbf{J}_s)\, d\Omega

This time-integrated VOLT option is described in the ANSYS theory/element documentation as a
way to obtain symmetric coupling (versus the true-voltage nonsymmetric form). [9]_

Discretization: SOLID236 (edge-based) vs SOLID97 (nodal)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edge-based (SOLID236): edge flux DOFs + nodal scalar potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Magnetic unknown (edge fluxes).** Each magnetic DOF is a line integral of
:math:`\mathbf{A}` along an element edge :math:`L_k`:

.. math::

   a_k = \int_{L_k}\mathbf{A}\cdot d\mathbf{l}

The FE approximation uses edge (vector) shape functions :math:`\mathbf{W}_k(\mathbf{x})`
(Nedelec type) such that

.. math::

   \mathbf{A}_h(\mathbf{x}) = \sum_{k=1}^{n_e}\mathbf{W}_k(\mathbf{x})\, a_k
   \equiv \mathbf{W}(\mathbf{x})\, \mathbf{a}_e

and the scalar potential uses nodal Lagrange shapes :math:`N_i`:

.. math::

   \varphi_h = \mathbf{N}\, \boldsymbol{\varphi}_e,
   \qquad
   \eta_h = \mathbf{N}\, \boldsymbol{\eta}_e

Define

.. math::

   \mathbf{C}(\mathbf{x}) \equiv \nabla\times\mathbf{W}(\mathbf{x}),
   \qquad
   \mathbf{G}(\mathbf{x}) \equiv \nabla\mathbf{N}(\mathbf{x})

so

.. math::

   \mathbf{B}_h = \nabla\times\mathbf{A}_h = \mathbf{C}\, \mathbf{a}_e,
   \qquad
   \nabla\varphi_h = \mathbf{G}\, \boldsymbol{\varphi}_e,
   \qquad
   \nabla\dot{\eta}_h = \mathbf{G}\, \dot{\boldsymbol{\eta}}_e

**Gauge:** SOLID236 uses tree gauging rather than a volumetric penalty term. [4]_ So for
SOLID236 set :math:`g(\cdot,\cdot)=0` in the element matrix expressions.

Nodal-based (SOLID97): nodal vector potential + scalar potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Interpolate :math:`\mathbf{A}` with standard nodal Lagrange shapes:

.. math::

   \mathbf{A}_h(\mathbf{x}) = \sum_{i=1}^{n} N_i(\mathbf{x})\, \mathbf{A}_i
   \equiv \mathbf{N}_A(\mathbf{x})\, \mathbf{a}_e

where :math:`\mathbf{a}_e=[A_{x1},A_{y1},A_{z1},\ldots,A_{xn},A_{yn},A_{zn}]^T` and

.. math::

   \mathbf{N}_A(\mathbf{x}) =
   \begin{bmatrix}
   N_1 & 0 & 0 & \cdots & N_n & 0 & 0 \\
   0 & N_1 & 0 & \cdots & 0 & N_n & 0 \\
   0 & 0 & N_1 & \cdots & 0 & 0 & N_n
   \end{bmatrix}

Then

.. math::

   \mathbf{B}_h = (\nabla\times\mathbf{N}_A)\, \mathbf{a}_e \equiv \mathbf{C}_A\, \mathbf{a}_e,
   \qquad
   \nabla\varphi_h = \mathbf{G}\, \boldsymbol{\varphi}_e,
   \qquad
   \nabla\dot{\eta}_h = \mathbf{G}\, \dot{\boldsymbol{\eta}}_e

**Gauge:** nodal :math:`\mathbf{A}` formulations typically require explicit
gauge/stabilization; a Coulomb gauge term
:math:`\propto \int (\nabla\cdot\mathbf{A})(\nabla\cdot\mathbf{w})` is commonly used in
literature discussing SOLID97. [5]_

**Important practical note:** comparisons report SOLID97's VOLT DOF is **time-integrated**
(not true :math:`\varphi`). [6]_ The non-integrated variant below is mainly for theoretical
parity with SOLID236.

Element matrices (quadrature form)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`g` index Gauss points with weights :math:`\omega_g`, Jacobian determinants
:math:`|J_g|`. Material tensors :math:`\boldsymbol{\nu}_g` and :math:`\boldsymbol{\sigma}_g` are
evaluated at Gauss points. The 3D MQS form is written (no dielectric/displacement current).

SOLID236 (edge-based) element matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Core magnetic stiffness
^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{K}_{AA}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{C}(g)^T\, \boldsymbol{\nu}_g\, \mathbf{C}(g)

(\ :math:`\mathbf{C}=\nabla\times\mathbf{W}`, mapping edge flux DOFs to :math:`\mathbf{B}`.)

Conductive (eddy-current) damping on :math:`\partial_t\mathbf{A}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{C}_{AA}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{W}(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{W}(g)

Magneto-electric coupling (A-V)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From :math:`\int \mathbf{w}^T \boldsymbol{\sigma}\nabla\varphi`:

.. math::

   \mathbf{K}_{A\varphi}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{W}(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{G}(g)

with :math:`\mathbf{G}=\nabla\mathbf{N}`.

Electric (V-V) conductivity matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{K}_{\varphi\varphi}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{G}(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{G}(g)

Source vector
^^^^^^^^^^^^^

.. math::

   \mathbf{f}_A^e
   =
   \sum_g \omega_g |J_g|\, \mathbf{W}(g)^T\, \mathbf{J}_s(g)

(For stranded/imposed current density regions, :math:`\mathbf{J}_s` is the imposed excitation;
for solid conductors you typically excite via terminals/total current, plus possibly
:math:`\mathbf{J}_s` in source regions.) [4]_

Block element equations: non-integrated VOLT (:math:`\varphi`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Element unknown vector :math:`\mathbf{x}_e=[\mathbf{a}_e;\boldsymbol{\varphi}_e]`. The
semi-discrete system is

.. math::

   \begin{bmatrix}
   \mathbf{C}_{AA}^e & \mathbf{0} \\
   (\mathbf{K}_{A\varphi}^e)^T & \mathbf{0}
   \end{bmatrix}
   \begin{bmatrix}
   \dot{\mathbf{a}}_e \\
   \dot{\boldsymbol{\varphi}}_e
   \end{bmatrix}
   +
   \begin{bmatrix}
   \mathbf{K}_{AA}^e & \mathbf{K}_{A\varphi}^e \\
   \mathbf{0} & \mathbf{K}_{\varphi\varphi}^e
   \end{bmatrix}
   \begin{bmatrix}
   \mathbf{a}_e \\
   \boldsymbol{\varphi}_e
   \end{bmatrix}
   =
   \begin{bmatrix}
   \mathbf{f}_A^e \\
   \mathbf{f}_\varphi^e
   \end{bmatrix}

This is the **strong coupling, generally nonsymmetric** true electric potential form
described for edge-based elements. [9]_

Block element equations: time-integrated VOLT (:math:`\eta`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now :math:`\boldsymbol{\varphi}_e = \dot{\boldsymbol{\eta}}_e`. With
:math:`\mathbf{x}_e=[\mathbf{a}_e;\boldsymbol{\eta}_e]`,

.. math::

   \begin{bmatrix}
   \mathbf{C}_{AA}^e & \mathbf{K}_{A\varphi}^e \\
   (\mathbf{K}_{A\varphi}^e)^T & \mathbf{K}_{\varphi\varphi}^e
   \end{bmatrix}
   \begin{bmatrix}
   \dot{\mathbf{a}}_e \\
   \dot{\boldsymbol{\eta}}_e
   \end{bmatrix}
   +
   \begin{bmatrix}
   \mathbf{K}_{AA}^e & \mathbf{0} \\
   \mathbf{0} & \mathbf{0}
   \end{bmatrix}
   \begin{bmatrix}
   \mathbf{a}_e \\
   \boldsymbol{\eta}_e
   \end{bmatrix}
   =
   \begin{bmatrix}
   \mathbf{f}_A^e \\
   \mathbf{f}_\eta^e
   \end{bmatrix}

This is the **symmetric time-integrated** coupling option referenced for SOLID236/237. [9]_

SOLID97 (nodal) element matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here :math:`\mathbf{A}` is nodal-interpolated, so replace :math:`\mathbf{W}` by
:math:`\mathbf{N}_A` in the conductive terms, and use
:math:`\mathbf{C}_A=\nabla\times\mathbf{N}_A` in the magnetic stiffness.

Magnetic stiffness (+ optional gauge stabilization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{K}_{AA}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{C}_A(g)^T\, \boldsymbol{\nu}_g\, \mathbf{C}_A(g)
   +
   \sum_g \omega_g |J_g|\,
   \gamma\, \mathbf{D}_A(g)^T \mathbf{D}_A(g)

where :math:`\mathbf{D}_A(g)` maps nodal :math:`\mathbf{a}_e` to
:math:`(\nabla\cdot\mathbf{A}_h)(g)`. (For edge/tree-gauged SOLID236 this term is not used; for
nodal formulations it is a common way to enforce a Coulomb-type gauge.) [5]_

Conductive damping on :math:`\partial_t\mathbf{A}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{C}_{AA}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{N}_A(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{N}_A(g)

A-V coupling and V-V matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf{K}_{A\varphi}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{N}_A(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{G}(g)

.. math::

   \mathbf{K}_{\varphi\varphi}^e
   =
   \sum_g \omega_g |J_g|\,
   \mathbf{G}(g)^T\, \boldsymbol{\sigma}_g\, \mathbf{G}(g)

Source vector
^^^^^^^^^^^^^

.. math::

   \mathbf{f}_A^e
   =
   \sum_g \omega_g |J_g|\, \mathbf{N}_A(g)^T\, \mathbf{J}_s(g)

Block forms (same algebra as SOLID236)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Non-integrated** :math:`\varphi` block form is identical in structure to the SOLID236 form
  above with :math:`\mathbf{W}\to\mathbf{N}_A` and :math:`\mathbf{C}\to\mathbf{C}_A`.
* **Time-integrated** :math:`\eta` block form is identical in structure to the SOLID236 form
  above, yielding a symmetric time-integrated coupling matrix.

**But:** in practical APDL usage, SOLID97 is widely reported as having **time-integrated VOLT
only**, while SOLID236/237 allow choosing true vs time-integrated VOLT via options. [6]_

Interpretation vs conductor cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Imposed** :math:`\mathbf{J}_s` (stranded/source regions): :math:`\mathbf{J}_s` is applied as
  an excitation; the formulation returns fields and induced effects (if conductive coupling is
  active).
* **Solid conductors (massive):** the FE system computes

  .. math::

     \mathbf{J}_c = -\boldsymbol{\sigma}(\partial_t\mathbf{A}+\nabla\varphi)
     \quad \text{or} \quad
     \mathbf{J}_c = -\boldsymbol{\sigma}(\partial_t\mathbf{A}+\nabla\dot{\eta}),

  so eddy currents, shielding, and terminal conditions come out **self-consistently**. This is
  exactly what the A-VOLT option is for. [8]_

.. [1] https://ansyshelp.ansys.com/public//Views/Secured/corp/v251/en/ans_elem/Hlp_E_PLANE233.html
   "PLANE233"
.. [2] https://ansyshelp.ansys.com/public///////Views/Secured/corp/v251/en/ans_lof/Hlp_G_ELE2_2a.html
   "2.2. Older vs. Current 2D Magnetic Element Technologies"
.. [3] https://ansyshelp.ansys.com/public/////////Views/Secured/corp/v242/en/ans_thry/thy_el233.html
   "13.233. PLANE233 - 2D 8-Node Electromagnetic Solid"
.. [4] https://ansyshelp.ansys.com/public//Views/Secured/corp/v251/en/ans_elem/Hlp_E_SOLID236.html
   "SOLID236"
.. [5] https://arxiv.org/pdf/1802.01701
   "3D Numerical Modeling of Liquid Metal Turbulent Flow in ..."
.. [6] https://pure.mpg.de/rest/items/item_3291644/component/file_3292241/content
   "Specific Features of Eddy Current Analysis with ANSYS ..."
.. [7] https://oss.jishulink.com/caenet/forums/upload/2013/11/16/389/156790308906700.pdf
   "ANSYS Low-Frequency Electromagnetic Analysis Guide"
.. [8] https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/pdf/ANSYS_Mechanical_APDL_Low-Frequency_Electromagnetic_Analysis_Guide.pdf
   "Low-Frequency Electromagnetic Analysis Guide"
.. [9] https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_thry/thy_emg2.html
   "5.2. Derivation of Electromagnetic Matrices"
