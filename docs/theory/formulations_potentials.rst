.. _theory-potentials:

Scalar potential, A, T-A and H-phi formulations
================================================

Besides the A-V formulations (:doc:`formulation_2d`, :doc:`formulation_3d`), the element library has four
formulations that are common in the literature and in applied superconductivity: the magnetic scalar potential for
magnetostatics, the vector potential alone in the conductors (A, or A*), the current vector potential in the
conductors coupled with the magnetic vector potential (T-A), and the magnetic field in the conductors coupled with a
scalar potential outside them (H-phi). This page derives their weak forms; the element types are listed in
:doc:`../user_manual/element_library`. The notation is that of :doc:`notation`: :math:`\nu = \mu^{-1}` the reluctivity,
:math:`\rho = \sigma^{-1}` the resistivity, and the material law :math:`\mathbf{B} = \mu(\mathbf{H} + \mathbf{M})`.

.. contents:: On this page
   :local:
   :depth: 1

Scalar potential: ``PHI3D_nodal``, ``PHI2D_nodal``
--------------------------------------------------

Where there is no current, :math:`\nabla\times\mathbf{H} = \mathbf{0}`, the field is a gradient. With a source field
:math:`\mathbf{H}_s` of the currents (Biot-Savart, :math:`\nabla\times\mathbf{H}_s = \mathbf{J}_s`):

.. math::

   \mathbf{H} = \mathbf{H}_s - \nabla\phi, \qquad \nabla\cdot\mathbf{B} = \nabla\cdot\mu(\mathbf{H}_s + \mathbf{M} -
   \nabla\phi) = 0.

With :math:`\mathbf{H}_s = 0`, :math:`\phi` is the **total** potential; with :math:`\mathbf{H}_s \ne 0`, the **reduced**
potential [SimkinTrowbridge1979]_. The Galerkin weak form, test function :math:`q`:

.. math::

   \int \nabla q\cdot\mu\nabla\phi = \int \nabla q\cdot\mu\,\mathbf{H}_s + \int\nabla q\cdot\mu\,\mathbf{M}
   - \oint q\,\mathbf{n}\cdot\mathbf{B}.

The boundary term is the natural condition: a free boundary is a flux wall (:math:`\mathbf{n}\cdot\mathbf{B} = 0`),
and a fixed :math:`\phi` imposes the tangential field (:math:`\phi` constant: :math:`\mathbf{H}\times\mathbf{n} =
\mathbf{H}_s\times\mathbf{n}`, a flux-normal plane or an ideal iron pole). A uniform applied field :math:`\mathbf{H}_0`
is :math:`\phi = -\mathbf{H}_0\cdot\mathbf{x}` on the outer boundary.

* The reduced potential suffers from cancellation in iron (:math:`\mathbf{H}_s` and :math:`\nabla\phi` large and
  nearly opposite, [Mayergoyz1987]_); the remedy is the two-potential method (total potential in the iron, reduced
  elsewhere), which needs two element types and is not in the library.
* In a region that is not simply connected around a current (a coil threading the domain), the total potential is
  multivalued: use the reduced potential there.
* The formulation is static (``time_order`` 0) and linear; the B-H law of the library replaces :math:`\nu`, not
  :math:`\mu`, so it does not apply.

``PHI2D_nodal`` is the same weak form on a 2D cell, planar or axisymmetric about y: the out-of-plane length
(``thickness``, times :math:`2\pi r`) weights the integrals as in ``AV2D_nodal``, and :math:`\phi` varies in the plane.
It solves an axisymmetric magnetostatic problem whose fields lie in the meridian plane.

Vector potential alone in the conductors: ``A3D_edge``
------------------------------------------------------

The A* (modified vector potential) formulation [EmsonSimkin1983]_, [BiroPreis1989]_ absorbs the electric potential
into the vector potential: :math:`\mathbf{E} = -\partial_t\mathbf{A}`, so that in the conductors

.. math::

   \nabla\times\nu\nabla\times\mathbf{A} + \sigma\,\partial_t\mathbf{A} = \mathbf{J}_s + \nabla\times\mathbf{M},

and the conduction current :math:`-\sigma\partial_t\mathbf{A}` is divergence free by the equation itself (its
divergence is that of a curl). The weak form is ``AV3D_edge`` without V:

.. math::

   \int\nabla\times\mathbf{w}\cdot\nu\nabla\times\mathbf{A} + \int_{\Omega_c}\mathbf{w}\cdot\sigma\partial_t\mathbf{A}
   = \int\mathbf{w}\cdot\mathbf{J}_s + \int\nabla\times\mathbf{w}\cdot(\mathbf{M} + \mathbf{T}_0).

* **Gauge.** In the conductors the :math:`\sigma` term fixes A completely (the gauge is implicit); in the
  non-conducting region A is defined up to gradients of functions that vanish on the conductor surface. The
  tree-cotree gauge must therefore span the **non-conducting edges only**, rooted on the Dirichlet boundary and on the
  conductor surfaces (``tree_gauge(mesh, non_conducting_point_ids, start_surface_point_ids)`` with the points of both
  in ``start_surface_point_ids``; a surface edge then never enters the tree). A tree through the conductors would
  over-constrain A there.
* **Edge elements only.** A nodal A* is wrong at conductor interfaces (the normal component of A jumps with the
  surface charges); the library has no nodal variant.
* **Currents.** Only the induced currents and imposed sources: no terminal voltage, no injected current. A conductor
  carrying a net current needs the stranded coil domain or the A-V types.

Source currents as a current vector potential: ``T0``
-----------------------------------------------------

``A3D_edge`` and ``TA3D_edge`` take, besides :math:`\mathbf{J}_s`, a source **current vector potential**
:math:`\mathbf{T}_0` (cell array ``T0``, every domain, default 0), :math:`\nabla\times\mathbf{T}_0 = \mathbf{J}_s`, in
the term :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{T}_0`. Its point is compatibility: with a tree gauge the
curl-curl system is solvable only if its right-hand side is orthogonal to the gradients, which
:math:`\int\mathbf{w}\cdot\mathbf{J}_s` is only when the discrete :math:`\mathbf{J}_s` is exactly divergence free,
whereas :math:`\int\nabla\times\mathbf{w}\cdot\mathbf{T}_0` vanishes for every gradient test function
:math:`\mathbf{w} = \nabla q`, whatever :math:`\mathbf{T}_0` is [Ren1996]_. A piecewise
constant :math:`\mathbf{T}_0` is enough; the Biot-Savart field of the coil, :math:`\mathbf{H}_s`, is one.

* Give a source either as ``Js`` or as ``T0``, not both (they would add).
* Integrated by parts, the term imposes :math:`\nabla\times\mathbf{T}_0` inside and the surface current
  :math:`\mathbf{T}_0\times\mathbf{n}` where :math:`\mathbf{T}_0` jumps; on a boundary with the natural condition it
  also imposes :math:`\mathbf{n}\times\mathbf{H} = \mathbf{n}\times\mathbf{T}_0`: :math:`\mathbf{T}_0` should vanish
  there, or be normal to it (a symmetry plane of the coil).
* It enters exactly as the magnetization: a uniform :math:`\mathbf{T}_0` in a ball is a uniformly magnetized ball,
  :math:`\mathbf{B} = \tfrac23\mu_0\mathbf{T}_0` inside (the ``t0_sphere`` benchmark).

Current vector potential in the conductors: ``TA3D_edge``
---------------------------------------------------------

The T-A formulation [Zhang2017]_, [Liang2017]_, [Huber2022]_ describes the current in the conductors by a current
vector potential, :math:`\mathbf{J} = \nabla\times\mathbf{T}` (divergence free by construction, no electric potential),
and the field everywhere by A. Ampère's law in the whole domain and Faraday's law in the conductors,
:math:`\nabla\times\mathbf{E} = -\partial_t\mathbf{B}` with :math:`\mathbf{E} = \rho\mathbf{J}`:

.. math::

   \nabla\times\nu\nabla\times\mathbf{A} = \nabla\times\mathbf{T} + \mathbf{J}_s + \ldots, \qquad
   \nabla\times\rho\nabla\times\mathbf{T} + \nabla\times\partial_t\mathbf{A} = \mathbf{0}\ \text{in}\ \Omega_c.

No current crosses the conductor surface, :math:`\mathbf{n}\cdot\mathbf{J} = 0`: :math:`\mathbf{n}\times\mathbf{T} =
\mathbf{0}` there. The weak form, test functions :math:`\mathbf{w}` (A) and :math:`\mathbf{t}` (T, zero tangential
trace on the surface):

.. math::

   \int\nabla\times\mathbf{w}\cdot\nu\nabla\times\mathbf{A} - \int_{\Omega_c}\mathbf{w}\cdot\nabla\times\mathbf{T}
   &= \int\mathbf{w}\cdot\mathbf{J}_s + \int\nabla\times\mathbf{w}\cdot(\mathbf{M} + \mathbf{T}_0), \\
   \int_{\Omega_c}\nabla\times\mathbf{t}\cdot\rho\nabla\times\mathbf{T} + \int_{\Omega_c}\nabla\times\mathbf{t}\cdot
   \partial_t\mathbf{A} &= 0.

The second row is :math:`\int\nabla\times\mathbf{t}\cdot\mathbf{E} = 0` with :math:`\mathbf{E} = \rho\mathbf{J} =
-\partial_t\mathbf{A} - \nabla V`: the gradient of the electric potential drops out because
:math:`\int\nabla\times\mathbf{t}\cdot\nabla V = \oint V\,\mathbf{n}\cdot\nabla\times\mathbf{t} = 0` for a
:math:`\mathbf{t}` with zero tangential trace. Both fields are first-order Nedelec, so
:math:`\mathbf{J} = \nabla\times\mathbf{T}` is face-conforming.

* **Gauges.** A appears only through its curl and :math:`\partial_t\mathbf{A}` only against curls: A is defined up to
  gradients **everywhere**, conductors included, so its tree spans all the edges (rooted on the Dirichlet boundary).
  T is defined up to gradients of functions vanishing on the surface: fix its surface edges (``_T`` = 0) and gauge its
  inner edges with a tree of the conductor rooted on its surface.
* **Topology.** With :math:`\mathbf{n}\times\mathbf{T} = 0` on the whole surface, the net current through any
  cross-section of a conductor is zero (it is the circulation of T along the section boundary). The volumetric T-A
  of this library models induced currents in simply connected conductors; a net transport current or a ring
  (multiply connected) needs cuts or a T jump across a cut, as in the thin-shell T-A of the literature, which embeds
  T on surfaces and is not available here.
* **Superconductors.** The formulation is meant for nonlinear :math:`\rho(\mathbf{J})` (power law); the library
  has the linear resistivity only, the nonlinear law being a later addition (:doc:`nonlinear`).

Magnetic field in the conductors: ``HPHI3D_edge``
-------------------------------------------------

The H-phi formulation [Bossavit1998]_, [Dular2020]_, [Arsenault2021]_ solves for the magnetic field. In the conductors H
is an edge field; outside, where :math:`\nabla\times\mathbf{H} = \mathbf{0}`, it is :math:`-\nabla\phi`. Faraday's law with
:math:`\mathbf{E} = \rho\nabla\times\mathbf{H}`, over the whole domain:

.. math::

   \int\partial_t(\mu\mathbf{H})\cdot\mathbf{H}' + \int_{\Omega_c}\rho\,\nabla\times\mathbf{H}\cdot\nabla\times\mathbf{H}'
   = 0 \qquad \forall\,\mathbf{H}',

where the test fields :math:`\mathbf{H}'` are of the same kind as H (in the non-conducting region,
:math:`-\nabla q`: the equation is then :math:`\nabla\cdot\partial_t\mathbf{B} = 0`).

**Discretization.** The library writes :math:`\mathbf{H} = \mathbf{H}_c - \nabla\phi`, with :math:`\phi` a Lagrange field
over the whole mesh and :math:`\mathbf{H}_c` (``Hc``) a Nedelec field of the conductors with zero tangential trace on
their surface (``_Hc`` = 0 on the surface edges). In the conductors :math:`\phi` is fixed to zero at the inner nodes
(``_PHI`` = 0), so that it only extends the outer potential into the first layer of cells; :math:`\mathbf{H}_c` adds
the rest. This decomposition is unique and spans the H-phi space; it is the T-Omega decomposition
[Carpenter1977]_ with :math:`\mathbf{T} = \mathbf{H}_c`. Inserting it, with the test functions
:math:`\mathbf{h}_c` and :math:`-\nabla q`:

.. math::

   \int_{\Omega_c}\rho\nabla\times\mathbf{H}_c\cdot\nabla\times\mathbf{h}_c + \int_{\Omega_c}\mu\,\partial_t
   (\mathbf{H}_c - \nabla\phi)\cdot\mathbf{h}_c &= 0, \\
   \int\mu\,\partial_t(\nabla\phi - \mathbf{H}_c)\cdot\nabla q &= -\oint q\,\mathbf{n}\cdot\partial_t\mathbf{B}.

The system is symmetric (K holds the resistive term, C the magnetic one).

* :math:`\phi` is **first order** (``"order": "first"`` in the library, whatever ``option2``): the gradients of
  first-order Lagrange functions lie in the first-order Nedelec space, which makes the sum :math:`\mathbf{H}_c -
  \nabla\phi` consistent (the Whitney complex). With second-order :math:`\phi` the uniform field is not even
  represented inside a conductor.
* **Transient only.** :math:`\phi` has no stiffness term: a static step pins it to zero (:doc:`implementation`) and
  :math:`\mathbf{H}_c` is then undetermined. Start from the zero field and drive the model through time.
* **Sources.** An applied field through :math:`\phi` on the outer boundary (:math:`\phi = -\mathbf{H}_0(t)\cdot
  \mathbf{x}`); a free boundary is a flux wall. Source currents (a reduced field :math:`\mathbf{H}_s`) and permanent
  magnets are not available: in this time-differentiated form they would need :math:`\partial_t\mathbf{H}_s`.
* **Topology.** As for the total scalar potential, a conductor carrying a net current makes :math:`\phi`
  multivalued around it: cuts (thick or thin, [Arsenault2021b]_) or cohomology basis functions [Lahtinen2015]_ are
  needed and not available. Induced currents in simply connected conductors need neither.
* Linear :math:`\mu` and :math:`\rho` only (the B-H and power laws are not attached to this type).

Validation
----------

The ``formulations`` suite of ``pyTEC3D_APDL_benchmarks`` checks every type against analytical solutions and, where
MAPDL has the formulation, against MAPDL on the same mesh (:doc:`../user_manual/validation`): a ball of
permeability :math:`\mu_r` in a uniform field (3D, total and reduced potential; 2D axisymmetric), a cylinder (2D
planar), a uniformly magnetized ball, a ball with a uniform :math:`\mathbf{T}_0`, and a conducting ball in a field
switched on at :math:`t = 0`, whose induced moment

.. math::

   m(t) = -2\pi a^3 H_0\,\frac{6}{\pi^2}\sum_{k\ge1}\frac{e^{-k^2\pi^2t/\tau}}{k^2}, \qquad \tau = \mu_0\sigma a^2,

is reproduced by ``AV3D_edge``, ``A3D_edge``, ``TA3D_edge`` and ``HPHI3D_edge``, the box carrying the exact exterior
field. The edge DOF of a Dirichlet array is the circulation of the field along the edge **from its lexicographically
smaller end point to the larger one** (x, then y, then z: libMesh's edge orientation), which is how the benchmarks
impose a non-zero tangential A.
