.. _theory-domains:

Domains
-------

Every cell of a model belongs to one **domain**, given by the ``option1`` of the row of its ``TYPE`` in the model data
(:doc:`../user_manual/model_data`). The domain says which physics the cell carries: whether a current can flow in it
and how it is determined, and whether it can be magnetized. The element type then gives the unknowns and the terms of
that physics. pyTEC3D has seven domains (``src/pytec3d/resources/domains.json``):

.. figure:: ../_static/field_domains.svg
   :alt: Domains of an eddy-current model: air, magnetic materials, conductors with imposed or computed current.
   :align: center
   :width: 85%

   The regions of an eddy-current model: air and magnetic materials, conductors with an imposed current density,
   massive conductors where the current is computed, windings fed through terminals.

.. list-table::
   :header-rows: 1
   :widths: 5 25 32 38

   * - #
     - Domain
     - Current
     - Magnetic behaviour
   * - 0
     - ``NON_CONDUCTING``
     - none: :math:`\mathbf{J}=0`
     - linear, :math:`\nu=(\mu_0\mu_r)^{-1}` (air, vacuum, insulators)
   * - 1
     - ``SOFT_MAGNETIC``
     - none
     - linear or saturating (``bh_curve`` law), not magnetized (iron, steel without currents)
   * - 2
     - ``HARD_MAGNETIC``
     - none
     - linear recoil with a magnetization :math:`\mathbf{M}` (permanent magnets)
   * - 3
     - ``STRANDED_CONDUCTOR``
     - imposed: :math:`\mathbf{J}=\mathbf{J}_s` (cell load ``Js``)
     - may be magnetized
   * - 4
     - ``CIRCUIT_COUPLED_CONDUCTOR``
     - computed, eddy currents plus a **terminal** (2D only)
     - may be magnetized
   * - 5
     - ``SOLID_CONDUCTOR``
     - computed: massive conductor with eddy currents
     - may be magnetized
   * - 6
     - ``STRANDED_COIL``
     - imposed shape, amplitude from the **circuit equation** of the winding
     - may be magnetized

The current and the magnetic behaviour are two independent attributes. The current partitions the model into
regions without current (:math:`\Omega_0`: domains 0, 1, 2), with an imposed current density (:math:`\Omega_s`: 3 and,
with a circuit-determined amplitude, 6) and with a computed, eddy current density (:math:`\Omega_c`: 4 and 5). Any
domain may have a saturating B-H curve (the law replaces the reluctivity of the material, :doc:`constitutive`), and
every domain except ``NON_CONDUCTING`` and ``SOFT_MAGNETIC`` may carry a magnetization: a magnetized structure that
also conducts, or a coil wound on a magnetized core, are modelled directly.

Conductors: imposed, computed, through terminals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Stranded conductor.** A region where the current density is known: a coil fed by a current source, the plasma of
the equilibrium. :math:`\mathbf{J}=\mathbf{J}_s` is a load; no conductivity, no eddy current, no electric unknown.

**Solid conductor.** A massive conductor (vacuum vessel, passive rings, plates): the current density follows from the
electric field, :math:`\mathbf{J}=\sigma\mathbf{E}`, :math:`\mathbf{E}=-\partial_t\mathbf{A}-\nabla V`.

* In **2D** the current flows out of plane and a conductor without a terminal has no applied voltage along that
  direction: :math:`\mathbf{J}=-\sigma\,\partial_t\mathbf{A}`, the current induced in a closed ring (axisymmetric) or
  in an infinitely long conductor whose ends are connected at infinity (planar). There is **no electric unknown**.
* In **3D** charge conservation :math:`\nabla\cdot\mathbf{J}=0` is an equation of the model: the electric potential
  :math:`V` is an unknown of every solid conductor. Currents can be injected at its nodes (nodal load ``AMPS``), and
  :math:`V` must be fixed on one node of every conductor not connected to a fixed potential.

**Circuit-coupled conductor (2D).** A massive conductor with a **terminal**: in 2D a voltage :math:`V` applied across
the out-of-plane length of the conductor (the ends of a planar bar, the gap of an axisymmetric ring). :math:`V` is a
DOF, coupled over the whole cross-section (one value per conductor, :doc:`../user_manual/couplings`), and the
conductor is driven by its total current (nodal load ``AMPS`` on the master node) or by its voltage (``_V``). In 3D the
potential exists in every solid conductor, so the two domains coincide: ``CIRCUIT_COUPLED_CONDUCTOR`` is not a domain
of the 3D element types, and a 3D massive conductor is a ``SOLID_CONDUCTOR`` (with its terminals made by ``AMPS`` and
``_V``).

**Stranded coil.** A winding of many thin insulated turns: the current density is uniform over its cross-section,
:math:`\mathbf{J}=(N/S)\,i\,\mathbf{T}`, and the current :math:`i` follows from the circuit equation
:math:`\mathrm{VOLT}=R\,i+\mathrm{EMF}`. The circuit unknowns ``VOLT`` (terminal voltage) and ``EMF`` (induced voltage)
are reserved to this domain (:ref:`theory-stranded-coil`).

Unknowns of each domain
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Domain
     - ``AV2D_nodal`` (2D)
     - ``AV3D_nodal``, ``AV3D_edge`` (3D)
   * - ``NON_CONDUCTING``, ``SOFT_MAGNETIC``, ``HARD_MAGNETIC``
     - :math:`A_z`
     - :math:`\mathbf{A}`
   * - ``STRANDED_CONDUCTOR``
     - :math:`A_z`, load ``Js``
     - :math:`\mathbf{A}`, load ``Js``
   * - ``CIRCUIT_COUPLED_CONDUCTOR``
     - :math:`A_z`, :math:`V` (terminal), load ``AMPS``
     - not a domain of the 3D types
   * - ``SOLID_CONDUCTOR``
     - :math:`A_z`
     - :math:`\mathbf{A}`, :math:`V`, load ``AMPS``
   * - ``STRANDED_COIL``
     - :math:`A_z`, ``VOLT``, ``EMF``, load ``CURR``
     - :math:`\mathbf{A}`, ``VOLT``, ``EMF``, load ``CURR``

In the time-integrated element types (``_timint``) the electric unknowns :math:`V`, ``VOLT`` and ``EMF`` are replaced
by their time integrals (:ref:`theory-timint`). Each element type lists the domains it accepts (``"domains"`` in
``element_types.json``); a model with a row whose domain the element type does not have is rejected when it is
checked.
