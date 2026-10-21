Loads
=====

A load is an array value, read at every solve. Loads can change freely between steps (the system matrix is reused).

.. list-table::
   :header-rows: 1
   :widths: 12 10 30 48

   * - Array
     - Where
     - Element types / domains
     - Meaning
   * - ``Js``
     - cells, 3 comp.
     - ``STRANDED_CONDUCTOR``
     - Imposed source current density [A/m^2], in the element coordinate system of the cell. In 2D only component 2
       (out of plane / toroidal) acts.
   * - ``CURR``
     - nodes
     - ``STRANDED_COIL``
     - Load of the VOLT equation: coil current [A]. Put it on the master node of the coupling.
   * - ``AMPS``
     - nodes
     - 2D ``CIRCUIT_COUPLED_CONDUCTOR``; 3D ``SOLID_CONDUCTOR``
     - Load of the V equation. 2D: total current through the conductor, on the master node of its coupling. 3D:
       current injected into the conductor at the node [A] (a negative value extracts it).
   * - ``M``
     - material
     - all but ``NON_CONDUCTING``, ``SOFT_MAGNETIC``
     - Magnetization [A/m] (a material property, constant in time): permanent magnets, magnetized conductors.

Imposed current density (stranded conductors)
---------------------------------------------

The simplest source: a region where the current density is known and no eddy current flows (a coil of thin turns
fed by a current source, a plasma current density given by the plasma algorithm).

.. code-block:: python

   selection.select(mesh, "MATERIAL", 23, location="elem")
   selection.set_value_on_selection(mesh, "Js", 1.0e6, 2, location="elem")    # J_z = 1 MA/m^2

* For a coil of ``N`` turns carrying ``i``, :math:`J_s = N i / S` (S the cross-section of the winding).
* The density must be divergence-free (a closed loop, or a conductor crossing the model boundary): a source with
  :math:`\nabla\cdot\mathbf{J}_s\neq0` has no solution in a static magnetic problem.
* Curved coils: give the components in a cylindrical **element coordinate system** (``COORDINATE_SYSTEM`` of the coil
  cells; Js[1] is then the toroidal component), or compute the global components yourself
  (``model_edit.map_vector_field_to_global_components`` converts an array given in the components of a coordinate
  system; or write them with NumPy from ``cell_centroids``). See ``pyTEC3D_APDL_benchmarks/suites/simple/cases/sphere_coil_3D*``.

Current or voltage of a 2D conductor with a terminal and of a coil
------------------------------------------------------------------

A 2D ``CIRCUIT_COUPLED_CONDUCTOR`` has a terminal, the field V (the voltage across its out-of-plane length); a
``STRANDED_COIL`` has the terminal voltage VOLT and the induced voltage EMF. Each is coupled to a single value over the
conductor / coil (:doc:`couplings`), and the terminal is driven by **one** of:

* the **current**: ``AMPS`` (conductor) or ``CURR`` (coil) on the master node, the voltage free (``_V`` / ``_VOLT`` =
  NaN);
* the **voltage**: ``_V`` / ``_VOLT`` on the master node (the current load is then irrelevant: set it to 0);
* a **short** (conductor): ``_V = 0``, which is the same as a ``SOLID_CONDUCTOR``.

.. code-block:: python

   selection.select(mesh, "ID", master)
   selection.set_value_on_selection(mesh, "AMPS", 100.0)                    # current driven
   ...
   selection.set_value_on_selection(mesh, "AMPS", 0.0)
   selection.set_value_on_selection(mesh, "_V", 1.0e4)                      # voltage driven from now on

A 2D ``SOLID_CONDUCTOR`` has no terminal: it carries the induced current of a closed conductor, nothing to set. For a
coil modelled on a sector (``symmetry`` > 1), ``CURR`` and ``VOLT`` are the current and the voltage of the **whole**
coil. With a ``_timint`` element type ``_V`` and ``_VOLT`` are time integrals of the voltages.

Injected currents (3D conductors)
---------------------------------

In 3D, a current is passed through a conductor by injecting it at one terminal and fixing the potential at the
other:

.. code-block:: python

   selection.select(mesh, "TERMINAL_IN", 1)                  # a point array of the mesh: nodes of a face
   ids = selection.selected_ids(mesh, "node")
   selection.select(mesh, "ID", ids[0])
   selection.set_value_on_selection(mesh, "AMPS", 1000.0)     # 1 kA in at one node (or spread over the face)
   selection.select(mesh, "TERMINAL_OUT", 1)
   selection.set_value_on_selection(mesh, "_V", 0.0)          # out through the grounded face

Couple the V of the terminal face (one CP id, field ``V``) to inject through the whole face as an equipotential. The
plasma algorithm writes the halo currents as ``AMPS`` on the vacuum vessel wall.

Time-dependent loads
--------------------

There are no load tables: the script sets the values before each solve, e.g. from a waveform:

.. code-block:: python

   import numpy as np
   times = np.linspace(0.0, 0.1, 101)[1:]
   for t in times:
       solver.time_step["time"] = t
       selection.select(mesh, "ID", master)
       selection.set_value_on_selection(mesh, "_VOLT", 1.0e3 * np.sin(2 * np.pi * 50 * t))
       solver.solve()

Before each solve, a warning is logged when every load array of the model is zero (nothing drives the problem); it
is harmless when the analysis is driven by Dirichlet values only (e.g. a fixed coil voltage).
