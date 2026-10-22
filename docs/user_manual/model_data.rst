The Model Data
==============

The model data says what the cells are: which formulation and physical region each ``TYPE`` id stands for, the
material properties of each ``MATERIAL`` id, the real constants of each ``CONSTANT`` id, which DOFs are coupled, which
nonlinear laws apply and which coordinate systems exist. It is a JSON object, usually kept in ``model_data.json``
next to the mesh and read with ``model.read_json()``; in Python it is the dict ``model.data``.

.. code-block:: json

   {
     "element_types":      {"21": {"element_type": "AV2D_nodal", "option1": "NON_CONDUCTING", "option2": 2, "option3": 1},
                            "23": {"element_type": 0, "option1": 4, "option2": 2, "option3": 1}},
     "materials":          {"21": {"mur": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},
                            "23": {"mur": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "rho": [[1e-7, 0, 0], [0, 1e-7, 0], [0, 0, 1e-7]]}},
     "constants":          {"1": {"thickness": 1.0}},
     "dof_couplings":      {"1": ["V"]},
     "nonlinear_laws":     {},
     "coordinate_systems": {}
   }

Every section is an object keyed by **positive integer ids** (written as strings in JSON). The cells select their
rows through the cell arrays ``TYPE``, ``MATERIAL``, ``CONSTANT`` and ``COORDINATE_SYSTEM``, independently: two cells
with the same ``TYPE`` may have different materials.

In Python the model data can be edited before the model is loaded:

.. code-block:: python

   data = solver.model.data
   data["materials"]["30"] = {"mur": [[1000, 0, 0], [0, 1000, 0], [0, 0, 1000]]}
   solver.model.write_json("model_data_new.json")     # indented JSON
   solver.model.data = {...}                          # replace everything (missing sections are added, empty)

element_types
-------------

.. code-block:: text

   "element_types": {"<TYPE id>": {"element_type": ..., "option1": ..., "option2": 1, "option3": 0}}

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Key
     - Meaning
   * - ``element_type``
     - The formulation: a name of the :doc:`element_library`, or its index (required).

       ===== ===================== ====================================================================
       Index Name                  Formulation
       ===== ===================== ====================================================================
       0     ``AV2D_nodal``        2D planar / axisymmetric, A_z (+ V, VOLT, EMF)
       1     ``AV3D_nodal``        3D, nodal A + V (+ VOLT, EMF), Coulomb gauge term
       2     ``AV3D_edge``         3D, edge (Nedelec) A + nodal V (+ VOLT, EMF)
       3     ``AV2D_nodal_timint`` as 0, with time-integrated electric unknowns
       4     ``AV3D_nodal_timint`` as 1, with time-integrated electric unknowns
       5     ``AV3D_edge_timint``  as 2, with time-integrated electric unknowns
       ===== ===================== ====================================================================
   * - ``option1``
     - The **domain** (physical region): a name or its index (required).

       ===== ============================== =========================================================
       Index Name                           Region
       ===== ============================== =========================================================
       0     ``NON_CONDUCTING``             air, vacuum, insulators: no current
       1     ``SOFT_MAGNETIC``              magnetic material without current (B-H curve possible)
       2     ``HARD_MAGNETIC``              permanent magnet: magnetization M
       3     ``STRANDED_CONDUCTOR``         imposed source current density Js, no eddy currents
       4     ``CIRCUIT_COUPLED_CONDUCTOR``  2D massive conductor with a terminal (V); not a 3D domain
       5     ``SOLID_CONDUCTOR``            massive conductor (2D: no electric DOF; 3D: potential V)
       6     ``STRANDED_COIL``              multi-turn coil with circuit equation (VOLT, EMF)
       ===== ============================== =========================================================
   * - ``option2``
     - FE order of the Lagrange fields, 1 or 2 (default 1). Order 2 needs quadratic cells. Edge fields are always
       first-order Nedelec; ``option2`` then sets the order of V, VOLT and EMF.
   * - ``option3``
     - 1 for a 2D axisymmetric model (about the y axis), 0 otherwise (default 0). Must be the same in the rows used by
       the cells (rows no cell uses are free), and 0 for 3D cells.

Rules checked when the model is loaded:

* every ``TYPE`` value of the cells has a row;
* **all the cells of a model have the same** ``element_type`` (the domains, orders and ids may differ): a model is
  either 2D nodal, 3D nodal or 3D edge, linear or time-integrated;
* the dimension of the element type matches the cells (2D cells, 2D type);
* the domain is one of the domains of the element type (e.g. ``CIRCUIT_COUPLED_CONDUCTOR`` exists only in 2D);
* rows not used by any cell are reported as warnings; keys other than ``element_type`` and ``option1`` - ``option3``
  are errors (nonlinear laws are selected in ``nonlinear_laws``).

materials
---------

.. code-block:: text

   "materials": {"<MATERIAL id>": {"mur": [[...], [...], [...]], "rho": [[...], [...], [...]], "M": [0, 0, 0]}}

.. list-table::
   :header-rows: 1
   :widths: 12 18 70

   * - Name
     - Shape
     - Meaning
   * - ``mur``
     - 3x3 tensor
     - Relative permeability (dimensionless). Needed by every element. The reluctivity is
       ``nu = (mu0 mur)^-1``. An anisotropic tensor is given in the element coordinate system of the cell.
   * - ``rho``
     - 3x3 tensor
     - Resistivity [Ohm m] of the massive conductors (``SOLID_CONDUCTOR``, ``CIRCUIT_COUPLED_CONDUCTOR``);
       conductivity ``sigma = rho^-1``. Ignored in the other domains.
   * - ``M``
     - number or 3-vector
     - Magnetization [A/m] (default 0): ``B = mu0 mur (H + M)``, i.e. ``H = nu B - M``. Read in every domain but
       ``NON_CONDUCTING`` and ``SOFT_MAGNETIC``: permanent magnets (``HARD_MAGNETIC``) and magnetized conductors or
       coils. A number stands for equal components.

Rules: tensors must be full 3x3 arrays (a scalar is not accepted: write ``[[m, 0, 0], [0, m, 0], [0, 0, m]]``); values
finite; unknown names are errors. A property that the element of a cell needs and that its material leaves out is
reported by the solver when the model is loaded (``'rho' is not defined for element type ... in domain ...``). Every
``MATERIAL`` value of the cells must have a row.

constants (real constants)
--------------------------

.. code-block:: json

   "constants": {"<CONSTANT id>": {"thickness": 1.0, "cross_section": 0.04, "turns": 1000, "radius": 0.5,
                                   "resistance": 50.0, "symmetry": 1.0, "T": [0, 0, 1]}}

.. list-table::
   :header-rows: 1
   :widths: 16 12 72

   * - Name
     - Shape
     - Meaning
   * - ``thickness``
     - number
     - 2D only: out-of-plane length factor (default 1). Planar: the depth of the model [m]. Axisymmetric: the
       fraction of the full turn (1 = 360 degrees); the length is ``thickness * 2 pi r``. **0 is used as given** and
       makes the model empty.
   * - ``cross_section``
     - number
     - Stranded coil: cross-section of the winding S [m^2].
   * - ``turns``
     - number
     - Stranded coil: number of turns N.
   * - ``resistance``
     - number
     - Stranded coil: resistance R of the whole coil [Ohm].
   * - ``radius``
     - number
     - 2D stranded coil: mean radius of the coil [m]; the coil volume is ``(axisymmetric ? 2 pi radius : 1) * S *
       thickness``.
   * - ``volume``
     - number
     - 3D stranded coil: volume of the **whole** real coil Vc [m^3] (also when only a sector is modelled).
   * - ``symmetry``
     - number
     - Stranded coil: ratio between the whole coil and the modelled part (1 for a full coil, 12 for a 30 degree
       sector, 2 for a half model). The loads ``CURR`` and ``VOLT`` then refer to the whole coil.
   * - ``T``
     - 3-vector
     - Stranded coil: winding direction (unit vector, in the element coordinate system of the cell).

**Coil volume and symmetry** (:ref:`theory-stranded-coil`, [MAPDLTheory2024]_). The current density of the coil
is ``(N / S) i T`` with ``i = (VOLT - EMF) / R`` and the resistance term of the circuit equation is
``symmetry / (R Vc)`` integrated over the modelled coil cells, Vc the coil volume above. The coil current is the
imposed ``CURR`` exactly when ``symmetry * V_model = Vc``, V_model the volume of the modelled coil cells: ``symmetry``
accounts for the part of the coil that is not modelled. If the mesh does not represent the coil volume exactly (a
curved coil meshed with flat facets), the current density differs from ``N CURR / S`` by the factor
``Vc / (symmetry * V_model)`` (1.8e-3 for a 30 degree sector in 5 straight divisions).

A coil made of straight and curved pieces uses one row per piece (with its own ``T``) and the same terminal coupling
for all of them. ``CONSTANT`` ids used by the cells but absent from the section are warned about; the solver reports
any constant a cell actually needs.

dof_couplings
-------------

.. code-block:: json

   "dof_couplings": {"1": ["VOLT", "EMF"], "2": ["V"], "100": ["A", "V"]}

A CP id and the DOF fields it couples. The nodes that have the same CP id in a point array ``CP_k`` share the value
of those fields (:doc:`couplings`). The fields must be DOF fields of the element type of the model (``A``, ``V``,
``VOLT``, ``EMF``). Every positive CP id found in a ``CP_k`` array must have an entry. In Python:
``model.add_dof_couplings({1: ["VOLT", "EMF"]})``; ``model.periodic_boundaries(...)`` adds its own.

.. note::

   A coupling of a field that does not exist on the coupled nodes does nothing. ``VOLT`` and ``EMF`` exist only on
   stranded coils; the potential of the massive conductors is ``V`` (2D circuit-coupled conductors, 3D solid
   conductors).

nonlinear_laws
--------------

.. code-block:: json

   "nonlinear_laws": {"26": {"type": "bh_curve", "B": [0.0, 0.5, 1.2, 1.8], "H": [0.0, 100.0, 1000.0, 10000.0]}}

Selects a law of the library (``src/pytec3d/resources/nonlinear_laws.json``) on the cells of an id, with its data.
Whether the id is a ``MATERIAL``, ``CONSTANT`` or ``TYPE`` id is given by the *source* of the law (``bh_curve``:
MATERIAL). An id can have a list of entries for several laws. See :doc:`nonlinear`. In Python:
``model.set_bh_curve(26, [[H0, B0], [H1, B1], ...])`` or ``model.set_nonlinear_law(26, "bh_curve", B=[...], H=[...])``.

coordinate_systems
------------------

.. code-block:: json

   "coordinate_systems": {"1": {"type": "y_axis_cylindrical", "origin": [0, 0, 0]},
                          "2": {"type": "cylindrical", "origin": [1, 2, 3], "angles": [0, 0, 90]}}

Element coordinate systems, referenced by the cell array ``COORDINATE_SYSTEM`` (0, the global Cartesian system,
needs no entry). See :doc:`coordinate_systems`. In Python: ``model.set_coordinate_system(1,
YAxisCylindricalCoordinateSystem())``.
