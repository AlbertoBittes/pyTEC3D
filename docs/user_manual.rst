User Manual
===========

This manual documents the practical input rules enforced by model validation.
Use it as a checklist before loading a model into the solver and before solve.

Validation Philosophy
---------------------

pyTEC3D validation has two stages:

1. Model validation (structure and consistency): validates mesh arrays,
	 metadata tables, and cross-references.
2. Pre-solve validation (runtime readiness): validates time-step settings and
	 warns when loads or boundary conditions appear to be missing.

Model validation raises errors for invalid inputs. Some non-fatal situations
are reported as warnings.

Required Mesh Content
---------------------

A loaded mesh must be a valid unstructured grid with non-zero points and cells.

Required cell-data arrays
~~~~~~~~~~~~~~~~~~~~~~~~~

The following cell-data arrays must exist:

- TYPE
- MATERIAL
- CONSTANT
- COORDINATE_SYSTEM

If any are missing, validation fails.

Required point-data arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~

The array ``angles`` must exist in point data.

Required initialized arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The arrays a model needs follow from its element types (``src/resources/element_types.json``, read by
``pytec3d.element_library``). ``ModelEdit.update_grid`` (called when a model is read, when its model data is
read, and before a solve) adds the missing ones:

- TYPE, MATERIAL, CONSTANT (elem, default 1), COORDINATE_SYSTEM (elem, default 0), CP_0 (node, default 0);
- cell_centroids and cell_volume, when not computed yet;
- for each DOF field F: F (node) and its Dirichlet values _F (node, NaN where free); 3 components for a vector
  field, 1 for a scalar;
- for each nodal field: F (node); for each element field: F (elem); the same name can be both;
- the nodal loads (``mn.`` parameters, node) and element loads (``me.`` parameters, elem).

With the current library: A, _A, V, _V, VOLT, _VOLT, EMF, _EMF, B, H (node and elem), Jt (elem), AMPS, CURR
(node), Js (elem). The edge (NEDELEC_ONE) unknowns of A read their Dirichlet values from component 2 of _A on
the midside nodes; DOF couplings name them A.

Coordinate System
-----------------------

The coordinate system is used to define a transormation of coordinate from a global
cartesian coordinate system to a user-defined one, to help with many mathematical operations during analysis setup.
E.g. in tokamak simulations, the user can define a cylindrical coordinate system to simplify the definition of loads and boundary conditions.

The user can define a new coordinate system in the ``coordinate_systems`` dictionary of the model object.
To define a new coordinate system, the user must create a new entry in the ``coordinate_systems`` dictionary,
by instantiating a ``CoordinateSystem`` object with the appropriate type, orientation and position (see next).
A coordinate system is identified by a unique ID in the dictionary, which must be an integer greater than 0.
By default, coordinate system 0 is assumed to be the global cartesian coordinate system, and it is not present in the dictionary.



Elemental vector/tensor fields associated to the element are transformed from the global coordinate system
to the local coordinate system associated to the element defined by the user.
This includes:
- Material properties: mu, rho, M
- Elemental loads: Js
- Real constants: T

A coordinate system can be:
- Cartesian (Assumed as the defualt, if an element has ``COORDINATE_SYSTEM = 0`` it is assumed to be in the global Cartesian coordinate system)

- COORDINATE_SYSTEM values in mesh cells must be ``>= 0``.
- If any cell uses ``COORDINATE_SYSTEM > 0``:
	- ``model.coordinate_systems`` must exist.
	- It must be a dictionary.
	- Every positive ID used in mesh must be present in that dictionary.

Solver Settings Rules
---------------------

``solver_settings`` must satisfy:

- ``ksp_type``: string and one of supported PETSc KSP names.
- ``pc_type``: string and one of supported PETSc PC names.
- ``initial_linear_tolerance``: finite and ``> 0``.
- ``minimum_linear_tolerance``: finite and ``> 0`` and
	``<= initial_linear_tolerance``.
- ``max_linear_iterations``: integer ``>= 1``.

DataFrame Rules
---------------

General rule for all main model tables:

- ``element_types``, ``materials``, ``constants``, and ``dof_couplings``
	cannot be empty.
- Required columns must exist.
- Required columns cannot contain empty/NaN values.
- String columns cannot contain blank strings.

Element Types
~~~~~~~~~~~~~

The element type is used to associate the physics to a determined element.
A type represents the finite elements formulation used to build
the element matrix and the relative Degrees-Of-Freedom (DOFs, the electric/magnetic potentials of the
electromagnetic field in Maxwell equations) will be associated, as well as
the physical region it represents within the simulation, the loads and boundary conditions that can be applied on it.
All formulations are defined for both static and transient analysis.

The user can define the element types in the ``element_types`` table of the model object.
To define an element type, the user must create a row in the ``materials`` table
with the appropriate values for each property. An element type is identified by a unique ``TYPE_ID``, which
must be an integer greater than 0.

In the column ``ELEMENT_TYPE``, the user specifies the
the finite element formulation associated with that element type. The available options are 0, 1, and 2:

- 0: 2D A-V nodal formulation.
- 1: 3D A-V nodal formulation.
- 2: 3D A-V edge formulation.

.. 2D A-V nodal formulation (``ELEMENT_TYPE = 0``)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Use this formulation for 2D cells. For element types used by 2D cells,
.. ``ELEMENT_TYPE`` must be 0.

.. 3D A-V nodal formulation (``ELEMENT_TYPE = 1``)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Use this formulation for 3D cells with nodal unknowns. For element types used by
.. 3D cells, ``ELEMENT_TYPE`` can be 1 (or 2 for the edge-based formulation).

.. 3D A-V edge-based formulation (``ELEMENT_TYPE = 2``)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Use this formulation for 3D cells with edge-based unknowns. For element types
.. used by 3D cells, ``ELEMENT_TYPE`` can be 2 (or 1 for the nodal formulation).

Clearly, the user must assign appropriately a finite element formulation depending
on the dimensionality of the problem. The theoretical details of each formulation
are reported in the theory section.

To each element type, beyond the formulation, there are five options
(OPTION1 to OPTION5) that define specific characteristics or behaviors of the element:

- OPTION1: it is the physical region represented by the element in the problem. The available options are:
  - 0: Non-conducting air
  - 1: Soft-magnetic
  - 2: Hard magnetic
  - 3: Stranded conductor
  - 4: Short-circuit conductor
  - 5: Solid conductor
  - 6: Stranded coil
- OPTION2: defines the finite elements approximation order of the element. The available options are:
  - 1: First order: linear approximation
  - 2: Second order: quadratic approximation
- OPTION3: 2D behavior mode. The available options are:
  - 0: Planar
  - 1: Axisymmetric
- LAW (optional): the nonlinear laws of source ``me`` selected on the cells of this TYPE (see
  `Nonlinear laws`_). The formulation with time-integrated electric potentials is an element type of its own
  (``*_timint``), not an option.

Summary of required columns in ``element_types`` table:

- TYPE_ID, ELEMENT_TYPE, OPTION1, OPTION2, OPTION3

Summary of input rules:

- TYPE_ID must be ``> 0`` and unique.
- ELEMENT_TYPE must be in ``{0, 1, 2}``.
- OPTION1 must be in ``[0, 6]``.
- OPTION2 must be in ``{1, 2}``.
- OPTION3 must be in ``{0, 1}``.

Mesh compatibility rules:

- Every mesh TYPE value must exist in TYPE_ID (all elements must have a type assigned).
- For a TYPE_ID used by 2D cells, ELEMENT_TYPE must be 0 (planar or axisymmetric formulation).
- For a TYPE_ID used by 3D cells, ELEMENT_TYPE must be 1 or 2 (3D formulations, nodal or edge).
- For a TYPE_ID used by 3D cells, OPTION3 must be 0 (axisymmetric is undefined for 3D analysis).

Materials
~~~~~~~~~

The user can define the material properties sets in the ``materials`` table of the model object.
To define a set of material properties, the user must create a row in the ``materials`` table
with the appropriate values for each property. A set of material properties is
identified by a unique ``MATERIAL_ID``, which must be an integer ``> 0`` and unique.

The material properties define the electromagnetic material properties:
- magnetic permeability: MUR_X, MUR_Y, MUR_Z. They represent the three diagonal
     components of the magnetic permeability tensor in the X, Y, and Z directions.
- electrical resistivity: RHO_X, RHO_Y, RHO_Z. They represent the three diagonal
     components of the electrical resistivity tensor in the X, Y, and Z directions.
- magnetization: M_X, M_Y, M_Z. They represent the three components of the magnetization vector in the X, Y, and Z directions.

Not all the material properties are applicable to all element types
(e.g. a non-conducting element does not have electrical resistivity).
If a material property is not applicable to a certain element type, it will be ignored.
Details on which materials properties are applicable to which element types can be found in the individual element type documentation pages.

Summary of required columns in ``materials`` table:

- MATERIAL_ID, MUR_X, MUR_Y, MUR_Z, RHO_X, RHO_Y, RHO_Z, M_X, M_Y, M_Z

Summary Rules:

- MATERIAL_ID must be ``> 0`` and unique.
- Each property must be one of the material properties of the element types, with their shape (number,
  3-vector or 3x3 tensor), and finite.
- Every mesh MATERIAL value must exist in MATERIAL_ID.

Which properties an element needs depends on its element type and domain: the solver checks that they are
defined when the model is loaded.

A material may also select nonlinear laws of source ``mat`` with its ``law`` key, e.g. ``"law": "ferromagnetic"``
(see `Nonlinear laws`_).

Nonlinear laws
~~~~~~~~~~~~~~

The nonlinear laws are defined in ``src/resources/laws.json`` and compiled into the solver. A law replaces one or
more coefficients of the element types (e.g. ``ferromagnetic`` replaces the reluctivity ``nu`` and ``nu_gauge``)
by values computed at each quadrature point from a field (``B``) and, optionally, a curve of the model data.

A law is selected by the ``law`` key of the rows of its source:

- ``mat``: the ``materials`` rows (the law follows the material);
- ``rc``: the ``constants`` rows;
- ``me``: the ``element_types`` rows (by TYPE).

``law`` is a law name (or index) or a list of them. On a cell, the laws selected by its TYPE, MATERIAL and CONSTANT
rows apply wherever its element type uses a coefficient they replace; elsewhere the coefficients are the linear
material properties and real constants. A model with at least one such cell is solved with Newton's method. The
curve of a law is keyed by the id of the row that selects it, e.g. ``bh_curves[<MATERIAL id>]`` for
``ferromagnetic``, and every row selecting a law with a curve must have one.

Constants
~~~~~~~~~

Required columns:

- CONSTANTS_ID, THK, CS, TURNS, VOL, RAD, T_X, T_Y, T_Z,
	R, S, V_X, V_Y, V_Z, OMEGA_X, OMEGA_Y, OMEGA_Z, LOC_X, LOC_Y, LOC_Z

Rules:

- CONSTANTS_ID must be ``> 0`` and unique.
- All numeric constant fields must be finite.
- Every mesh CONSTANT value must exist in CONSTANTS_ID.

Which real constants an element needs depends on its element type and domain: the solver checks that they
are defined when the model is loaded.

dof_couplings
~~~~~~~~~~~~~

Required columns:

- CP_ID, VARIABLE

Rules:

- CP_ID must be ``> 0``.
- VARIABLE must be one of: A, Az, VOLT, V, EMF.

CP Array to Coupling Table Consistency
--------------------------------------

Point-data arrays matching ``CP_<xx>`` (for integer ``xx``) are scanned.

- All finite values ``> 0`` found in these CP arrays must exist in
	``dof_couplings.CP_ID``.
- If a CP ID is used in mesh but not declared in ``dof_couplings``,
	validation fails.

Curves of the Laws
------------------

The curve of a law (e.g. ``bh_curves`` for ``ferromagnetic``) is a section of the model data keyed by the id of
the row selecting the law. Each curve is an object with the two arrays the law names (``B`` and ``H`` for
``bh_curves``):

- same length ``>= 2``, finite values;
- the x array (``B``) strictly increasing.

Beyond its last point a curve is extended along its last segment: give the table up to saturation.

Pre-Solve Validation
--------------------

Pre-solve validation runs before solver execution.

Warnings
~~~~~~~~

1. If all values in Js, AMPS, and CURR are zero:
	 warning that no loads are being applied.
2. If all values in the Dirichlet arrays (_A, _V, _VOLT, _EMF, ...) are NaN:
	 warning that no boundary conditions are applied.

Errors
~~~~~~

``time_step`` must satisfy:

- ``time > 0``
- ``0 <= theta <= 1`` (first-order formulations)
- ``0 < beta <= 0.5`` and ``0 <= gamma <= 1`` (Newmark, second-order formulations; default 0.25 and 0.5)

If not, pre-solve validation raises an error and solve is blocked.

Recommended Input Workflow
--------------------------

1. Load mesh.
2. Ensure required VTK arrays exist (or run array initialization workflow).
3. Fill ``element_types``, ``materials``, ``constants``, and
	 ``dof_couplings`` with no empty/NaN fields.
4. Verify mesh IDs (TYPE/MATERIAL/CONSTANT and CP IDs) are declared in tables.
5. Configure ``solver_settings`` and ``time_step``.
6. Run validation.
7. Solve.

Troubleshooting Tips
--------------------

- Missing-array errors usually mean the mesh was loaded without initialization
	or arrays were removed/renamed by an external tool.
- CP consistency errors indicate couplings were painted into mesh arrays but
	not declared in ``dof_couplings``.
- No-load and no-BC warnings do not block solve, but they often indicate an
	unintended setup that produces trivial results.
