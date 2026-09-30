User Manual
===========

This manual documents the practical input rules enforced by model validation.
Use it as a checklist before loading a model into the solver and before solve.

Validation Philosophy
---------------------

APEC-M validation has two stages:

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

All arrays listed in ``VTK_ARRAYS_TO_INIT`` must exist in the expected
location (node or elem), with expected tuple count and number of components.

Current required arrays are:

- angles (node, 3)
- A (node, 3)
- Az (node, 1)
- VOLT (node, 1)
- V (node, 1)
- EMF (node, 1)
- _A (node, 3)
- _Az (node, 1)
- _VOLT (node, 1)
- _V (node, 1)
- _EMF (node, 1)
- B (node, 3)
- H (node, 3)
- Be (elem, 3)
- He (elem, 3)
- Jt (elem, 3)
- AMPS (node, 1)
- CURR (node, 1)
- Js (elem, 3)

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
- OPTION4: defines if a formulation with time integrated electric potentials is used. The available options are:
  - 0: Formulation with non time integrated electric scalar potentials
  - 1: Formulation with time integrated electric scalar potentials
- OPTION5: reserved for future use.

Summary of required columns in ``element_types`` table:

- TYPE_ID, ELEMENT_TYPE, OPTION1, OPTION2, OPTION3, OPTION4, OPTION5

Summary of input rules:

- TYPE_ID must be ``> 0`` and unique.
- ELEMENT_TYPE must be in ``{0, 1, 2}``.
- OPTION1 must be in ``[0, 6]``.
- OPTION2 must be in ``{1, 2}``.
- OPTION3 must be in ``{0, 1}``.
- OPTION4 must be in ``{0, 1}``.

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
- MUR_* values must be finite and ``> 0``.
- RHO_* and M_* values must be finite and ``>= 0``.
- Every mesh MATERIAL value must exist in MATERIAL_ID.

Constants
~~~~~~~~~

Required columns:

- CONSTANTS_ID, THK, CS, TURNS, VOL, RAD, T_X, T_Y, T_Z,
	R, S, V_X, V_Y, V_Z, OMEGA_X, OMEGA_Y, OMEGA_Z, LOC_X, LOC_Y, LOC_Z

Rules:

- CONSTANTS_ID must be ``> 0`` and unique.
- All numeric constant fields must be finite.
- Every mesh CONSTANT value must exist in CONSTANTS_ID.

Context-specific rules (based on element type/option):

- For 2D type-0 cells:
	- THK must be ``> 0``.
	- If OPTION3=1 then THK must be ``<= 1``.
- For 2D type-0 with OPTION1=6:
	- CS>0, TURNS>0, RAD>0, T_Z in {-1, 1}, R>0, S>=1.
- For 3D type-1 or type-2 with OPTION1=6:
	- CS>0, TURNS>0, VOL>0, R>0, S>=1.
	- T_X, T_Y, T_Z each in [-1, 1].
	- (T_X, T_Y, T_Z) must be unit norm (within tolerance).

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

BH Curve Rules
--------------

``bh_curves`` must be a dictionary from material ID to a numeric ``(N,2)``
array of ``[H, B]`` pairs.

Rules for each curve:

- Key must be integer material ID > 0.
- Curve must be convertible to float array.
- Shape must be ``(N, 2)`` with ``N >= 2``.
- All values must be finite.
- First point must be (0, 0).
- H and B values must be non-negative and strictly increasing.

Additional solver-consistency notes may be emitted as warnings.

Pre-Solve Validation
--------------------

Pre-solve validation runs before solver execution.

Warnings
~~~~~~~~

1. If all values in Js, AMPS, and CURR are zero:
	 warning that no loads are being applied.
2. If all values in _A, _Az, _VOLT, _V, and _EMF are NaN:
	 warning that no boundary conditions are applied.

Errors
~~~~~~

``time_step`` must satisfy:

- ``time > 0``
- ``0 <= theta <= 1``

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
