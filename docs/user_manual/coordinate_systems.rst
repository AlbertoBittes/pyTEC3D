Coordinate Systems
==================

Coordinate systems (:mod:`pytec3d.coordinate_system`) are used in three places:

* **selections and geometry** in Python: ``select_pos(..., coordinate_system=cs)``, ``periodic_boundaries(cs, ...)``,
  ``model_edit.to_local_coordinates`` / ``map_vector_field_to_*_components``;
* **element coordinate systems**: the cell array ``COORDINATE_SYSTEM`` refers a cell to an entry of
  ``data["coordinate_systems"]``; the vector and tensor parameters of the cell (``mur``, ``rho``, ``M``, ``T``, ``Js``)
  are then given in that system;
* **nodal frames** (``node_angles``), a separate mechanism for the Dirichlet values and the couplings of each node
  (:doc:`boundary_conditions`).

Types
-----

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Class / ``type``
     - Coordinates
     - Default orientation
   * - :class:`~pytec3d.coordinate_system.CartesianCoordinateSystem` / ``"cartesian"``
     - x, y, z
     - the global axes
   * - :class:`~pytec3d.coordinate_system.CylindricalCoordinateSystem` / ``"cylindrical"``
     - r, theta [deg], z
     - local z axis = global -y
   * - :class:`~pytec3d.coordinate_system.XAxisCylindricalCoordinateSystem` / ``"x_axis_cylindrical"``
     - r, theta, z
     - axis = global x (fixed orientation)
   * - :class:`~pytec3d.coordinate_system.YAxisCylindricalCoordinateSystem` / ``"y_axis_cylindrical"``
     - r, theta, z
     - axis = global y: theta from +x towards -z, z = global y (the axis of the 2D axisymmetric models)
   * - :class:`~pytec3d.coordinate_system.ZAxisCylindricalCoordinateSystem` / ``"z_axis_cylindrical"``
     - r, theta, z
     - axis = global z: theta from +x towards +y
   * - :class:`~pytec3d.coordinate_system.SphericalCoordinateSystem` / ``"spherical"``
     - r, theta, phi [deg]
     - theta from the local z axis, phi around it

A system is an origin and a rotation (the matrix whose columns are the local x, y, z axes in global coordinates).
The orientation is given by one of:

.. code-block:: python

   from pytec3d import CylindricalCoordinateSystem, YAxisCylindricalCoordinateSystem

   YAxisCylindricalCoordinateSystem()                                   # default orientation, origin (0, 0, 0)
   YAxisCylindricalCoordinateSystem(origin=[0.0, 0.5, 0.0])
   CylindricalCoordinateSystem.from_angles([1, 2, 3], [0, 0, 90])       # thxy, thyz, thzx (APDL LOCAL)
   CylindricalCoordinateSystem.from_points([0, 0, 0], nxax=[1, 0, 0], nxypl=[0, 1, 0])
   CylindricalCoordinateSystem(local_to_global_cart=[[...], [...], [...]])

Methods: ``to_local(points)`` / ``to_global(points)`` for coordinates, ``vector_components_to_local(vectors,
points)`` / ``vector_components_to_global`` for vector components at points, ``local_to_global_cart`` (the rotation),
``axis_direction()``.

Element coordinate systems
--------------------------

.. code-block:: python

   solver.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())   # data["coordinate_systems"]["1"]
   selection.select(mesh, "MATERIAL", 23, location="elem")
   selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")

In JSON (``matrix``, ``angles`` or ``points``; the ``*_axis_cylindrical`` types accept only ``matrix``):

.. code-block:: json

   "coordinate_systems": {"1": {"type": "y_axis_cylindrical", "origin": [0, 0, 0]},
                          "2": {"type": "cylindrical", "origin": [1, 2, 3], "angles": [0, 0, 90]},
                          "3": {"type": "cartesian", "origin": [0, 0, 0],
                                "points": {"nxax": [1, 1, 0], "nxypl": [0, 1, 0]}}}

The solver evaluates the rotation of the system **at the centre of each cell**: in a cylindrical system the local
basis :math:`(\mathbf{e}_r,\mathbf{e}_\theta,\mathbf{e}_z)` of the cell is used. So:

* ``T = (0, 1, 0)`` is the toroidal direction in every cell of a coil about the axis of the system;
* ``Js = (0, J, 0)`` is a toroidal current density;
* an orthotropic ``mur`` or ``rho`` tensor is given in the local axes (e.g. a laminated core).

Use one system per curved part with its own centre (e.g. the two arcs of a racetrack coil).
``COORDINATE_SYSTEM = 0`` (default) is the global Cartesian system. Every positive id used by the cells must exist in
``data["coordinate_systems"]``.
