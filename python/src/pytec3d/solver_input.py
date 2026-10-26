"""What the C++ solver receives: typed views of the mesh and of the arrays it needs, the rotations, the model data.

The solver (``pytec3d_solver``) does not know VTK nor this package. :func:`model_input` hands it the memory of the VTK
arrays of the model as NumPy *views* (no copy), so the solver reads the loads and the Dirichlet values straight from the
mesh at every solve, and writes the solution straight into its field arrays::

    {
        "mesh": {"points": (n_points, 3) float64,
                 "offsets": (n_cells + 1,) int64, "connectivity": (offsets[-1],) int64,
                 "cell_types": (n_cells,) uint8},                # VTK cell type numbers and node order
        "arrays": {                                             # model_arrays(): used until the session closes
            "ids": {"TYPE", "MATERIAL", "CONSTANT", "COORDINATE_SYSTEM": (n_cells,) int32},
            "couplings": {"CP_<k>": (n_points,) int32},
            "point_fields": {name: (n_points, 1 or 3) float64},  # DOF fields F and _F, nodal fields, nodal loads
            "cell_fields": {name: (n_cells, 1 or 3) float64},    # element fields, element loads
        },
        "element_rotations": (n_cells, 3, 3) float64 or None,   # COORDINATE_SYSTEM of each cell, at its vertex average
        "node_rotations": (n_points, 9) float64 or None,        # the node_angles frame of each point
        "model_data": dict,                                     # Model.data, plain JSON types
    }

Only the arrays the solver needs are handed over: the fields are those of the element types of the model data
(:func:`pytec3d.element_library.grid_arrays`, the same list the solver derives from the model data and checks). Every
array has the type of its role; :func:`pytec3d.mesh.update_grid` (called when the model is loaded) stores the mesh
with these types, so the views need no conversion.

The rotations are computed here, from the coordinate systems of the model data and the ``node_angles`` array: they
are the only arrays that are not views of the mesh. The solver copies them, with the mesh, when the model is loaded
and keeps no reference to them: they are freed with the model input once the load returns. The views stay valid while
the VTK arrays they come from are alive and not resized: :func:`model_arrays` is called again before every solve
(:meth:`pytec3d.Solver.solve`), so that an array replaced in the mesh after the load (e.g. by
:func:`pytec3d.mesh.add_array` with ``replace=True``) is the one the solver uses.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

from . import element_library
from . import mesh as meshes
from ._logging import get_logger
from ._typing import Location
from .coordinate_system import rotation_coefficients_from_angles

if TYPE_CHECKING:
    from .model import Model

LOGGER = get_logger("pytec3d.solver_input")

#: Corner nodes (listed first in the cell) of each VTK cell type the solver reads.
CELL_VERTICES: dict[int, int] = {
    3: 2, 21: 2,                # line, quadratic edge
    5: 3, 22: 3,                # triangle, quadratic triangle
    9: 4, 23: 4, 28: 4,         # quad, quadratic quad, biquadratic quad
    10: 4, 24: 4,               # tetra, quadratic tetra
    12: 8, 25: 8, 29: 8,        # hexahedron, quadratic hexahedron, triquadratic hexahedron
    13: 6, 26: 6, 32: 6,        # wedge, quadratic wedge, biquadratic-quadratic wedge
    14: 5, 27: 5,               # pyramid, quadratic pyramid
}

#: NumPy types of the views, by role.
REAL = np.dtype(np.float64)
ID = np.dtype(np.int32)
INDEX = np.dtype(np.int64)
CELL_TYPE = np.dtype(np.uint8)


def _typed_view(array: Any, dtype: np.dtype, what: str) -> np.ndarray:
    """NumPy view of a VTK data array, which must have ``dtype`` and be contiguous.

    Raises:
        ValueError: If it has another type (the mesh was not stored with the types of the solver) or no contiguous
            memory.
    """
    try:
        values = vtknp.vtk_to_numpy(array)
    except Exception as exc:  # noqa: BLE001 - bit arrays and other non-AOS arrays
        raise ValueError(f"{what} has no contiguous memory the solver can view.") from exc
    if values.dtype != dtype or not values.flags.c_contiguous:
        raise ValueError(
            f"{what} is {values.dtype}{'' if values.flags.c_contiguous else ', not contiguous'}; the solver views "
            f"it as {dtype}: call model.update_grid() first (pytec3d.Solver.load_model_into_solver does)."
        )
    return values


def _cell_types_array(mesh: vtkUnstructuredGrid) -> Any:
    """The VTK cell type array of the mesh (``GetCellTypes()`` from VTK 9.6, ``GetCellTypesArray()`` before)."""
    try:
        array = mesh.GetCellTypes()
    except TypeError:  # VTK < 9.6: GetCellTypes(vtkCellTypes) fills its argument
        array = None
    return array if array is not None and hasattr(array, "GetNumberOfTuples") else mesh.GetCellTypesArray()


def mesh_views(mesh: vtkUnstructuredGrid) -> dict[str, np.ndarray]:
    """Views of the points (float64), the cell offsets and connectivity (int64) and the cell types (uint8).

    Raises:
        ValueError: If the mesh has no points, or one of them does not have its type (see :func:`_typed_view`).
    """
    if mesh.GetPoints() is None:
        raise ValueError("The mesh has no points.")
    cells = mesh.GetCells()
    return {
        "points": _typed_view(mesh.GetPoints().GetData(), REAL, "The points of the mesh"),
        "offsets": _typed_view(cells.GetOffsetsArray(), INDEX, "The cell offsets of the mesh"),
        "connectivity": _typed_view(cells.GetConnectivityArray(), INDEX, "The cell connectivity of the mesh"),
        "cell_types": _typed_view(_cell_types_array(mesh), CELL_TYPE, "The cell types of the mesh"),
    }


def _array_view(mesh: vtkUnstructuredGrid, name: str, location: Location, dtype: np.dtype) -> np.ndarray:
    where = "point" if location == "node" else "cell"
    return _typed_view(meshes.get_array(mesh, name, location), dtype, f"The {where} array '{name}'")


def model_arrays(model: Model) -> dict[str, dict[str, np.ndarray]]:
    """Views of the arrays the solver needs (see the module documentation): the ids, the couplings, and the fields of
    the element types of the model data.

    Raises:
        MissingArrayError: If one of them is missing (:meth:`pytec3d.Model.update_grid` adds them).
        ValueError: If one of them does not have its type.
    """
    mesh = model.mesh
    point_data = mesh.GetPointData()
    fields = element_library.grid_arrays(model.element_type_names() or None)
    return {
        "ids": {name: _array_view(mesh, name, "elem", ID) for name in meshes.MODEL_CELL_ARRAYS},
        "couplings": {
            name: _array_view(mesh, name, "node", ID)
            for name in (point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays()))
            if name and meshes.COUPLING_ARRAY.match(name)
        },
        "point_fields": {a.name: _array_view(mesh, a.name, "node", REAL) for a in fields if a.location == "node"},
        "cell_fields": {a.name: _array_view(mesh, a.name, "elem", REAL) for a in fields if a.location == "elem"},
    }


def vertex_averages(mesh_arrays: dict[str, np.ndarray], cells: np.ndarray) -> np.ndarray:
    """Mean of the corner points of each of ``cells`` (the element centre the solver uses), shape (len(cells), 3).

    Raises:
        ValueError: If a cell type is not one of :data:`CELL_VERTICES`.
    """
    points = mesh_arrays["points"]
    offsets = mesh_arrays["offsets"]
    connectivity = mesh_arrays["connectivity"]
    types = mesh_arrays["cell_types"][cells]
    table = np.zeros(256, dtype=np.int64)
    for cell_type, n_vertices in CELL_VERTICES.items():
        table[cell_type] = n_vertices
    vertices = table[types]
    if np.any(vertices == 0):
        unknown = sorted(set(types[vertices == 0].tolist()))
        raise ValueError(f"VTK cell types {unknown} are not supported by the solver.")

    centres = np.empty((cells.size, 3), dtype=float)
    for n_vertices in np.unique(vertices).tolist():
        rows = np.flatnonzero(vertices == n_vertices)
        first = offsets[cells[rows]]
        total = points[connectivity[first]].copy()
        for corner in range(1, n_vertices):  # summed in node order, as libMesh does
            total += points[connectivity[first + corner]]
        centres[rows] = total / n_vertices
    return centres


def element_rotations(model: Model, mesh_arrays: dict[str, np.ndarray]) -> np.ndarray | None:
    """Local-to-global rotation of the COORDINATE_SYSTEM of each cell at its vertex average, shape (n_cells, 3, 3).

    None when every cell is in the global system (COORDINATE_SYSTEM 0).

    Raises:
        KeyError: If a cell uses a coordinate system the model data does not define.
    """
    ids = meshes.array_values(model.mesh, "COORDINATE_SYSTEM", "elem")
    cells = np.flatnonzero(ids != 0)
    if cells.size == 0:
        return None
    rotations = np.empty((ids.shape[0], 3, 3), dtype=float)
    rotations[:] = np.eye(3)
    centres = vertex_averages(mesh_arrays, cells)
    cell_ids = ids[cells]
    for cs_id in np.unique(cell_ids).tolist():
        rows = np.flatnonzero(cell_ids == cs_id)
        rotations[cells[rows]] = model.coordinate_system(int(cs_id)).local_to_global_matrices(centres[rows])
    return rotations


def node_rotations(mesh: vtkUnstructuredGrid) -> np.ndarray | None:
    """Local-to-global rotation of the ``node_angles`` frame of each point, flattened, shape (n_points, 9).

    None when the mesh has no ``node_angles`` array or all its angles are zero.
    """
    if mesh.GetPointData().GetArray(meshes.NODE_ANGLES_ARRAY) is None:
        return None
    angles = np.asarray(meshes.array_values(mesh, meshes.NODE_ANGLES_ARRAY, "node"), dtype=float).reshape(-1, 3)
    if not np.any(angles):
        return None
    return rotation_coefficients_from_angles(angles)


def model_input(model: Model) -> dict[str, Any]:
    """The model input of the solver (see the module documentation): views of the mesh and of the arrays it needs.

    Raises:
        SolverStateError: If the model has no mesh.
        MissingArrayError: If an array the solver needs is missing.
        ValueError: If an array does not have its type, or the model data cannot be written as JSON.
        KeyError: If a cell uses a coordinate system the model data does not define.
    """
    mesh = model.mesh
    mesh_arrays = mesh_views(mesh)
    return {
        "mesh": mesh_arrays,
        "arrays": model_arrays(model),
        "element_rotations": element_rotations(model, mesh_arrays),
        "node_rotations": node_rotations(mesh),
        # Plain JSON types, finite numbers (as checked by Model.data_json).
        "model_data": json.loads(model.data_json()),
    }
