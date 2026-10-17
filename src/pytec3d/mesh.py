"""The VTK unstructured grid of a model: reading and writing it, its arrays and its geometry.

Arrays live on the points (``location="node"``) or on the cells (``location="elem"``). The functions here return
NumPy *views* of VTK arrays where they can, so writing into the result changes the mesh; call ``Modified()`` on the
VTK array afterwards (the functions of the package do).
"""

from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import VTK_DOUBLE, VTK_INT, vtkDataArray
from vtkmodules.vtkCommonDataModel import vtkCellTypes, vtkDataSetAttributes, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkCellCenters
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter
from vtkmodules.vtkIOXML import (
    vtkXMLPUnstructuredGridReader,
    vtkXMLUnstructuredGridReader,
    vtkXMLUnstructuredGridWriter,
)

from . import element_library
from ._logging import get_logger
from ._typing import LOCATIONS, Location, PathType
from .errors import MissingArrayError

LOGGER = get_logger("pytec3d.mesh")

#: Cell arrays every model has, with their default value: the TYPE, MATERIAL, CONSTANT and COORDINATE_SYSTEM ids.
MODEL_CELL_ARRAYS: dict[str, int] = {"TYPE": 1, "MATERIAL": 1, "CONSTANT": 1, "COORDINATE_SYSTEM": 0}
#: Point arrays every model has, with their default value: the first DOF coupling ids.
MODEL_POINT_ARRAYS: dict[str, int] = {"CP_0": 0}
#: Point array of the orientation angles (thxy, thyz, thzx, degrees) of each node: the local frame of its DOFs.
NODE_ANGLES_ARRAY = "node_angles"
#: Cell arrays of the cell centroids and of the cell sizes (area of 2D cells, volume of 3D cells).
CENTROIDS_ARRAY = "cell_centroids"
SIZES_ARRAY = "cell_volume"


# ---------------------------------------------------------------------------------------------------- checks


def check_grid(mesh: Any, name: str = "mesh") -> vtkUnstructuredGrid:
    """Return ``mesh`` if it is a ``vtkUnstructuredGrid``.

    Raises:
        TypeError: If it is not.
    """
    if not isinstance(mesh, vtkUnstructuredGrid):
        raise TypeError(f"{name} must be a vtkUnstructuredGrid, got {type(mesh).__name__}.")
    return mesh


def check_location(location: str) -> Location:
    """Return ``location`` if it is ``"node"`` or ``"elem"``.

    Raises:
        ValueError: If it is not.
    """
    if location not in LOCATIONS:
        raise ValueError(f"location must be 'node' or 'elem', got {location!r}.")
    return location  # type: ignore[return-value]


# ---------------------------------------------------------------------------------------------------- access


def attributes(mesh: vtkUnstructuredGrid, location: Location) -> vtkDataSetAttributes:
    """The point data (``"node"``) or the cell data (``"elem"``) of the mesh."""
    return mesh.GetPointData() if check_location(location) == "node" else mesh.GetCellData()


def size(mesh: vtkUnstructuredGrid, location: Location) -> int:
    """Number of points (``"node"``) or cells (``"elem"``) of the mesh."""
    return mesh.GetNumberOfPoints() if check_location(location) == "node" else mesh.GetNumberOfCells()


def get_array(mesh: vtkUnstructuredGrid, name: str, location: Location) -> vtkDataArray:
    """The VTK array ``name`` of the point or cell data.

    Raises:
        MissingArrayError: If the mesh has no such array.
    """
    array = attributes(mesh, location).GetArray(name)
    if array is None:
        raise MissingArrayError(name, location)
    return array


def array_values(mesh: vtkUnstructuredGrid, name: str, location: Location) -> np.ndarray:
    """NumPy view of the array ``name``: shape (N,) for one component, (N, components) otherwise.

    Raises:
        MissingArrayError: If the mesh has no such array.
    """
    return vtknp.vtk_to_numpy(get_array(mesh, name, location))


def points(mesh: vtkUnstructuredGrid) -> np.ndarray:
    """NumPy view of the point coordinates, shape (N, 3); empty when the mesh has no points."""
    vtk_points = mesh.GetPoints()
    return np.empty((0, 3)) if vtk_points is None else vtknp.vtk_to_numpy(vtk_points.GetData())


def cell_connectivity(mesh: vtkUnstructuredGrid) -> tuple[np.ndarray, np.ndarray]:
    """``(offsets, connectivity)`` of the cells: cell ``i`` has the points ``connectivity[offsets[i]:offsets[i+1]]``."""
    cells = mesh.GetCells()
    offsets = np.asarray(vtknp.vtk_to_numpy(cells.GetOffsetsArray()), dtype=np.int64)
    connectivity = np.asarray(vtknp.vtk_to_numpy(cells.GetConnectivityArray()), dtype=np.int64)
    return offsets, connectivity


def cell_types(mesh: vtkUnstructuredGrid) -> np.ndarray:
    """VTK cell type of each cell (``VTK_TRIANGLE``, ``VTK_QUADRATIC_HEXAHEDRON``...)."""
    return np.asarray(vtknp.vtk_to_numpy(mesh.GetCellTypesArray()), dtype=np.int64)


def cell_dimensions(mesh: vtkUnstructuredGrid) -> np.ndarray:
    """Topological dimension (0 to 3) of each cell."""
    types = cell_types(mesh)
    unique, inverse = np.unique(types, return_inverse=True)
    dimensions = np.array([vtkCellTypes.GetDimension(int(t)) for t in unique], dtype=np.int8)
    return dimensions[inverse]


# ---------------------------------------------------------------------------------------------------- arrays


def add_array(
    mesh: vtkUnstructuredGrid,
    name: str,
    location: Location = "node",
    components: int = 1,
    fill: float = 0,
    dtype: int | None = VTK_DOUBLE,
    values: np.ndarray | None = None,
    replace: bool = False,
) -> vtkDataArray:
    """Add an array to the point or cell data, unless it is already there.

    Args:
        mesh: The mesh.
        name: Array name.
        location: ``"node"`` (point data) or ``"elem"`` (cell data).
        components: Number of components (at least 1).
        fill: Value of every entry when ``values`` is not given.
        dtype: VTK type of the array (``VTK_DOUBLE``, ``VTK_INT``...); None to keep the NumPy type.
        values: Initial values, shape (N,) or (N, components).
        replace: Replace an existing array of that name instead of keeping it.

    Returns:
        The new array, or the existing one when ``replace`` is False.

    Raises:
        ValueError: If ``location`` is invalid or ``values`` has the wrong shape.
    """
    data = attributes(mesh, location)
    existing = data.GetArray(name)
    if existing is not None:
        if not replace:
            return existing
        data.RemoveArray(name)

    count = size(mesh, location)
    components = max(1, int(components))
    shape = (count,) if components == 1 else (count, components)
    if values is None:
        stored = np.full(shape, fill)
    else:
        stored = np.asarray(values)
        if stored.shape == (count, 1) and components == 1:
            stored = stored.reshape(count)
        if stored.shape != shape:
            raise ValueError(f"Values for array '{name}' must have shape {shape}, got {stored.shape}.")

    array = vtknp.numpy_to_vtk(stored, deep=True, array_type=dtype)
    array.SetName(name)
    array.SetNumberOfComponents(components)
    data.AddArray(array)
    array.Modified()
    return array


def remove_array(mesh: vtkUnstructuredGrid, name: str, location: Location | None = None) -> bool:
    """Remove an array from the point or cell data.

    Args:
        mesh: The mesh.
        name: Array name.
        location: Where to look; None removes it from the cell data, or else from the point data.

    Returns:
        True when an array was removed.
    """
    for where in LOCATIONS[::-1] if location is None else (check_location(location),):
        data = attributes(mesh, where)
        if data.HasArray(name):
            data.RemoveArray(name)
            return True
    LOGGER.debug("remove_array: no array '%s' to remove.", name)
    return False


# ---------------------------------------------------------------------------------------------------- geometry


def cell_centroids(mesh: vtkUnstructuredGrid) -> np.ndarray:
    """Centroid of each cell (``vtkCellCenters``), shape (cells, 3)."""
    centers = vtkCellCenters()
    centers.SetInputData(mesh)
    centers.Update()
    return np.asarray(vtknp.vtk_to_numpy(centers.GetOutput().GetPoints().GetData()))


def cell_sizes(mesh: vtkUnstructuredGrid) -> np.ndarray:
    """Area of each 2D cell and volume of each 3D cell (0 for lower dimensions), shape (cells,)."""
    n_cells = mesh.GetNumberOfCells()
    if n_cells == 0:
        return np.zeros(0)

    size_filter = vtkCellSizeFilter()
    size_filter.SetInputData(mesh)
    size_filter.SetComputeVertexCount(False)
    size_filter.SetComputeLength(False)
    size_filter.SetComputeArea(True)
    size_filter.SetComputeVolume(True)
    size_filter.SetComputeSum(False)
    size_filter.SetAreaArrayName("__area__")
    size_filter.SetVolumeArrayName("__volume__")
    size_filter.Update()
    cell_data = size_filter.GetOutput().GetCellData()
    areas = np.asarray(vtknp.vtk_to_numpy(cell_data.GetArray("__area__")), dtype=float).reshape(-1)
    volumes = np.asarray(vtknp.vtk_to_numpy(cell_data.GetArray("__volume__")), dtype=float).reshape(-1)

    dimensions = cell_dimensions(mesh)
    sizes = np.zeros(n_cells)
    sizes[dimensions >= 3] = volumes[dimensions >= 3]
    sizes[dimensions == 2] = areas[dimensions == 2]
    return sizes


def store_cell_centroids(mesh: vtkUnstructuredGrid, name: str = CENTROIDS_ARRAY) -> vtkDataArray:
    """Compute the cell centroids and store them in the cell array ``name`` (replacing it)."""
    return add_array(mesh, name, "elem", components=3, dtype=VTK_DOUBLE, values=cell_centroids(mesh), replace=True)


def store_cell_sizes(mesh: vtkUnstructuredGrid, name: str = SIZES_ARRAY) -> vtkDataArray:
    """Compute the cell areas / volumes and store them in the cell array ``name`` (replacing it)."""
    return add_array(mesh, name, "elem", components=1, dtype=VTK_DOUBLE, values=cell_sizes(mesh), replace=True)


# ---------------------------------------------------------------------------------------------------- model arrays


def update_grid(
    mesh: vtkUnstructuredGrid,
    element_types: Iterable[str | int] | None = None,
    reset: bool = False,
) -> None:
    """Give the mesh the arrays a model needs.

    * the integer cell arrays TYPE, MATERIAL, CONSTANT (default 1) and COORDINATE_SYSTEM (default 0), the integer
      point array CP_0 (default 0) and the point array ``node_angles`` (3 components, default 0), when missing;
    * the cell centroids and sizes, when not computed yet;
    * the field arrays of the element types (:func:`pytec3d.element_library.grid_arrays`): DOF fields and their
      Dirichlet arrays ``_F`` (NaN), nodal and element fields, nodal and element loads (0). Missing ones are added, and
      an existing one with another number of components or another type is replaced. With ``element_types`` given,
      the field arrays of the other element types of the library are removed; other arrays are never touched.

    Args:
        mesh: The mesh to update in place.
        element_types: The element types (names or indices) of the model; all of the library when None.
        reset: Also set the field arrays of the element types back to their default values.
    """
    check_grid(mesh)
    for name, default in MODEL_CELL_ARRAYS.items():
        add_array(mesh, name, "elem", fill=default, dtype=VTK_INT)
    for name, default in MODEL_POINT_ARRAYS.items():
        add_array(mesh, name, "node", fill=default, dtype=VTK_INT)
    add_array(mesh, NODE_ANGLES_ARRAY, "node", components=3, fill=0.0, dtype=VTK_DOUBLE)

    cell_data = mesh.GetCellData()
    if cell_data.GetArray(CENTROIDS_ARRAY) is None:
        store_cell_centroids(mesh)
    if cell_data.GetArray(SIZES_ARRAY) is None:
        store_cell_sizes(mesh)

    needed = element_library.grid_arrays(element_types)
    added: list[str] = []
    for grid_array in needed:
        existing = attributes(mesh, grid_array.location).GetArray(grid_array.name)
        mismatched = existing is not None and (
            existing.GetNumberOfComponents() != grid_array.components or existing.GetDataType() != VTK_DOUBLE
        )
        if existing is None or mismatched or reset:
            if existing is None or mismatched:
                added.append(f"{grid_array.name} ({'point' if grid_array.location == 'node' else 'cell'})")
            add_array(
                mesh,
                grid_array.name,
                grid_array.location,
                components=grid_array.components,
                fill=grid_array.fill,
                dtype=VTK_DOUBLE,
                replace=True,
            )

    if element_types is not None:
        kept = {(a.location, a.name) for a in needed}
        removed = [a for a in element_library.grid_arrays() if (a.location, a.name) not in kept]
        removed = [a for a in removed if attributes(mesh, a.location).HasArray(a.name)]
        for grid_array in removed:
            attributes(mesh, grid_array.location).RemoveArray(grid_array.name)
        if removed:
            LOGGER.info(
                "Removed the field arrays not used by the element types: %s.",
                ", ".join(f"{a.name} ({'point' if a.location == 'node' else 'cell'})" for a in removed),
            )
    if added:
        LOGGER.info("Added the field arrays: %s.", ", ".join(added))
    if reset:
        LOGGER.info("Field arrays reset to their defaults.")


# ---------------------------------------------------------------------------------------------------- files


def read_mesh(path: PathType, reset_arrays: bool = False) -> vtkUnstructuredGrid:
    """Read an unstructured grid from a ``.vtu`` or ``.pvtu`` file and add the arrays a model needs.

    Args:
        path: File to read.
        reset_arrays: Set the arrays of the element types back to their default values (:func:`update_grid`).

    Returns:
        The mesh.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the extension is not ``.vtu`` / ``.pvtu`` or the mesh has no points or no cells.
    """
    file_path = os.fspath(path)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Mesh file '{file_path}' not found.")
    extension = os.path.splitext(file_path)[1].lower()
    readers = {".vtu": vtkXMLUnstructuredGridReader, ".pvtu": vtkXMLPUnstructuredGridReader}
    if extension not in readers:
        raise ValueError(f"Unsupported mesh extension '{extension}': expected '.vtu' or '.pvtu'.")

    reader: Any = readers[extension]()
    reader.SetFileName(file_path)
    reader.Update()
    mesh = check_grid(reader.GetOutput())
    if mesh.GetNumberOfPoints() == 0 or mesh.GetNumberOfCells() == 0:
        raise ValueError(
            f"Mesh file '{file_path}' has {mesh.GetNumberOfPoints()} points and {mesh.GetNumberOfCells()} cells."
        )
    update_grid(mesh, reset=reset_arrays)
    LOGGER.info(
        "Mesh read from '%s': %d points, %d cells.", file_path, mesh.GetNumberOfPoints(), mesh.GetNumberOfCells()
    )
    return mesh


def save_mesh(mesh: vtkUnstructuredGrid, path: PathType, binary: bool = True) -> None:
    """Write the mesh to a ``.vtu`` file, creating its directory if needed.

    Args:
        mesh: The mesh.
        path: File to write.
        binary: Binary (True) or ASCII (False) data.

    Raises:
        OSError: If VTK fails to write the file.
    """
    check_grid(mesh)
    file_path = os.fspath(path)
    directory = os.path.dirname(file_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(file_path)
    writer.SetInputData(mesh)
    if binary:
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAscii()
    if writer.Write() != 1:
        raise OSError(f"Failed to write the mesh to '{file_path}'.")
    LOGGER.debug("Mesh written to '%s'.", file_path)


# ---------------------------------------------------------------------------------------------------- summary


def summary(mesh: vtkUnstructuredGrid, name: str, path: str | None = None) -> str:
    """A text table describing the mesh: geometry and arrays, for the log.

    Args:
        mesh: The mesh.
        name: Name of the model.
        path: File the mesh was read from, if any.
    """
    width = 70
    rule = " " + "-" * width
    bounds = mesh.GetBounds()

    def array_rows(data: Any) -> list[str]:
        rows = []
        for i in range(data.GetNumberOfArrays()):
            array = data.GetAbstractArray(i)
            if array is None:
                continue
            n = array.GetNumberOfComponents()
            components = f"{n} comp" + ("s" if n != 1 else " ")
            rows.append(f"    {array.GetName() or '(unnamed)':<28}{components:<10}{array.GetDataTypeAsString()}")
        return rows

    def span(axis: int) -> str:
        return f"[{bounds[2 * axis]:.6g}, {bounds[2 * axis + 1]:.6g}]"

    lines = [
        "",
        rule,
        f"  {'pyTEC3D  Mesh Summary':^{width}}",
        rule,
        f"  {'Job name':<18}: {name}",
        f"  {'Working dir':<18}: {os.getcwd()}",
        f"  {'Mesh file':<18}: {path or '(in-memory grid)'}",
        rule,
        "  Geometry",
        f"    {'Points':<20}: {mesh.GetNumberOfPoints():,}",
        f"    {'Cells':<20}: {mesh.GetNumberOfCells():,}",
        f"    {'Spatial dimension':<20}: {'2D' if bounds[5] - bounds[4] < 1e-10 else '3D'}",
        f"    {'Bounds X':<20}: {span(0)}",
        f"    {'Bounds Y':<20}: {span(1)}",
        f"    {'Bounds Z':<20}: {span(2)}",
    ]
    for title, data in (
        ("Point arrays", mesh.GetPointData()),
        ("Cell arrays", mesh.GetCellData()),
        ("Field arrays", mesh.GetFieldData()),
    ):
        lines += [rule, f"  {title} ({data.GetNumberOfArrays()})", *array_rows(data)]
    lines.append(rule)
    return "\n".join(lines)
