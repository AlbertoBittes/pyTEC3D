"""Selection of points and cells, in the style of the APDL selection commands.

The current selection is stored on the mesh, as the uint8 arrays ``selection_point`` (point data) and
``selection_cell`` (cell data), 1 where selected. Every selecting function takes a ``mode``:

* ``"s"``: replace the selection;
* ``"a"``: add to it;
* ``"r"``: reselect within it;
* ``"u"``: remove from it.

A typical sequence selects cells, then their nodes, then sets a value on them::

    from pytec3d import selection

    selection.select(mesh, "MATERIAL", 23, location="elem")
    selection.nsle(mesh)
    selection.set_value_on_selection(mesh, "CP_0", 1)
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import VTK_UNSIGNED_CHAR
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkSelection, vtkSelectionNode, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkFeatureEdges
from vtkmodules.vtkFiltersExtraction import vtkExtractSelection
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter

from . import mesh as meshes
from ._logging import get_logger
from ._typing import LOCATIONS, SELECTION_MODES, Location, SelectionMode

if TYPE_CHECKING:
    from .coordinate_system import CoordinateSystem

LOGGER = get_logger("pytec3d.selection")

#: Names of the selection arrays, by location.
SELECTION_ARRAYS: dict[Location, str] = {"node": "selection_point", "elem": "selection_cell"}
#: Pseudo array name that selects by point or cell index.
ID = "ID"

_ORIGINAL_IDS = "__pytec3d_original_point_ids__"


# ---------------------------------------------------------------------------------------------------- arrays


def ensure_selection_arrays(mesh: vtkUnstructuredGrid) -> None:
    """Add the selection arrays to the mesh, nothing selected, where they are missing."""
    meshes.check_grid(mesh)
    for location, name in SELECTION_ARRAYS.items():
        data = meshes.attributes(mesh, location)
        if data.GetArray(name) is None:
            array = vtknp.numpy_to_vtk(
                np.zeros(meshes.size(mesh, location), dtype=np.uint8), deep=True, array_type=VTK_UNSIGNED_CHAR
            )
            array.SetName(name)
            data.AddArray(array)


def selection_mask(mesh: vtkUnstructuredGrid, location: Location = "node") -> np.ndarray:
    """Boolean mask of the selected points or cells (a copy)."""
    ensure_selection_arrays(mesh)
    return meshes.array_values(mesh, SELECTION_ARRAYS[location], location) == 1


def selected_ids(mesh: vtkUnstructuredGrid, location: Location = "node") -> np.ndarray:
    """Indices of the selected points or cells, increasing."""
    return np.flatnonzero(selection_mask(mesh, location))


def apply_mask(mesh: vtkUnstructuredGrid, mask: np.ndarray, location: Location = "node",
               mode: SelectionMode = "s") -> None:
    """Combine a boolean mask with the current selection.

    Args:
        mesh: The mesh.
        mask: One boolean per point or cell.
        location: ``"node"`` or ``"elem"``.
        mode: ``"s"``, ``"a"``, ``"r"`` or ``"u"`` (see the module documentation).

    Raises:
        ValueError: If ``mode`` is unknown or ``mask`` has the wrong length.
    """
    if mode not in SELECTION_MODES:
        raise ValueError(f"Unknown selection mode {mode!r}: use one of {SELECTION_MODES}.")
    ensure_selection_arrays(mesh)
    array = meshes.get_array(mesh, SELECTION_ARRAYS[location], location)
    current = vtknp.vtk_to_numpy(array)
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != current.shape:
        raise ValueError(f"The selection mask has shape {mask.shape}, expected {current.shape}.")
    if mode == "a":
        mask = mask | (current == 1)
    elif mode == "r":
        mask = mask & (current == 1)
    elif mode == "u":
        mask = (current == 1) & ~mask
    current[:] = mask
    array.Modified()
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Selection (%s, mode %s): %d %s selected.", location, mode, int(np.count_nonzero(mask)),
                     "points" if location == "node" else "cells")


def _values(mesh: vtkUnstructuredGrid, array_name: str, component: int | str, location: Location) -> np.ndarray:
    """One value per point or cell: the index for ``"ID"``, a component, or the magnitude of the array."""
    if array_name == ID:
        return np.arange(meshes.size(mesh, location), dtype=np.int64)
    values = meshes.array_values(mesh, array_name, location)
    if component == "magnitude":
        return np.abs(values) if values.ndim == 1 else np.linalg.norm(values, axis=1)
    return values if values.ndim == 1 else values[:, int(component)]


def _matches(values: np.ndarray, value: Any) -> np.ndarray:
    """Mask of ``values`` matching ``value``: a range (tuple), a set of values (list, set, 1D array) or one value."""
    if isinstance(value, tuple):
        if len(value) != 2:
            raise ValueError("A range must be a tuple (low, high).")
        return (values >= value[0]) & (values <= value[1])
    if isinstance(value, (list, set, frozenset)):
        return np.isin(values, list(value))
    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            raise ValueError("An array of values must be 1D.")
        return np.isin(values, value)
    return values == value


# ---------------------------------------------------------------------------------------------------- select


def select(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    value: Any,
    component: int | str = 0,
    location: Location = "node",
    mode: SelectionMode = "s",
) -> None:
    """Select the points or cells whose array value matches ``value``.

    Args:
        mesh: The mesh.
        array_name: Array to compare, or ``"ID"`` for the point / cell index.
        value: One value (equality), a tuple ``(low, high)`` (inclusive range), or a list, set or 1D array of values.
        component: Component of a vector array, or ``"magnitude"``.
        location: ``"node"`` or ``"elem"``.
        mode: ``"s"``, ``"a"``, ``"r"`` or ``"u"``.

    Raises:
        MissingArrayError: If the mesh has no array ``array_name`` at ``location``.
        ValueError: If ``value`` or ``mode`` is invalid.
    """
    apply_mask(mesh, _matches(_values(mesh, array_name, component, location), value), location, mode)


def select_pos(
    mesh: vtkUnstructuredGrid,
    axis: int,
    value: Any,
    mode: SelectionMode = "s",
    coordinate_system: CoordinateSystem | None = None,
    location: Location = "node",
) -> None:
    """Select the points or cells (by centroid) whose coordinate along ``axis`` matches ``value``.

    Args:
        mesh: The mesh.
        axis: Coordinate index, 0, 1 or 2 (x, y, z, or the local coordinates of ``coordinate_system``).
        value: As in :func:`select`; usually a range ``(low, high)``.
        mode: ``"s"``, ``"a"``, ``"r"`` or ``"u"``.
        coordinate_system: Select on the local coordinates of this system (e.g. r, theta, z) instead of x, y, z.
        location: ``"node"`` or ``"elem"``.

    Raises:
        MissingArrayError: For cells, if the mesh has no ``cell_centroids`` array.
    """
    coordinates = (
        meshes.points(mesh) if location == "node"
        else np.asarray(meshes.array_values(mesh, meshes.CENTROIDS_ARRAY, "elem"), dtype=float)
    )
    if coordinate_system is not None:
        coordinates = coordinate_system.to_local(coordinates)
    apply_mask(mesh, _matches(coordinates[:, axis], value), location, mode)


def select_all(mesh: vtkUnstructuredGrid, location: Location | None = None) -> None:
    """Select every point, every cell, or both (``location=None``)."""
    _fill(mesh, location, lambda current: np.ones_like(current))


def unselect_all(mesh: vtkUnstructuredGrid, location: Location | None = None) -> None:
    """Unselect every point, every cell, or both (``location=None``)."""
    _fill(mesh, location, lambda current: np.zeros_like(current))


def inv(mesh: vtkUnstructuredGrid, location: Location | None = None) -> None:
    """Invert the selection of the points, of the cells, or of both (``location=None``)."""
    _fill(mesh, location, lambda current: 1 - current)


def _fill(mesh: vtkUnstructuredGrid, location: Location | None, update: Any) -> None:
    ensure_selection_arrays(mesh)
    for where in LOCATIONS if location is None else (meshes.check_location(location),):
        array = meshes.get_array(mesh, SELECTION_ARRAYS[where], where)
        current = vtknp.vtk_to_numpy(array)
        current[:] = update(current)
        array.Modified()


# ---------------------------------------------------------------------------------------------------- topology


def nsle(mesh: vtkUnstructuredGrid, mode: SelectionMode = "s") -> None:
    """Select the points of the selected cells (APDL ``NSLE``). Does nothing when no cell is selected."""
    cells = selection_mask(mesh, "elem")
    if not cells.any():
        return
    offsets, connectivity = meshes.cell_connectivity(mesh)
    if offsets.size < 2 or connectivity.size == 0:
        return
    point_mask = np.zeros(mesh.GetNumberOfPoints(), dtype=bool)
    point_mask[connectivity[np.repeat(cells, np.diff(offsets))]] = True
    if point_mask.any():
        apply_mask(mesh, point_mask, "node", mode)


def esln(mesh: vtkUnstructuredGrid, all_nodes: bool = True, mode: SelectionMode = "s") -> None:
    """Select the cells whose points are selected (APDL ``ESLN``). Does nothing when no point is selected.

    Args:
        mesh: The mesh.
        all_nodes: Select a cell when all its points are selected (True) or when any is (False).
        mode: ``"s"``, ``"a"``, ``"r"`` or ``"u"``.
    """
    selected_points = selection_mask(mesh, "node")
    if not selected_points.any():
        return
    offsets, connectivity = meshes.cell_connectivity(mesh)
    if offsets.size < 2:
        return
    lengths = np.diff(offsets)
    counts = np.add.reduceat(selected_points[connectivity].astype(np.int64), offsets[:-1])
    cell_mask = (lengths > 0) & (counts == lengths) if all_nodes else counts > 0
    if cell_mask.any():
        apply_mask(mesh, cell_mask, "elem", mode)


def _extract_cells(mesh: vtkUnstructuredGrid, cell_ids: np.ndarray) -> vtkUnstructuredGrid:
    node = vtkSelectionNode()
    node.SetFieldType(vtkSelectionNode.CELL)
    node.SetContentType(vtkSelectionNode.INDICES)
    node.SetSelectionList(vtknp.numpy_to_vtkIdTypeArray(np.asarray(cell_ids, dtype=np.int64), deep=True))
    cell_selection = vtkSelection()
    cell_selection.AddNode(node)
    extract = vtkExtractSelection()
    extract.SetInputData(0, mesh)
    extract.SetInputData(1, cell_selection)
    extract.Update()
    return extract.GetOutput()


def _external_boundary(mesh: vtkUnstructuredGrid, cell_ids: np.ndarray) -> vtkPolyData:
    """Boundary of the cells ``cell_ids``: the external surface (3D) or the boundary edges (2D).

    The points carry the point ids of ``mesh`` in the array ``_ORIGINAL_IDS``.
    """
    point_data = mesh.GetPointData()
    added_ids = point_data.GetArray(_ORIGINAL_IDS) is None
    if added_ids:
        ids = vtknp.numpy_to_vtkIdTypeArray(np.arange(mesh.GetNumberOfPoints(), dtype=np.int64), deep=True)
        ids.SetName(_ORIGINAL_IDS)
        point_data.AddArray(ids)
    try:
        surface = vtkGeometryFilter()
        surface.SetInputData(_extract_cells(mesh, cell_ids))
        surface.Update()
        if mesh.GetMaxSpatialDimension() > 2:
            boundary = vtkPolyData()
            boundary.DeepCopy(surface.GetOutput())
            return boundary
        edges = vtkFeatureEdges()
        edges.SetInputData(surface.GetOutput())
        edges.BoundaryEdgesOn()
        edges.FeatureEdgesOff()
        edges.ManifoldEdgesOff()
        edges.NonManifoldEdgesOff()
        edges.Update()
        boundary = vtkPolyData()
        boundary.DeepCopy(edges.GetOutput())
        return boundary
    finally:
        if added_ids:
            point_data.RemoveArray(_ORIGINAL_IDS)


def _original_point_ids(boundary: vtkPolyData) -> np.ndarray | None:
    data = boundary.GetPointData()
    array = data.GetArray(_ORIGINAL_IDS) or data.GetArray("vtkOriginalPointIds")
    return None if array is None else np.asarray(vtknp.vtk_to_numpy(array), dtype=np.int64)


def nsel_ext(mesh: vtkUnstructuredGrid, mode: SelectionMode = "s") -> None:
    """Select the points on the external boundary of the selected cells.

    The boundary is the external surface of 3D cells, or the boundary edges of 2D cells. Does nothing (and logs it)
    when no cell is selected.
    """
    cell_ids = selected_ids(mesh, "elem")
    if cell_ids.size == 0:
        LOGGER.info("nsel_ext: no cells are selected.")
        return
    point_ids = _original_point_ids(_external_boundary(mesh, cell_ids))
    if point_ids is None or point_ids.size == 0:
        LOGGER.info("nsel_ext: the selected cells have no external points.")
        return
    point_mask = np.zeros(mesh.GetNumberOfPoints(), dtype=bool)
    point_mask[point_ids] = True
    apply_mask(mesh, point_mask, "node", mode)


def external_node_adjacency(mesh: vtkUnstructuredGrid) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """The boundary of the selected cells of a 2D mesh, as a graph.

    Returns:
        ``(node_ids, adjacency)``: the boundary point ids (increasing), and for each the ids of its neighbours along
        the boundary edges (increasing).

    Raises:
        ValueError: If the mesh is 3D.
    """
    if mesh.GetMaxSpatialDimension() > 2:
        raise ValueError("external_node_adjacency() is for 2D meshes only.")
    cell_ids = selected_ids(mesh, "elem")
    if cell_ids.size == 0:
        return np.empty(0, dtype=np.int64), {}

    boundary = _external_boundary(mesh, cell_ids)
    original_ids = _original_point_ids(boundary)
    lines = boundary.GetLines()
    if original_ids is None or lines is None:
        return np.empty(0, dtype=np.int64), {}
    node_ids = np.unique(original_ids)
    if node_ids.size == 0:
        return node_ids, {}

    neighbours: dict[int, set[int]] = {int(node): set() for node in node_ids}
    offsets = np.asarray(vtknp.vtk_to_numpy(lines.GetOffsetsArray()), dtype=np.int64)
    connectivity = np.asarray(vtknp.vtk_to_numpy(lines.GetConnectivityArray()), dtype=np.int64)
    for start, stop in zip(offsets[:-1], offsets[1:], strict=True):
        path = original_ids[connectivity[start:stop]]
        for left, right in zip(path[:-1].tolist(), path[1:].tolist(), strict=True):
            if left != right:
                neighbours[left].add(right)
                neighbours[right].add(left)
    adjacency = {node: np.array(sorted(nodes), dtype=np.int64) for node, nodes in neighbours.items()}
    return node_ids, adjacency


# ---------------------------------------------------------------------------------------------------- values


def set_value_on_selection(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    value: Any,
    component: int = 0,
    location: Location = "node",
) -> None:
    """Set ``value`` in the array at the selected points or cells.

    Args:
        mesh: The mesh.
        array_name: Array to write.
        value: Value to set.
        component: Component of a vector array.
        location: ``"node"`` or ``"elem"``.

    Raises:
        MissingArrayError: If the mesh has no array ``array_name`` at ``location``.
    """
    array = meshes.get_array(mesh, array_name, location)
    values = vtknp.vtk_to_numpy(array)
    ids = selected_ids(mesh, location)
    if values.ndim == 1:
        values[ids] = value
    else:
        values[ids, component] = value
    array.Modified()
    LOGGER.debug("%s[%s] = %r on %d %s.", array_name, component, value, ids.size,
                 "points" if location == "node" else "cells")


def get_selection_values(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    component: int | str = 0,
    location: Location = "node",
) -> np.ndarray:
    """Values of the array at the selected points or cells, in increasing id order.

    Args:
        mesh: The mesh.
        array_name: Array to read, or ``"ID"`` for the ids themselves.
        component: Component of a vector array, or ``"magnitude"``.
        location: ``"node"`` or ``"elem"``.

    Raises:
        MissingArrayError: If the mesh has no array ``array_name`` at ``location``.
    """
    return _values(mesh, array_name, component, location)[selected_ids(mesh, location)]


__all__ = [
    "ID",
    "SELECTION_ARRAYS",
    "apply_mask",
    "ensure_selection_arrays",
    "esln",
    "external_node_adjacency",
    "get_selection_values",
    "inv",
    "nsel_ext",
    "nsle",
    "select",
    "select_all",
    "select_pos",
    "selected_ids",
    "selection_mask",
    "set_value_on_selection",
    "unselect_all",
]
