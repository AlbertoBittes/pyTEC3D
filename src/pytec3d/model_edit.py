"""Geometric operations on the mesh of a model: coordinate mappings and derived fields.

Periodic couplings are made by :meth:`pytec3d.Model.periodic_boundaries`, which also records them in the model data.

Points are the mesh points (``location="node"``) or the cell centroids (``location="elem"``, the ``cell_centroids``
array of :func:`pytec3d.mesh.update_grid`). With ``selected_only=True`` an operation only applies to the current
selection (:mod:`pytec3d.selection`).
"""

from __future__ import annotations

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import VTK_DOUBLE, VTK_UNSIGNED_CHAR, vtkDataArray
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkPointDataToCellData

from . import mesh as meshes
from . import selection
from ._logging import get_logger
from ._typing import Location
from .coordinate_system import CoordinateSystem, CylindricalCoordinateSystem

LOGGER = get_logger("pytec3d.model_edit")


def _positions(mesh: vtkUnstructuredGrid, location: Location) -> np.ndarray:
    """Point coordinates or cell centroids, shape (N, 3)."""
    if meshes.check_location(location) == "node":
        return meshes.points(mesh)
    return meshes.array_values(mesh, meshes.CENTROIDS_ARRAY, "elem")


def _selected_rows(mesh: vtkUnstructuredGrid, location: Location, selected_only: bool) -> np.ndarray:
    """Ids of the selected points / cells, or of all of them."""
    if selected_only:
        return selection.selected_ids(mesh, location)
    return np.arange(meshes.size(mesh, location), dtype=int)


# ---------------------------------------------------------------------------------------------------- coordinates


def to_local_coordinates(
    mesh: vtkUnstructuredGrid,
    coordinate_system: CoordinateSystem,
    location: Location = "node",
    selected_only: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates of the points (or cell centroids) in a coordinate system.

    Args:
        mesh: The mesh.
        coordinate_system: The system (its own coordinates, e.g. r, theta, z).
        location: ``"node"`` for the points, ``"elem"`` for the cell centroids.
        selected_only: Only the selected points / cells.

    Returns:
        ``(coordinates, ids)``: shape (N, 3), and the point / cell id of each row.
    """
    positions = _positions(mesh, location)
    if positions.size == 0:
        return np.empty((0, 3)), np.empty(0, dtype=int)
    ids = _selected_rows(mesh, location, selected_only)
    if ids.size == 0:
        return positions[:0], ids
    return coordinate_system.to_local(positions[ids] if selected_only else positions), ids


def to_global_coordinates(
    mesh: vtkUnstructuredGrid,
    local_positions: np.ndarray,
    coordinate_system: CoordinateSystem,
    location: Location = "node",
    selected_only: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Global points of coordinates in a coordinate system, one row per point / cell (or per selected one).

    Args:
        mesh: The mesh (gives the ids).
        local_positions: Coordinates in ``coordinate_system``, one row per (selected) point or cell.
        coordinate_system: The system.
        location: ``"node"`` or ``"elem"``.
        selected_only: The rows are those of the selected points / cells.

    Returns:
        ``(global_points, ids)``.

    Raises:
        ValueError: If the number of rows does not match.
    """
    ids = _selected_rows(mesh, meshes.check_location(location), selected_only)
    if ids.size == 0:
        return np.empty((0, 3)), ids
    local_positions = np.asarray(local_positions, dtype=float)
    if local_positions.shape[0] != ids.size:
        raise ValueError(f"local_positions has {local_positions.shape[0]} rows, expected {ids.size}.")
    return coordinate_system.to_global(local_positions), ids


# ---------------------------------------------------------------------------------------------------- vector fields


def _mapped_vector_field(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    coordinate_system: CoordinateSystem,
    location: Location,
    ids: np.ndarray,
    to_local: bool,
) -> np.ndarray:
    """The rows ``ids`` of a 3-component array, mapped to local (or global) components.

    Rows that are entirely NaN stay NaN; NaN components of the other rows count as 0.
    """
    vtk_array = meshes.get_array(mesh, array_name, location)
    if vtk_array.GetNumberOfComponents() != 3:
        raise ValueError(f"Array '{array_name}' must have 3 components, it has {vtk_array.GetNumberOfComponents()}.")
    vectors = np.asarray(vtknp.vtk_to_numpy(vtk_array), dtype=float)[ids].copy()
    positions = _positions(mesh, location)[ids]
    valid = np.any(~np.isnan(vectors), axis=1)
    if np.any(valid):
        values = np.nan_to_num(vectors[valid], nan=0.0)
        if to_local:
            vectors[valid] = coordinate_system.vector_components_to_local(values, positions[valid])
        else:
            vectors[valid] = coordinate_system.vector_components_to_global(values, positions[valid])
    return vectors


def _map_vector_field(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    coordinate_system: CoordinateSystem,
    location: Location,
    selected_only: bool,
    to_local: bool,
) -> np.ndarray:
    vtk_array = meshes.get_array(mesh, array_name, location)
    values = vtknp.vtk_to_numpy(vtk_array)
    ids = _selected_rows(mesh, location, selected_only)
    if ids.size:
        values[ids] = _mapped_vector_field(mesh, array_name, coordinate_system, location, ids, to_local)
        vtk_array.Modified()
    return values


def map_vector_field_to_local_components(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    coordinate_system: CoordinateSystem,
    location: Location = "node",
    selected_only: bool = False,
) -> np.ndarray:
    """Replace the global Cartesian components of a vector array by its components in a coordinate system.

    E.g. x y z -> r theta z for a cylindrical system. Rows that are all NaN are left as they are.

    Args:
        mesh: The mesh, modified in place.
        array_name: A 3-component array.
        coordinate_system: The system.
        location: ``"node"`` or ``"elem"``.
        selected_only: Only map the selected points / cells.

    Returns:
        NumPy view of the updated array.

    Raises:
        MissingArrayError: If the array (or, for cells, ``cell_centroids``) is missing.
        ValueError: If the array does not have 3 components.
    """
    return _map_vector_field(mesh, array_name, coordinate_system, location, selected_only, to_local=True)


def map_vector_field_to_global_components(
    mesh: vtkUnstructuredGrid,
    array_name: str,
    coordinate_system: CoordinateSystem,
    location: Location = "node",
    selected_only: bool = False,
) -> np.ndarray:
    """Replace the components of a vector array in a coordinate system by its global Cartesian components.

    The inverse of :func:`map_vector_field_to_local_components`; same arguments, result and errors.
    """
    return _map_vector_field(mesh, array_name, coordinate_system, location, selected_only, to_local=False)


# ---------------------------------------------------------------------------------------------------- derived fields


def element_solution(mesh: vtkUnstructuredGrid, point_array_name: str) -> vtkDataArray:
    """Average a point array onto the cells (``vtkPointDataToCellData``) as the cell array ``<name>e``.

    Args:
        mesh: The mesh, modified in place.
        point_array_name: Point array to average.

    Returns:
        The cell array.

    Raises:
        MissingArrayError: If the point array is missing.
    """
    source = meshes.get_array(mesh, point_array_name, "node")
    averaging = vtkPointDataToCellData()
    averaging.SetInputData(mesh)
    averaging.SetProcessAllArrays(False)
    averaging.AddPointDataArray(point_array_name)
    averaging.PassPointDataOff()
    averaging.Update()
    averaged = averaging.GetOutput().GetCellData().GetArray(point_array_name)
    values = vtknp.vtk_to_numpy(averaged)

    target_name = f"{point_array_name}e"
    target = mesh.GetCellData().GetArray(target_name)
    if target is not None and target.GetNumberOfComponents() == source.GetNumberOfComponents():
        target_values = vtknp.vtk_to_numpy(target)
        if target_values.shape == values.shape:
            target_values[:] = values
            target.Modified()
            return target
    return meshes.add_array(
        mesh, target_name, "elem", components=averaged.GetNumberOfComponents(), dtype=averaged.GetDataType(),
        values=values, replace=True,
    )


def poloidal_flux(mesh: vtkUnstructuredGrid, coordinate_system: CoordinateSystem) -> dict[str, np.ndarray]:
    """Poloidal flux r * A_theta from the vector potential: point array ``PSI`` from ``A``, cell ``PSIe`` from ``Ae``.

    On a 2D (axisymmetric) mesh r is the x coordinate and A_theta the z component; on a 3D mesh both are taken in
    ``coordinate_system``.

    Args:
        mesh: The mesh, modified in place.
        coordinate_system: Cylindrical system of the torus.

    Returns:
        ``{"node": PSI values, "elem": PSIe values}``, for the source arrays present.

    Raises:
        ValueError: If the system is not cylindrical, or the mesh has neither ``A`` nor ``Ae``.
    """
    meshes.check_grid(mesh)
    if not isinstance(coordinate_system, CylindricalCoordinateSystem):
        raise ValueError("poloidal_flux needs a cylindrical coordinate system.")
    planar = mesh.GetMaxSpatialDimension() != 3

    def store(location: Location, source_name: str, target_name: str) -> np.ndarray:
        source = meshes.get_array(mesh, source_name, location)
        if source.GetNumberOfComponents() != 3:
            raise ValueError(f"Array '{source_name}' must have 3 components.")
        vectors = np.asarray(vtknp.vtk_to_numpy(source), dtype=float)
        positions = np.asarray(_positions(mesh, location), dtype=float)
        if planar:
            flux = positions[:, 0] * vectors[:, 2]
        else:
            flux = coordinate_system.to_local(positions)[:, 0] * coordinate_system.vector_components_to_local(
                vectors, positions
            )[:, 1]
        target = meshes.attributes(mesh, location).GetArray(target_name)
        if target is None or target.GetNumberOfComponents() != 1:
            target = meshes.add_array(mesh, target_name, location, components=1, dtype=VTK_DOUBLE, replace=True)
        target_values = vtknp.vtk_to_numpy(target)
        target_values[:] = flux
        target.Modified()
        return target_values

    results: dict[str, np.ndarray] = {}
    if mesh.GetPointData().GetArray("A") is not None:
        results["node"] = store("node", "A", "PSI")
    if mesh.GetCellData().GetArray("Ae") is not None:
        results["elem"] = store("elem", "Ae", "PSIe")
    if not results:
        raise ValueError("poloidal_flux needs the point array 'A' or the cell array 'Ae'.")
    return results


def add_component(mesh: vtkUnstructuredGrid, component_name: str, location: Location = "elem") -> vtkDataArray:
    """Store the current selection as a named component: a uint8 array, 1 on the selected points / cells.

    An existing array of that name is kept and overwritten with the selection.
    """
    mask = selection.selection_mask(mesh, location)
    array = meshes.add_array(mesh, component_name, location, components=1, fill=0, dtype=VTK_UNSIGNED_CHAR)
    vtknp.vtk_to_numpy(array)[:] = mask
    array.Modified()
    return array
