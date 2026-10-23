"""Plasma equilibrium and evolution (APEC-M): Picard iterations of the free-boundary Grad-Shafranov problem.

The plasma lives on the *plasma-region mesh*: a 2D (R, Z) mesh of the PLASMA_REGION cells. The solver mesh is
that same 2D mesh or a 3D torus built by revolving it about the y axis; toroidal maps then carry the plasma current
density from the plasma-region mesh to the solver mesh (projection) and the poloidal flux back (toroidal average).

Each Picard iteration:

1. finds the magnetic axis (PSI maximum), the X point (PSI saddle) and the limiter point, and from them the last closed
   flux surface; the cells inside it, connected to the axis, are the plasma (``PLASMA`` cell array);
2. computes the plasma current density ``Js`` from the Grad-Shafranov profile (``alpha_mp``, ``alpha_np``,
   ``beta_p``), scaled to the target ``plasma_current``, and, in a limiter configuration with a halo, the halo
   currents ``AMPS`` along the vacuum-chamber wall;
3. for an equilibrium, corrects the control-coil currents so that the field at the target axis is zero;
4. projects ``Js`` and ``AMPS`` to the solver mesh, solves, and averages the new flux ``PSI`` back;
5. stops when the largest change of PSI is below ``tolerance_value``; otherwise relaxes PSI and PSIe towards the
   previous iterate, ``(1 - w) PSI_old + w PSI_new``, with ``w`` halved when the iterations oscillate or stall.

Arrays read: ``PLASMA_REGION``, ``PLASMA``, ``Js``, ``CONSTANT``, the control-coil component arrays and
``cell_centroids`` / ``cell_volume`` (cells); ``PSI`` and the vector potential ``A`` (points). ``AMPS`` (points) is
optional: the mesh only has it when an element type has the field it loads (the voltage ``V`` of a solid conductor);
without it there are no halo currents.
"""

from __future__ import annotations

import contextlib
import copy
import csv
import math
import os
import time
from collections import deque
from typing import TYPE_CHECKING, Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import reference, vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCellLocator,
    vtkDataObject,
    vtkPolyData,
    vtkStaticPointLocator,
    vtkUnstructuredGrid,
)
from vtkmodules.vtkFiltersCore import vtkPointDataToCellData
from vtkmodules.vtkFiltersGeneral import vtkGradientFilter

from . import mesh as meshes
from . import model_edit, selection
from ._logging import get_logger
from ._typing import Location
from .coordinate_system import CylindricalCoordinateSystem, YAxisCylindricalCoordinateSystem
from .errors import ConvergenceError, SolverStateError

if TYPE_CHECKING:
    from .solver import Solver

LOGGER = get_logger("pytec3d.algorithm")

MU_0 = 4.0e-7 * np.pi

#: Point arrays of the PSI gradient and Hessian computed on the plasma-region mesh.
PSI_GRADIENT = "__pytec3d_psi_gradient__"
PSI_HESSIAN = "__pytec3d_psi_hessian__"
#: Load arrays the algorithm writes only when the solver mesh has them, i.e. when an element type of the model has
#: the field they load (AMPS: the nodal current of the solid-conductor voltage V).
OPTIONAL_LOADS = frozenset({"AMPS"})
#: Directory (in the working directory) and file of the per-iteration diagnostics.
#: The diagnostics CSV of the Picard iterations, in the output folder of the solver: <mesh_name><suffix>.
DIAGNOSTICS_SUFFIX = "_iteration_diagnostics.csv"

#: Default input parameters; override them in :attr:`Algorithm.input_parameters`.
DEFAULT_INPUT_PARAMETERS: dict[str, Any] = {
    # Target position of the magnetic axis (x = R, y = Z).
    "target_central_axis_x": 0.0,
    "target_central_axis_y": 0.0,
    "target_central_axis_z": 0.0,
    # Major radius and Grad-Shafranov profile.
    "R0": 0.0,
    "alpha_mp": 1.0001,
    "alpha_np": 0.83197,
    "beta_p": 0.0,
    "AF0": 5.3 * 6.2,
    "plasma_current": 0.0,
    # Picard iterations.
    "max_picard_iterations": 1,
    "tolerance_value": 1e-6,
    # Relaxation of PSI between Picard iterations: PSI <- (1 - w) PSI_old + w PSI_new, starting from w = relaxation
    # at each time step; w is multiplied by relaxation_factor (down to relaxation_min) when the solved flux comes
    # back to the one of two iterations before (2-cycle), or when max |dPSI| has not reached a new minimum in
    # stall_iterations iterations (stall), and divided by it (up to relaxation) when max |dPSI| has decreased at
    # each of the last stall_iterations iterations (0: no stall test, no recovery). relaxation_min = relaxation keeps
    # w fixed.
    "relaxation": 1.0,
    "relaxation_min": 0.1,
    "relaxation_factor": 0.5,
    "stall_iterations": 5,
    # Field at the target axis of a unit radial / vertical correction current.
    "bradx_value": 0.0,
    "brady_value": 0.0,
    "bverx_value": 0.0,
    "bvery_value": 0.0,
    # Meshes and toroidal extrusion of a 3D solver mesh.
    "equilibrium_mesh": "",
    "static_mesh": "",
    "evolution_mesh": "",
    "extrusion_angle": 360.0,
    "n_toroidal_divisions": -1,
    # Halo: psi_2_min = psi_min - delta1 * (1 + (t - t0) * delta2) * (psi_max - psi_min), from t0_halo on.
    "halo_delta1": 0.0,
    "halo_delta2": 0.0,
    "t0_halo": 0.0,
    # Control coils: each coil gets I = radial_gain * I_radial + vertical_gain * I_vertical, as a current density
    # I / cross_section added to the Js of its cells (component array of that name).
    "control_coils": [
        {"component": "PFC_2", "radial_gain": 1.0, "vertical_gain": 0.0},
        {"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 1.0},
        {"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 1.0},
        {"component": "PFC_5", "radial_gain": -1.0, "vertical_gain": 0.0},
    ],
}

#: Initial status; :attr:`Algorithm.status` holds the results of the last iteration.
DEFAULT_STATUS: dict[str, Any] = {
    "psi_magnetic_axis": 0.0,
    "psi_max": 0.0,
    "psi_min": 0.0,
    "psi_2_min": 0.0,
    "magnetic_axis_node_id": -1,
    "limiter_node_id": -1,
    "x_point_node_id": -1,
    "magnetic_axis": (0.0, 0.0, 0.0),
    "b_target_axis": (0.0, 0.0),
    "last_psi_limiter": 0.0,
    "last_psi_x_point": 0.0,
    "separatrix_configuration": "uninitialized",
    "limiter_point": (np.nan, np.nan, np.nan),
    "x_point": (np.nan, np.nan, np.nan),
    "flux_norm": np.nan,
    "flux_tolerance": np.nan,
    "relaxation": 1.0,
    "oscillation": "",
    "radial_correction_current": 0.0,
    "vertical_correction_current": 0.0,
    "total_plasma_current_before_scaling": 0.0,
    "plasma_current_scaling_lambda": 0.0,
    "plasma_element_count": 0,
    "picard_iteration": 0,
    "converged": False,
}


def _has_array(mesh: vtkUnstructuredGrid, name: str, location: Location) -> bool:
    data = mesh.GetPointData() if location == "node" else mesh.GetCellData()
    return data.GetArray(name) is not None


class Algorithm:
    """Plasma equilibrium and evolution solves driving a :class:`~pytec3d.Solver`.

    The algorithm only uses the public API of the solver (its model, ``time_step``, ``solve()``)::

        solver = Solver()
        solver.read_model("ITER_2D.vtu")
        solver.model.read_json("model_data.json")
        solver.load_model_into_solver()

        algorithm = Algorithm(solver)
        algorithm.input_parameters["R0"] = 6.2
        algorithm.equilibrium_solve()

    Attributes:
        solver: The solver the algorithm drives.
        input_parameters: Parameters of the algorithm (:data:`DEFAULT_INPUT_PARAMETERS`).
        status: Results of the last iteration (:data:`DEFAULT_STATUS`).
    """

    def __init__(self, solver: Solver, input_parameters: dict[str, Any] | None = None) -> None:
        """Create the algorithm for a solver.

        Args:
            solver: The solver whose model the algorithm drives.
            input_parameters: Overrides of :data:`DEFAULT_INPUT_PARAMETERS`.
        """
        self.solver = solver
        self.input_parameters: dict[str, Any] = copy.deepcopy(DEFAULT_INPUT_PARAMETERS)
        if input_parameters:
            self.input_parameters.update(input_parameters)
        self.status: dict[str, Any] = copy.deepcopy(DEFAULT_STATUS)
        self._equilibrium_mesh: vtkUnstructuredGrid | None = None
        self.reset()

    def reset(self) -> None:
        """Forget the iterations: status, previous flux, cached maps and graphs (the plasma-region mesh stays)."""
        self.status = copy.deepcopy(DEFAULT_STATUS)
        self._first_absolute_iteration = True
        self._previous_flux: np.ndarray | None = None
        self._previous_cell_flux: np.ndarray | None = None
        self._solved_fluxes: deque[np.ndarray] = deque(maxlen=2)
        self._flux_norms: list[float] = []
        self._relaxation_time: float | None = None
        self._control_current_density_base: np.ndarray | None = None
        self._b_target_cache: tuple[tuple[Any, ...], np.ndarray, np.ndarray] | None = None
        self._toroidal_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_cells = np.zeros((0, 0), dtype=np.int64)
        self._halo_node_groups: tuple[np.ndarray, ...] = ()
        self._halo_boundary: tuple[np.ndarray, dict[int, np.ndarray]] | None = None

    # ------------------------------------------------------------------------------------------- vacuum chamber

    @property
    def plasma_region_mesh(self) -> vtkUnstructuredGrid:
        """The 2D mesh of the vacuum chamber, where the plasma is computed.

        Set it from a 2D mesh with a ``PLASMA_REGION`` cell array: a copy of its PLASMA_REGION cells is kept.
        :meth:`equilibrium_solve` sets it from the solver mesh when it is not set.

        Raises:
            SolverStateError: When reading it before it is set, or setting it twice.
            ValueError: When the mesh has no PLASMA_REGION cell.
        """
        if self._equilibrium_mesh is None:
            raise SolverStateError("The plasma-region mesh is not set.")
        return self._equilibrium_mesh

    @plasma_region_mesh.setter
    def plasma_region_mesh(self, value: vtkUnstructuredGrid) -> None:
        if self._equilibrium_mesh is not None:
            raise SolverStateError("The plasma-region mesh is already set.")
        source = vtkUnstructuredGrid()
        source.DeepCopy(meshes.check_grid(value))
        selection.select(source, "PLASMA_REGION", 1, location="elem")
        cells = selection.selected_ids(source, "elem")
        if cells.size == 0:
            raise ValueError("The mesh has no PLASMA_REGION cell.")
        extracted = vtkUnstructuredGrid.SafeDownCast(selection._extract_cells(source, cells))
        if extracted is None or extracted.GetNumberOfCells() == 0:
            raise ValueError("Failed to extract the PLASMA_REGION cells.")
        selection.select_all(extracted, location="elem")
        selection.nsle(extracted)
        self._equilibrium_mesh = extracted
        self._toroidal_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_cells = np.zeros((0, 0), dtype=np.int64)
        self._halo_boundary = None

    # ------------------------------------------------------------------------------------------- toroidal maps

    @staticmethod
    def build_toroidal_slices(
        input_mesh: vtkUnstructuredGrid,
        target_mesh: vtkUnstructuredGrid,
        location: Location = "node",
        extrusion_angle: float = 360.0,
        n_divisions: int | None = None,
        coordinate_system: CylindricalCoordinateSystem | None = None,
    ) -> np.ndarray:
        """Map the selected points / cells of a 2D mesh to their images on a target mesh around the torus.

        Each selected entity of ``input_mesh`` is revolved about the axis of ``coordinate_system`` to the toroidal
        stations of the target, and matched with the nearest selected entity of ``target_mesh``.

        Args:
            input_mesh: The 2D (R, Z) mesh.
            target_mesh: The target mesh: the same 2D mesh, or the 3D torus.
            location: ``"node"`` (points) or ``"elem"`` (cell centroids).
            extrusion_angle: Toroidal extent of a 3D target, degrees.
            n_divisions: Number of cells along the torus of a 3D target. A 3D target has ``n_divisions + 1`` point
                stations or ``n_divisions`` cell stations; a 2D target one.
            coordinate_system: The cylindrical system of the torus; :class:`YAxisCylindricalCoordinateSystem` by
                default.

        Returns:
            Shape (selected inputs, 1 + stations): the input id, then the matched target id at each station.

        Raises:
            ValueError: On an invalid argument, a 3D input mesh, or an empty target selection.
        """
        location = meshes.check_location(location)
        if extrusion_angle <= 0.0:
            raise ValueError("extrusion_angle must be positive.")
        if input_mesh.GetNumberOfCells() > 0 and input_mesh.GetMaxSpatialDimension() > 2:
            raise ValueError("The input mesh must be 2D.")
        coordinate_system = coordinate_system or YAxisCylindricalCoordinateSystem()
        if not isinstance(coordinate_system, CylindricalCoordinateSystem):
            raise ValueError("The toroidal maps need a cylindrical coordinate system.")

        if target_mesh.GetMaxSpatialDimension() > 2:
            if n_divisions is None or int(n_divisions) <= 0:
                raise ValueError("A 3D target mesh needs a positive n_divisions.")
            n_divisions = int(n_divisions)
            if location == "node":
                stations = np.linspace(0.0, float(extrusion_angle), n_divisions + 1)
            else:
                stations = (np.arange(n_divisions, dtype=np.float64) + 0.5) * (float(extrusion_angle) / n_divisions)
        else:
            stations = np.zeros(1)
        n_stations = stations.size

        inputs, input_ids = model_edit.to_local_coordinates(input_mesh, coordinate_system, location, True)
        targets, target_ids = model_edit.to_local_coordinates(target_mesh, coordinate_system, location, True)
        if inputs.size == 0:
            return np.empty((0, 1 + n_stations), dtype=np.int64)
        if targets.size == 0:
            raise ValueError("The target mesh has no selected points / cells to map to.")
        inputs = np.asarray(inputs, dtype=np.float64)
        targets = np.asarray(targets, dtype=np.float64)

        # Cartesian positions in the frame of the y-axis system: x = R cos(theta), y = Z, z = R sin(theta).
        target_theta = np.deg2rad(targets[:, 1])
        target_xyz = np.column_stack(
            (targets[:, 0] * np.cos(target_theta), targets[:, 2], targets[:, 0] * np.sin(target_theta))
        )
        locator_points = vtkPoints()
        locator_points.SetData(vtknp.numpy_to_vtk(target_xyz, deep=True))
        locator_data = vtkPolyData()
        locator_data.SetPoints(locator_points)
        locator = vtkStaticPointLocator()
        locator.SetDataSet(locator_data)
        locator.AutomaticOn()
        locator.BuildLocator()

        theta = np.deg2rad(stations)
        n_inputs = inputs.shape[0]
        queries = np.empty((n_inputs * n_stations, 3))
        queries[:, 0] = (inputs[:, 0, None] * np.cos(theta)[None, :]).ravel()
        queries[:, 1] = np.broadcast_to(inputs[:, 2, None], (n_inputs, n_stations)).ravel()
        queries[:, 2] = (inputs[:, 0, None] * np.sin(theta)[None, :]).ravel()
        nearest = np.fromiter((locator.FindClosestPoint(q) for q in queries), dtype=np.int64, count=len(queries))

        mapping = np.empty((n_inputs, 1 + n_stations), dtype=np.int64)
        mapping[:, 0] = input_ids
        mapping[:, 1:] = np.asarray(target_ids)[nearest].reshape(n_inputs, n_stations)
        return mapping

    def _ensure_toroidal_maps(self) -> None:
        """Build the point and cell toroidal maps between the plasma-region mesh and the solver mesh, once."""
        if self._toroidal_nodes.size and self._toroidal_cells.size:
            return
        LOGGER.debug("Building the toroidal maps.")
        chamber = self.plasma_region_mesh
        target = self.solver.model.mesh
        divisions = int(self.input_parameters.get("n_toroidal_divisions", -1))
        if target.GetMaxSpatialDimension() > 2 and divisions <= 0:
            raise ValueError("input_parameters['n_toroidal_divisions'] must be positive for a 3D solver mesh.")
        if divisions <= 0:
            divisions = 1
        angle = float(self.input_parameters.get("extrusion_angle", 360.0))

        selection.select(target, "PLASMA_REGION", 1, location="elem")
        selection.nsle(target)
        selection.select(chamber, "PLASMA_REGION", 1, location="elem")
        self._toroidal_nodes = self.build_toroidal_slices(chamber, target, "node", angle, divisions)
        self._toroidal_cells = self.build_toroidal_slices(chamber, target, "elem", angle, divisions)

    def _map(self, location: Location) -> tuple[np.ndarray, np.ndarray]:
        self._ensure_toroidal_maps()
        mapping = self._toroidal_nodes if location == "node" else self._toroidal_cells
        return np.asarray(mapping[:, 0], dtype=np.int64), np.asarray(mapping[:, 1:], dtype=np.int64)

    def toroidal_average(self, array_name: str, location: Location, component: int | None = None) -> None:
        """Average an array of the solver mesh over the toroidal stations into the plasma-region mesh.

        Args:
            array_name: Array present on both meshes.
            location: ``"node"`` or ``"elem"``.
            component: Component of a vector array; all of them when None.

        Raises:
            MissingArrayError: If the array is missing on a mesh.
            ValueError: If ``component`` is out of range.
        """
        location = meshes.check_location(location)
        sources, targets = self._map(location)
        chamber_array = meshes.get_array(self.plasma_region_mesh, array_name, location)
        chamber = vtknp.vtk_to_numpy(chamber_array)
        model = meshes.array_values(self.solver.model.mesh, array_name, location)
        if chamber.ndim == 1:
            chamber[sources] = np.mean(model[targets], axis=1)
        else:
            for c in self._components(chamber, component, array_name):
                chamber[sources, c] = np.mean(model[targets, c], axis=1)
        chamber_array.Modified()

    def toroidal_projection(self, array_name: str, location: Location, component: int | None = None) -> None:
        """Copy an array of the plasma-region mesh to every toroidal station of the solver mesh.

        A nodal load (:data:`OPTIONAL_LOADS`, the halo currents ``AMPS``) is the load of the whole sector: it is not
        copied but shared among the node stations (trapezoid weights), so that the sector receives it once.

        Same arguments and errors as :meth:`toroidal_average`, except for an optional load: when the solver mesh does
        not have it (no element type has the field it loads) nothing is projected.
        """
        location = meshes.check_location(location)
        if array_name in OPTIONAL_LOADS and not _has_array(self.solver.model.mesh, array_name, location):
            LOGGER.debug("No '%s' array on the solver mesh (no field takes it): not projected.", array_name)
            return
        sources, targets = self._map(location)
        chamber = meshes.array_values(self.plasma_region_mesh, array_name, location)
        model_array = meshes.get_array(self.solver.model.mesh, array_name, location)
        model = vtknp.vtk_to_numpy(model_array)
        flat_targets = targets.reshape(-1)
        valid = flat_targets >= 0
        if not np.any(valid):
            return
        flat_targets = flat_targets[valid]
        if array_name in OPTIONAL_LOADS and chamber.ndim == 1 and targets.shape[1] > 1:
            # A nodal load of the chamber is the current of the whole sector: shared among the node stations with
            # the trapezoid weights (1/(2n) on the two end stations, 1/n inside), so that the sector receives it
            # once, also when the end stations are coupled by cyclic symmetry (the load of a slave node is carried
            # to its master).
            n = targets.shape[1] - 1
            weights = np.full(targets.shape[1], 1.0 / n)
            weights[[0, -1]] = 0.5 / n
            shares = (chamber[sources, np.newaxis] * weights[np.newaxis, :]).reshape(-1)[valid]
            model[flat_targets] = 0.0
            np.add.at(model, flat_targets, shares)
            model_array.Modified()
            return
        if chamber.ndim == 1:
            model[flat_targets] = np.broadcast_to(chamber[sources, np.newaxis], targets.shape).reshape(-1)[valid]
        else:
            for c in self._components(chamber, component, array_name):
                projected = np.broadcast_to(chamber[sources, c][:, np.newaxis], targets.shape).reshape(-1)[valid]
                model[flat_targets, c] = projected
        model_array.Modified()

    @staticmethod
    def _components(values: np.ndarray, component: int | None, array_name: str) -> list[int]:
        """All the components of a vector array, or the one asked for."""
        if component is None:
            return list(range(values.shape[1]))
        if not 0 <= component < values.shape[1]:
            raise ValueError(f"Component {component} out of range for array '{array_name}'.")
        return [component]

    def _swap_projected_js_components(self) -> None:
        """On a 3D solver mesh, turn the projected (R, Z, phi) components of Js into global Cartesian ones."""
        target = self.solver.model.mesh
        if target.GetMaxSpatialDimension() <= 2:
            return
        _, targets = self._map("elem")
        js_array = meshes.get_array(target, "Js", "elem")
        js = vtknp.vtk_to_numpy(js_array)
        if js.ndim != 2 or js.shape[1] < 3:
            raise ValueError("Js must be a vector array.")
        cells = targets.reshape(-1)
        cells = cells[cells >= 0]
        if cells.size == 0:
            return
        js[cells, :] = js[cells][:, [0, 2, 1]]
        selection.select(target, selection.ID, np.unique(cells), location="elem")
        model_edit.map_vector_field_to_global_components(
            target, "Js", YAxisCylindricalCoordinateSystem(), "elem", selected_only=True
        )
        js_array.Modified()

    # ------------------------------------------------------------------------------------------- flux

    def calculate_poloidal_flux(self) -> None:
        """Update ``Ae``, ``PSI`` and ``PSIe`` on the solver mesh from the vector potential ``A``.

        Raises:
            MissingArrayError: If the point array ``A`` is missing.
        """
        target = self.solver.model.mesh
        meshes.get_array(target, "A", "node")
        model_edit.element_solution(target, "A")
        model_edit.poloidal_flux(target, YAxisCylindricalCoordinateSystem())

    def _compute_psi_derivatives(self) -> None:
        """Store the PSI gradient and Hessian (``vtkGradientFilter``) as point arrays of the plasma-region mesh."""
        chamber = self.plasma_region_mesh
        gradient = vtkGradientFilter()
        gradient.SetInputData(chamber)
        gradient.SetInputScalars(vtkDataObject.FIELD_ASSOCIATION_POINTS, "PSI")
        gradient.SetResultArrayName(PSI_GRADIENT)
        gradient.SetComputeGradient(True)
        gradient.Update()
        gradient_array = gradient.GetOutput().GetPointData().GetArray(PSI_GRADIENT)
        if gradient_array is None:
            raise RuntimeError("Failed to compute the PSI gradient.")

        hessian = vtkGradientFilter()
        hessian.SetInputConnection(gradient.GetOutputPort())
        hessian.SetInputScalars(vtkDataObject.FIELD_ASSOCIATION_POINTS, PSI_GRADIENT)
        hessian.SetResultArrayName(PSI_HESSIAN)
        hessian.SetComputeGradient(True)
        hessian.Update()
        hessian_array = hessian.GetOutput().GetPointData().GetArray(PSI_HESSIAN)
        if hessian_array is None:
            raise RuntimeError("Failed to compute the PSI Hessian.")

        point_data = chamber.GetPointData()
        for array in (gradient_array, hessian_array):
            point_data.RemoveArray(array.GetName())
            point_data.AddArray(array)

    # ------------------------------------------------------------------------------------------- separatrix

    def _set_psi_magnetic_axis(self) -> None:
        """The magnetic axis: the point of largest PSI."""
        chamber = self.plasma_region_mesh
        psi = meshes.array_values(chamber, "PSI", "node")
        axis = int(np.argmax(psi))
        psi_axis = float(psi[axis])
        point = tuple(float(v) for v in meshes.points(chamber)[axis])
        self.status.update(
            psi_magnetic_axis=psi_axis, psi_max=psi_axis, magnetic_axis_node_id=axis, magnetic_axis=point
        )
        LOGGER.debug("Magnetic axis at (%.4g, %.4g, %.4g), PSI = %.6g.", *point, psi_axis)

    def _check_x_point(self) -> None:
        """The X point: among the PSI saddle points (det(Hessian) < 0), the one of smallest |grad PSI|.

        Raises:
            ValueError: If PSI has no saddle point.
        """
        chamber = self.plasma_region_mesh
        psi = meshes.array_values(chamber, "PSI", "node")
        hessian = meshes.array_values(chamber, PSI_HESSIAN, "node")
        gradient = meshes.array_values(chamber, PSI_GRADIENT, "node")
        if hessian.ndim != 2 or hessian.shape[1] != 9:
            raise ValueError(f"The PSI Hessian must have 9 components, got shape {hessian.shape}.")
        # Row-major tensor: [xx, xy, xz, yx, yy, yz, zx, zy, zz].
        determinant = hessian[:, 0] * hessian[:, 4] - hessian[:, 1] * hessian[:, 3]
        saddles = np.flatnonzero(np.isfinite(determinant) & (determinant < 0.0))
        if saddles.size == 0:
            raise ValueError("PSI has no saddle point: no X point.")
        if gradient.ndim != 2 or gradient.shape[1] < 2:
            raise ValueError(f"The PSI gradient must have at least 2 components, got shape {gradient.shape}.")
        magnitude = np.sqrt(np.square(gradient[:, 0]) + np.square(gradient[:, 1]))
        saddles = saddles[np.isfinite(magnitude[saddles])]
        if saddles.size == 0:
            raise ValueError("The PSI gradient is not finite at any saddle point.")
        x_point = int(saddles[int(np.argmin(magnitude[saddles]))])
        point = tuple(float(v) for v in meshes.points(chamber)[x_point])
        self.status.update(last_psi_x_point=float(psi[x_point]), x_point_node_id=x_point, x_point=point)
        LOGGER.debug("X point at (%.4g, %.4g, %.4g), PSI = %.6g.", *point, psi[x_point])

    def _check_limiter(self) -> None:
        """The limiter point: the wall point of largest PSI, among those where the poloidal field turns as in the core.

        Raises:
            ValueError: If the plasma-region mesh has no boundary point.
        """
        chamber = self.plasma_region_mesh
        psi = meshes.array_values(chamber, "PSI", "node")
        selection.select_all(chamber, location="elem")
        selection.nsel_ext(chamber)
        candidates = selection.selected_ids(chamber, "node")
        if candidates.size == 0:
            raise ValueError("The vacuum chamber has no boundary point to be a limiter.")

        # Drop the divertor / private-flux points, where (r - r_axis) . grad(PSI) has the sign opposite to the core's.
        gradient = meshes.array_values(chamber, PSI_GRADIENT, "node")
        coordinates = meshes.points(chamber)
        axis = int(self.status.get("magnetic_axis_node_id", -1))
        psi_x = self.status.get("last_psi_x_point")
        if 0 <= axis < chamber.GetNumberOfPoints():
            reference_psi = float(psi_x) if psi_x is not None else float(psi.max())
            core_sign = 1.0 if reference_psi >= psi[axis] else -1.0
            radial = coordinates[candidates, :2] - coordinates[axis, :2]
            alignment = radial[:, 0] * gradient[candidates, 0] + radial[:, 1] * gradient[candidates, 1]
            outward = core_sign * alignment > 0.0
            if np.any(outward):
                candidates = candidates[outward]

        limiter = int(candidates[int(np.argmax(psi[candidates]))])
        point = tuple(float(v) for v in coordinates[limiter])
        self.status.update(last_psi_limiter=float(psi[limiter]), limiter_node_id=limiter, limiter_point=point)
        LOGGER.debug("Limiter at (%.4g, %.4g, %.4g), PSI = %.6g.", *point, psi[limiter])

    def _halo_active(self) -> bool:
        """A halo is requested: limiter configuration, halo_delta1 > 0, from t0_halo on."""
        return (
            self.status["separatrix_configuration"] == "limiter"
            and float(self.input_parameters.get("halo_delta1", 0.0)) > 0.0
            and self.solver.time_step["time"] >= float(self.input_parameters.get("t0_halo", 0.0))
        )

    def find_separatrix(self) -> None:
        """Find the axis, X point and limiter, the last closed flux surface (``psi_min``), the halo boundary
        (``psi_2_min``), and mark the plasma cells.

        Raises:
            ValueError: If halo_delta1 is negative, or no X point / limiter is found.
        """
        self._reset_plasma_components()
        self._set_psi_magnetic_axis()
        self._check_x_point()
        self._check_limiter()
        psi_limiter = float(self.status["last_psi_limiter"])
        psi_x_point = float(self.status["last_psi_x_point"])
        if psi_limiter > psi_x_point:
            self.status.update(psi_min=psi_limiter, separatrix_configuration="limiter")
        else:
            self.status.update(psi_min=psi_x_point, separatrix_configuration="X-point")

        delta1 = float(self.input_parameters.get("halo_delta1", 0.0))
        if delta1 < 0.0:
            raise ValueError("input_parameters['halo_delta1'] must not be negative.")
        if self._halo_active():
            delta2 = float(self.input_parameters.get("halo_delta2", 0.0))
            t0 = float(self.input_parameters.get("t0_halo", 0.0))
            scale = delta1 * (1.0 + (self.solver.time_step["time"] - t0) * delta2)
            psi_min, psi_max = self.status["psi_min"], self.status["psi_max"]
            self.status["psi_2_min"] = psi_min - scale * (psi_max - psi_min)
        else:
            self.status["psi_2_min"] = self.status["psi_min"]
        self._identify_plasma_region()
        LOGGER.info(
            "Separatrix: %s, PSI_min = %.6g, PSI_2_min = %.6g, PSI_axis = %.6g.",
            self.status["separatrix_configuration"], self.status["psi_min"], self.status["psi_2_min"],
            self.status["psi_magnetic_axis"],
        )

    def _reset_plasma_components(self) -> None:
        """PLASMA = 1, AMPS = 0, Js = 0 on the whole plasma-region mesh."""
        chamber = self.plasma_region_mesh
        for name, location, value in (("PLASMA", "elem", 1), ("AMPS", "node", 0.0), ("Js", "elem", 0.0)):
            if name in OPTIONAL_LOADS and not _has_array(chamber, name, location):  # type: ignore[arg-type]
                continue
            array = meshes.get_array(chamber, name, location)  # type: ignore[arg-type]
            vtknp.vtk_to_numpy(array)[...] = value
            array.Modified()

    @staticmethod
    def _main_plasma_cells(
        chamber: vtkUnstructuredGrid, axis_node: int, psi: np.ndarray, psi_x: float
    ) -> np.ndarray:
        """The cells inside the flux surface ``psi_x`` connected to the magnetic axis.

        A candidate cell has all its points on the axis side of ``psi_x``. The connected region grows from the
        candidates touching the axis point through cells sharing at least two points, so it does not cross an X point
        (where two regions only share one point).

        Returns:
            The cell ids, increasing.
        """
        n_cells, n_points = chamber.GetNumberOfCells(), chamber.GetNumberOfPoints()
        if not 0 <= axis_node < n_points or not np.isfinite(psi[axis_node]):
            return np.zeros(0, dtype=np.int64)
        offsets, connectivity = meshes.cell_connectivity(chamber)
        lengths = np.diff(offsets)
        cell_of_entry = np.repeat(np.arange(n_cells), lengths)
        entry_psi = psi[connectivity]
        inside = entry_psi >= psi_x if psi[axis_node] >= psi_x else entry_psi <= psi_x
        outside_points = np.bincount(cell_of_entry, weights=~(np.isfinite(entry_psi) & inside), minlength=n_cells)
        candidate = outside_points == 0

        # Point -> candidate cells, as a CSR table.
        entries = candidate[cell_of_entry]
        entry_points, entry_cells = connectivity[entries], cell_of_entry[entries]
        order = np.lexsort((entry_cells, entry_points))
        cells_by_point = entry_cells[order].tolist()
        point_start = np.searchsorted(entry_points[order], np.arange(n_points + 1)).tolist()

        def candidates_at(point: int) -> list[int]:
            return cells_by_point[point_start[point]:point_start[point + 1]]

        seeds = candidates_at(axis_node)
        if not seeds:
            return np.zeros(0, dtype=np.int64)
        visited = np.zeros(n_cells, dtype=bool)
        visited[seeds] = True
        queue = deque(seeds)
        while queue:
            cell = queue.popleft()
            points = set(connectivity[offsets[cell]:offsets[cell + 1]].tolist())
            neighbours = {other for point in points for other in candidates_at(point)}
            for other in neighbours:
                if visited[other]:
                    continue
                shared = sum(1 for p in connectivity[offsets[other]:offsets[other + 1]].tolist() if p in points)
                if shared >= 2:
                    visited[other] = True
                    queue.append(other)
        return np.flatnonzero(visited).astype(np.int64)

    def _identify_plasma_region(self) -> None:
        """PLASMA = 1 on the cells that carry plasma current, 0 elsewhere; with a halo, find the halo wall paths.

        One search for the core and the halo: the cells with all their points inside the boundary of the current,
        connected to the magnetic axis (:meth:`_main_plasma_cells`). The boundary is ``psi_2_min`` (the limiter or
        X-point flux, or the halo boundary), but never beyond the X-point flux surface, which keeps the private flux
        region and the scrape-off layer beyond the separatrix out. The halo wall paths are the wall points of these
        cells: the interface between the plasma region and the wall.

        Raises:
            ValueError: If a halo is requested but psi_2_min >= psi_min.
        """
        chamber = self.plasma_region_mesh
        psi_2_min, psi_min = float(self.status["psi_2_min"]), float(self.status["psi_min"])
        halo_requested = self._halo_active()
        if halo_requested and psi_2_min >= psi_min:
            raise ValueError("psi_2_min must be below psi_min to identify the plasma region with a halo.")
        self._halo_node_groups = ()

        axis = int(self.status.get("magnetic_axis_node_id", -1))
        if axis < 0:
            return
        psi = meshes.array_values(chamber, "PSI", "node")
        psi_axis = float(psi[axis])
        boundary = psi_2_min
        psi_x = self.status.get("last_psi_x_point")
        if psi_x is not None:
            # The one of psi_2_min and the X-point flux nearer to the axis.
            inward = 1.0 if psi_axis >= float(psi_x) else -1.0
            boundary = inward * max(inward * psi_2_min, inward * float(psi_x))
        plasma = np.zeros(chamber.GetNumberOfCells(), dtype=bool)
        plasma[self._main_plasma_cells(chamber, axis, psi, boundary)] = True
        if not np.all(plasma):
            selection.apply_mask(chamber, ~plasma, "elem", "s")
            selection.set_value_on_selection(chamber, "PLASMA", 0, location="elem")

        if halo_requested:
            # Wall points of the plasma cells, as ordered paths along the wall.
            selection.select(chamber, "PLASMA", 1, location="elem")
            selection.nsle(chamber)
            wall_points, adjacency = self._halo_boundary_graph()
            halo_points = np.intersect1d(selection.selected_ids(chamber, "node"), wall_points, assume_unique=True)
            self._halo_node_groups = self._orient_halo_node_groups(
                self._connected_node_groups(halo_points, adjacency),
                meshes.points(chamber),
                self.status["magnetic_axis"],
                float(self.input_parameters["plasma_current"]),
            )
        LOGGER.debug(
            "Plasma region: %d cells (boundary PSI = %.6g), %d halo group(s).",
            int(np.count_nonzero(plasma)), boundary, len(self._halo_node_groups),
        )

    def _halo_boundary_graph(self) -> tuple[np.ndarray, dict[int, np.ndarray]]:
        """The wall of the vacuum chamber as a graph of points (computed once)."""
        if self._halo_boundary is None:
            chamber = self.plasma_region_mesh
            selection.select_all(chamber)
            self._halo_boundary = selection.external_node_adjacency(chamber)
        return self._halo_boundary

    @staticmethod
    def _connected_node_groups(nodes: np.ndarray, adjacency: dict[int, np.ndarray]) -> tuple[np.ndarray, ...]:
        """Split points into the paths they form along the wall graph, each ordered from one end to the other.

        Returns:
            The paths, sorted by their first point.

        Raises:
            ValueError: If a group branches or closes on itself.
        """
        nodes = np.unique(np.asarray(nodes, dtype=np.int64).reshape(-1))
        if nodes.size <= 1:
            return (nodes,) if nodes.size else ()

        def neighbours_of(node: int) -> list[int]:
            return adjacency[node].tolist() if node in adjacency else []

        selected = set(nodes.tolist())
        visited: set[int] = set()
        groups: list[np.ndarray] = []
        for start in nodes.tolist():
            if start in visited:
                continue
            stack, component = [start], []
            visited.add(start)
            while stack:
                node = stack.pop()
                component.append(node)
                for neighbour in neighbours_of(node):
                    if neighbour in selected and neighbour not in visited:
                        visited.add(neighbour)
                        stack.append(neighbour)
            if len(component) <= 2:
                groups.append(np.asarray(component, dtype=np.int64))
                continue

            members = set(component)
            neighbours = {node: [n for n in neighbours_of(node) if n in members] for node in component}
            if any(len(n) > 2 for n in neighbours.values()):
                raise ValueError("A halo point group branches: it cannot be ordered as one path along the wall.")
            ends = [node for node, n in neighbours.items() if len(n) == 1]
            if len(ends) != 2:
                raise ValueError("A halo point group is not an open path along the wall.")
            ordered = [min(ends)]
            previous = -1
            while True:
                following = [n for n in neighbours[ordered[-1]] if n != previous and n not in ordered]
                if not following:
                    break
                previous = ordered[-1]
                ordered.append(following[0])
            if len(ordered) != len(component):
                raise ValueError("A halo point group is not a simple path along the wall.")
            groups.append(np.asarray(ordered, dtype=np.int64))
        return tuple(sorted(groups, key=lambda group: int(group[0])))

    @staticmethod
    def _orient_halo_node_groups(
        groups: tuple[np.ndarray, ...],
        coordinates: np.ndarray,
        axis_point: tuple[float, float, float],
        plasma_current: float,
    ) -> tuple[np.ndarray, ...]:
        """Orient each path: clockwise around the axis for a positive plasma current, counterclockwise otherwise."""
        if not groups or plasma_current == 0.0:
            return groups
        axis_xy = np.asarray(axis_point[:2], dtype=float)
        clockwise_wanted = plasma_current > 0.0
        oriented = []
        for group in groups:
            if group.size > 1:
                xy = np.asarray(coordinates[group, :2], dtype=float) - axis_xy
                cross = xy[:-1, 0] * xy[1:, 1] - xy[:-1, 1] * xy[1:, 0]
                dot = np.einsum("ij,ij->i", xy[:-1], xy[1:])
                turn = float(np.sum(np.arctan2(cross, dot)))
                if abs(turn) > np.finfo(float).eps and (turn < 0.0) != clockwise_wanted:
                    group = group[::-1].copy()
            oriented.append(group)
        return tuple(oriented)

    # ------------------------------------------------------------------------------------------- plasma current

    def compute_jplasma_gradshafranov(self) -> None:
        """Compute the plasma current density Js of the PLASMA cells, and the halo currents AMPS on the wall.

        Toroidal: ``J_phi = lambda * (1 - psi_n^alpha_mp)^alpha_np * (beta_p R / R0 + (1 - beta_p) R0 / R)``, with
        ``psi_n`` the normalized flux of the cell and ``lambda`` scaling the total to ``plasma_current``. Poloidal:
        ``J_p = F'(psi) B_p`` from the same profile.
        """
        chamber = self.plasma_region_mesh
        js_array = meshes.get_array(chamber, "Js", "elem")
        js = vtknp.vtk_to_numpy(js_array)
        psi_point = meshes.array_values(chamber, "PSI", "node")
        psi_cell = meshes.array_values(chamber, "PSIe", "elem")
        meshes.get_array(chamber, PSI_GRADIENT, "node")

        averaging = vtkPointDataToCellData()
        averaging.SetInputData(chamber)
        averaging.SetProcessAllArrays(False)
        averaging.AddPointDataArray(PSI_GRADIENT)
        averaging.Update()
        cell_gradient = vtknp.vtk_to_numpy(averaging.GetOutput().GetCellData().GetArray(PSI_GRADIENT))
        radius = meshes.array_values(chamber, meshes.CENTROIDS_ARRAY, "elem")[:, 0]
        # Poloidal field: B_R = -dPSI/dZ / R, B_Z = dPSI/dR / R.
        b_r = -cell_gradient[:, 1] / radius
        b_z = cell_gradient[:, 0] / radius
        volume = meshes.array_values(chamber, meshes.SIZES_ARRAY, "elem")

        psi_axis = float(self.status["psi_magnetic_axis"])
        psi_min = float(self.status["psi_2_min"])
        psi_span = psi_min - psi_axis
        p = self.input_parameters
        if abs(psi_span) <= np.finfo(float).eps:
            self.status.update(total_plasma_current_before_scaling=0.0, plasma_current_scaling_lambda=0.0,
                               plasma_element_count=0)
            self._require_plasma(0)
            return

        beta = float(p["beta_p"])
        alpha_m, alpha_n = float(p["alpha_mp"]), float(p["alpha_np"])
        r0, af0 = float(p["R0"]), float(p["AF0"])
        target_current = float(p["plasma_current"])

        selection.select(chamber, "PLASMA", 1, location="elem")
        plasma = selection.selected_ids(chamber, "elem")
        self._require_plasma(plasma.size)

        psi_n = np.nan_to_num((psi_cell[plasma] - psi_axis) / psi_span, nan=0.0, posinf=1.0, neginf=0.0)
        psi_n = np.clip(psi_n, 0.0, 1.0)
        one_minus = np.maximum(1.0 - np.power(psi_n, alpha_m), 0.0)
        profile = np.power(one_minus, alpha_n)
        density = profile * (beta * radius[plasma] / r0 + (1.0 - beta) * (r0 / radius[plasma]))
        total = float(np.sum(density * volume[plasma]))
        scale = target_current / total if abs(total) > np.finfo(float).eps else 0.0

        c_const = 2.0 * MU_0 * r0 * (1.0 - beta) * (psi_axis - psi_min) * scale / (alpha_n + 1.0)
        f_value = np.sqrt(np.maximum(c_const * np.power(one_minus, alpha_n + 1.0) + af0 * af0, 0.0))
        f_prime = (profile * (1.0 - beta) * r0 * scale) / f_value
        js[plasma, 0] = f_prime * b_r[plasma]
        js[plasma, 1] = f_prime * b_z[plasma]
        js[plasma, 2] = density * scale
        js_array.Modified()

        amps_array = chamber.GetPointData().GetArray("AMPS")
        halo_currents = (
            self.solver.model.mesh.GetMaxSpatialDimension() > 2
            and self._halo_active()
            and float(self.status["psi_2_min"]) < float(self.status["psi_min"])
        )
        if halo_currents and amps_array is not None and self._halo_node_groups:
            amps = vtknp.vtk_to_numpy(amps_array)
            angle_factor = float(p.get("extrusion_angle", 360.0)) / 360.0
            for path in self._halo_node_groups:
                if path.size <= 1:
                    continue
                # Current between two wall points: the change of the poloidal current function F along the wall.
                psi_seg = np.clip((psi_point[path] - psi_axis) / psi_span, 0.0, 1.0)
                integrand = np.power(np.maximum(1.0 - np.power(psi_seg, alpha_m), 0.0), alpha_n)
                integral = np.concatenate(([0.0], np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(psi_seg))))
                f_const = 2.0 * MU_0 * r0 * (1.0 - beta) * psi_span * scale
                f_profile = np.sqrt(np.maximum(f_const * integral + af0 * af0, 0.0))
                half_currents = 0.5 * (2.0 * np.pi / MU_0) * (f_profile[:-1] - f_profile[1:]) * angle_factor
                np.add.at(amps, path[:-1], half_currents)
                np.add.at(amps, path[1:], half_currents)
            amps_array.Modified()

        self.status.update(
            total_plasma_current_before_scaling=total,
            plasma_current_scaling_lambda=float(scale),
            plasma_element_count=int(plasma.size),
        )
        LOGGER.info("Grad-Shafranov: %d plasma cells, I_raw = %.6g A, lambda = %.6g.", plasma.size, total, scale)

    def _require_plasma(self, cell_count: int) -> None:
        """Stop the Picard loop when a plasma current is requested but no cell carries it.

        Raises:
            ConvergenceError: If plasma_current is not zero and there is no plasma cell.
        """
        if cell_count == 0 and float(self.input_parameters["plasma_current"]) != 0.0:
            raise ConvergenceError(
                f"The plasma region is empty at Picard iteration {self.status['picard_iteration']} "
                f"(t = {self.solver.time_step['time']}, {self.status['separatrix_configuration']} configuration, "
                f"PSI_axis = {self.status['psi_magnetic_axis']}, PSI_2_min = {self.status['psi_2_min']})."
            )

    # ------------------------------------------------------------------------------------------- control

    def control_currents(self) -> None:
        """Correct the control-coil currents so that the poloidal field at the target axis vanishes."""
        self._set_b_target_axis()
        self._compute_correction_currents()
        self._apply_control_currents()
        LOGGER.info(
            "Control currents: radial %.6g A, vertical %.6g A.",
            self.status["radial_correction_current"], self.status["vertical_correction_current"],
        )

    def _set_b_target_axis(self) -> None:
        """The poloidal field (B_R, B_Z) at the target axis, from the PSI gradient interpolated in its cell.

        Raises:
            ValueError: If the target point is outside the vacuum chamber or at R <= 0.
        """
        chamber = self.plasma_region_mesh
        gradient = meshes.array_values(chamber, PSI_GRADIENT, "node")
        p = self.input_parameters
        target = tuple(float(p[f"target_central_axis_{axis}"]) for axis in "xyz")
        if target[0] <= 0.0:
            raise ValueError(f"The target axis must be at R > 0, got R = {target[0]}.")

        key = (id(chamber), chamber.GetNumberOfPoints(), chamber.GetNumberOfCells(), target)
        if self._b_target_cache is None or self._b_target_cache[0] != key:
            locator = vtkCellLocator()
            locator.SetDataSet(chamber)
            locator.BuildLocator()
            cell_id = int(locator.FindCell(list(target)))
            if cell_id < 0:
                raise ValueError(f"The target axis {target} is outside the vacuum chamber.")
            cell = chamber.GetCell(cell_id)
            n = int(cell.GetNumberOfPoints())
            weights = [0.0] * n
            inside = cell.EvaluatePosition(
                list(target), [0.0, 0.0, 0.0], reference(0), [0.0, 0.0, 0.0], reference(0.0), weights  # type: ignore[arg-type]
            )
            if int(inside) < 0:
                raise ValueError(f"Cannot interpolate at the target axis {target} in cell {cell_id}.")
            ids = cell.GetPointIds()
            point_ids = np.array([ids.GetId(i) for i in range(n)], dtype=np.int64)
            self._b_target_cache = (key, point_ids, np.asarray(weights, dtype=float))
        _, point_ids, weights_np = self._b_target_cache

        d_psi = weights_np @ gradient[point_ids, :2]
        b_target = (float(-d_psi[1] / target[0]), float(d_psi[0] / target[0]))
        self.status["b_target_axis"] = b_target
        LOGGER.debug("Field at the target axis: B_R = %.6g, B_Z = %.6g.", *b_target)

    def _compute_correction_currents(self) -> None:
        """Radial and vertical correction currents cancelling the field at the target axis (2x2 linear system)."""
        p = self.input_parameters
        bradx, brady = float(p["bradx_value"]), float(p["brady_value"])
        bverx, bvery = float(p["bverx_value"]), float(p["bvery_value"])
        btx, bty = self.status["b_target_axis"]
        determinant = bradx * bvery - bverx * brady
        if abs(determinant) <= np.finfo(float).eps:
            radial = vertical = 0.0
        else:
            radial = (bverx * bty - bvery * btx) / determinant
            vertical = (brady * btx - bradx * bty) / determinant
        self.status.update(radial_correction_current=float(radial), vertical_correction_current=float(vertical))
        LOGGER.debug("Correction currents: radial %.6g A, vertical %.6g A.", radial, vertical)

    def _control_coils(self) -> list[dict[str, Any]]:
        """``input_parameters["control_coils"]``, checked: component name and the two gains of each coil."""
        coils = self.input_parameters.get("control_coils") or []
        if not isinstance(coils, list):
            raise TypeError("input_parameters['control_coils'] must be a list of dicts.")
        out = []
        for index, coil in enumerate(coils, start=1):
            if not isinstance(coil, dict):
                raise TypeError(f"control_coils[{index}] must be a dict, got {type(coil).__name__}.")
            component = str(coil.get("component", "")).strip()
            missing = [k for k in ("radial_gain", "vertical_gain") if k not in coil]
            if not component or missing:
                raise ValueError(f"control_coils[{index}] needs a 'component', 'radial_gain' and 'vertical_gain'.")
            out.append({"component": component, "radial_gain": float(coil["radial_gain"]),
                        "vertical_gain": float(coil["vertical_gain"])})
        return out

    def _apply_control_currents(self) -> None:
        """Add the correction current density of each control coil to the Js of its cells on the solver mesh.

        Raises:
            ValueError: If a coil cell has a constants row without (or with a zero) ``cross_section``.
        """
        coils = self._control_coils()
        if not coils:
            LOGGER.debug("No control coils.")
            return
        target = self.solver.model.mesh
        js_array = meshes.get_array(target, "Js", "elem")
        js = vtknp.vtk_to_numpy(js_array)
        radial = float(self.status["radial_correction_current"])
        vertical = float(self.status["vertical_correction_current"])

        in_coil = np.zeros(js.shape[0], dtype=bool)
        current = np.zeros(js.shape[0])
        for coil in coils:
            cells = meshes.array_values(target, coil["component"], "elem") == 1
            in_coil |= cells
            current[cells] += coil["radial_gain"] * radial + coil["vertical_gain"] * vertical
        cells = np.flatnonzero(in_coil)
        if cells.size == 0:
            LOGGER.debug("The control coils have no cells.")
            return

        constants_ids = meshes.array_values(target, "CONSTANT", "elem")[cells].astype(int)
        areas = {int(i): self.solver.model.real_constant(int(i), "cross_section") for i in np.unique(constants_ids)}
        missing = sorted(i for i, area in areas.items() if area is None)
        if missing:
            raise ValueError(f"The constants rows {missing} of control-coil cells have no 'cross_section'.")
        cell_areas = np.array([float(areas[i]) for i in constants_ids.tolist()])
        zero = np.abs(cell_areas) <= np.finfo(float).eps
        if np.any(zero):
            raise ValueError(f"The constants rows {np.unique(constants_ids[zero]).tolist()} have cross_section = 0.")

        if self._control_current_density_base is None or self._control_current_density_base.shape[0] != js.shape[0]:
            self._control_current_density_base = np.asarray(js[:, 2], dtype=float).copy()
        js[cells, 2] = self._control_current_density_base[cells] + current[cells] / cell_areas
        js_array.Modified()

    # ------------------------------------------------------------------------------------------- iterations

    def check_convergence(self) -> bool:
        """Compare PSI with the PSI the iteration started from: converged when ``max |dPSI| <= tolerance_value``.

        ``max |dPSI|`` is the fixed-point residual, whatever the relaxation. When not converged, PSI and PSIe of the
        plasma-region mesh are relaxed towards the previous iterate (:meth:`_relax_flux`). The first call only stores
        PSI and returns False.
        """
        chamber = self.plasma_region_mesh
        psi = meshes.array_values(chamber, "PSI", "node").copy()
        psi_cell = meshes.array_values(chamber, "PSIe", "elem").copy() if _has_array(chamber, "PSIe", "elem") else None
        tolerance = float(self.input_parameters["tolerance_value"])
        self.status["flux_tolerance"] = tolerance
        if self._first_absolute_iteration or self._previous_flux is None:
            self._previous_flux, self._previous_cell_flux = psi, psi_cell
            self.status.update(flux_norm=np.nan, converged=False)
            return False
        if self._previous_flux.size != psi.size:
            raise ValueError("PSI changed size between two iterations.")
        flux_norm = float(np.max(np.abs(psi - self._previous_flux)))
        converged = flux_norm <= tolerance
        self.status.update(flux_norm=flux_norm, converged=converged)
        LOGGER.info("Convergence: max |dPSI| = %.6g, tolerance %.6g: %s.", flux_norm, tolerance,
                    "converged" if converged else "not converged")
        if converged:
            self.status["oscillation"] = ""
            self._previous_flux, self._previous_cell_flux = psi, psi_cell
        else:
            self._detect_oscillation(psi, flux_norm)
            self._relax_flux(psi, psi_cell)
        return converged

    def _reset_relaxation(self) -> None:
        """Start a time step from w = relaxation, with no oscillation history (kept across calls at the same time)."""
        p = self.input_parameters
        relaxation, minimum = float(p["relaxation"]), float(p["relaxation_min"])
        factor, stall = float(p["relaxation_factor"]), int(p["stall_iterations"])
        if not 0.0 < relaxation <= 1.0 or not 0.0 < minimum <= 1.0:
            raise ValueError("input_parameters['relaxation'] and ['relaxation_min'] must be in (0, 1].")
        if not 0.0 < factor < 1.0:
            raise ValueError("input_parameters['relaxation_factor'] must be in (0, 1).")
        if stall < 0:
            raise ValueError("input_parameters['stall_iterations'] must not be negative.")
        time_now = float(self.solver.time_step["time"])
        if self._relaxation_time == time_now:
            return
        self._relaxation_time = time_now
        self._solved_fluxes.clear()
        self._flux_norms = []
        self.status.update(relaxation=relaxation, oscillation="")

    def _detect_oscillation(self, psi: np.ndarray, flux_norm: float) -> None:
        """Reduce the relaxation factor when the iterations oscillate or stall, raise it back when they converge.

        2-cycle: the solved PSI comes back within half its last change to the one of two iterations before. Stall:
        ``max |dPSI|`` has not gone below its earlier minimum in the last ``stall_iterations`` iterations (a slow but
        monotone convergence is not a stall). Recovery: ``max |dPSI|`` has decreased at each of the last
        ``stall_iterations`` iterations, so a reduced factor does not slow down the rest of the convergence. After a
        change of the factor the history restarts, so the next one needs new evidence.
        """
        p = self.input_parameters
        oscillation = ""
        if len(self._solved_fluxes) == 2:
            change = float(np.max(np.abs(psi - self._solved_fluxes[-1])))
            if change > float(self.status["flux_tolerance"]) and \
                    float(np.max(np.abs(psi - self._solved_fluxes[0]))) < 0.5 * change:
                oscillation = "2-cycle"
        stall = int(p["stall_iterations"])
        self._flux_norms.append(flux_norm)
        if not oscillation and stall > 0 and len(self._flux_norms) > stall and \
                min(self._flux_norms[-stall:]) >= min(self._flux_norms[:-stall]):
            oscillation = "stall"
        self._solved_fluxes.append(psi)

        relaxation = float(self.status["relaxation"])
        minimum, initial, factor = float(p["relaxation_min"]), float(p["relaxation"]), float(p["relaxation_factor"])
        recent = self._flux_norms[-1 - stall:]
        converging = stall > 0 and len(recent) == stall + 1 and all(b < a for a, b in zip(recent, recent[1:]))
        if oscillation and relaxation > minimum:
            relaxation = max(minimum, relaxation * factor)
            self._solved_fluxes.clear()
            self._flux_norms = []
            LOGGER.warning("Picard %s detected (max |dPSI| = %.6g): relaxation reduced to %.4g.",
                           oscillation, flux_norm, relaxation)
        elif not oscillation and converging and relaxation < initial:
            relaxation = min(initial, relaxation / factor)
            self._solved_fluxes.clear()
            self._flux_norms = []
            LOGGER.info("Picard converging (max |dPSI| = %.6g): relaxation raised to %.4g.", flux_norm, relaxation)
        elif oscillation:
            LOGGER.info("Picard %s detected (max |dPSI| = %.6g) at the minimum relaxation %.4g.",
                        oscillation, flux_norm, relaxation)
        self.status.update(relaxation=relaxation, oscillation=oscillation)

    def _relax_flux(self, psi: np.ndarray, psi_cell: np.ndarray | None) -> None:
        """PSI <- (1 - w) PSI_old + w PSI on the plasma-region mesh (PSIe too), and store it as the next start.

        PSI and PSIe are linear in the field, so the relaxed PSIe is the cell flux of the relaxed field; the axis,
        X point, limiter and plasma region of the next iteration all come from the relaxed flux.
        """
        chamber = self.plasma_region_mesh
        w = float(self.status["relaxation"])
        if w < 1.0 and self._previous_flux is not None:
            psi = w * psi + (1.0 - w) * self._previous_flux
            values = [("PSI", "node", psi)]
            if psi_cell is not None and self._previous_cell_flux is not None:
                psi_cell = w * psi_cell + (1.0 - w) * self._previous_cell_flux
                values.append(("PSIe", "elem", psi_cell))
            for name, location, relaxed in values:
                array = meshes.get_array(chamber, name, location)  # type: ignore[arg-type]
                vtknp.vtk_to_numpy(array)[...] = relaxed
                array.Modified()
        self._previous_flux, self._previous_cell_flux = psi, psi_cell

    def _diagnostics_path(self) -> str:
        from pytec3d.solver import OUTPUT_DIRECTORY
        directory = os.path.join(os.getcwd(), OUTPUT_DIRECTORY)
        os.makedirs(directory, exist_ok=True)
        return os.path.join(directory, self.solver.model.mesh_name + DIAGNOSTICS_SUFFIX)

    def _clear_diagnostics(self) -> None:
        if self.solver.rank == 0:
            with contextlib.suppress(FileNotFoundError):
                os.remove(self._diagnostics_path())

    def _diagnostics_row(self) -> dict[str, Any]:
        """One CSV row of the status; points become one column per coordinate."""
        row: dict[str, Any] = {
            "time": float(self.solver.time_step["time"]),
            "picard_iteration": int(self.status["picard_iteration"]),
            "flux_metric": float(self.status["flux_norm"]),
            "flux_threshold": float(self.status["flux_tolerance"]),
            "converged": bool(self.status["converged"]),
        }
        expanded = {"magnetic_axis": "magnetic_axis", "limiter_point": "limiter", "x_point": "x_point",
                    "b_target_axis": "b_target_axis"}
        for key, value in self.status.items():
            if key in ("picard_iteration", "converged", "flux_norm", "flux_tolerance"):
                continue
            if key in expanded:
                for axis, coordinate in zip("xyz", value, strict=False):
                    row[f"{expanded[key]}_{axis}"] = float(coordinate)
            else:
                row[key] = value
        return row

    def _append_diagnostics(self) -> None:
        """Append the status to the diagnostics CSV (rank 0); NaN is written as an empty field."""
        if self.solver.rank != 0:
            return
        row = {k: "" if isinstance(v, float) and math.isnan(v) else v for k, v in self._diagnostics_row().items()}
        path = self._diagnostics_path()
        new_file = not os.path.exists(path) or os.path.getsize(path) == 0
        with open(path, "a", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            if new_file:
                writer.writeheader()
            writer.writerow(row)

    def _run_picard_loop(self, label: str, equilibrium: bool) -> bool:
        """The Picard iterations shared by the equilibrium and evolution solves.

        Raises:
            ConvergenceError: If PSI does not converge within ``max_picard_iterations``.
        """
        max_iterations = max(1, int(self.input_parameters["max_picard_iterations"]))
        self.status["converged"] = False
        if self._first_absolute_iteration:
            self._ensure_toroidal_maps()
            self._halo_boundary_graph()
        self._reset_relaxation()

        for iteration in range(1, max_iterations + 1):
            self.status["picard_iteration"] = iteration
            LOGGER.info("%s Picard iteration %d/%d (t = %g).", label, iteration, max_iterations,
                        self.solver.time_step["time"])
            started = time.perf_counter()
            if self._first_absolute_iteration:
                self.check_convergence()
            else:
                self._compute_psi_derivatives()
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()
                if equilibrium:
                    self.control_currents()

            self.toroidal_projection("Js", "elem")
            self._swap_projected_js_components()
            self.toroidal_projection("AMPS", "node")
            self.solver.solve()
            self.calculate_poloidal_flux()
            self.toroidal_average("PSI", "node")
            self.toroidal_average("PSIe", "elem")

            converged = self.check_convergence() or (not equilibrium and self._first_absolute_iteration)
            self.status["converged"] = converged
            self._append_diagnostics()
            self._first_absolute_iteration = False
            LOGGER.info("%s Picard iteration %d done in %.2f s.", label, iteration, time.perf_counter() - started)
            if converged:
                LOGGER.info("%s solve converged at iteration %d (max |dPSI| = %s).", label, iteration,
                            self.status["flux_norm"])
                return True

        raise ConvergenceError(
            f"The {label.lower()} solve did not converge in {max_iterations} Picard iterations "
            f"(last max |dPSI| = {self.status['flux_norm']})."
        )

    def _require_session(self, what: str) -> None:
        if not self.solver.is_loaded:
            raise SolverStateError(f"No solver session: call solver.load_model_into_solver() before {what}().")

    def equilibrium_solve(self) -> bool:
        """Iterate the free-boundary equilibrium, with control currents, until PSI converges.

        The plasma-region mesh is taken from the solver mesh if it is not set.

        Returns:
            True (it raises otherwise).

        Raises:
            SolverStateError: If the model is not loaded into a solver session.
            ConvergenceError: If PSI does not converge within ``max_picard_iterations``.
        """
        LOGGER.info("Starting the equilibrium solve.")
        self._require_session("equilibrium_solve")
        if self._first_absolute_iteration:
            self._clear_diagnostics()
        if self._equilibrium_mesh is None:
            self.plasma_region_mesh = self.solver.model.mesh
        return self._run_picard_loop("Equilibrium", equilibrium=True)

    def evolution_solve(self) -> bool:
        """Iterate one time step of the plasma evolution (no control currents) until PSI converges.

        The first call only initializes the flux and returns after one solve.

        Returns:
            True (it raises otherwise).

        Raises:
            SolverStateError: If the model is not loaded or the plasma-region mesh is not set.
            ConvergenceError: If PSI does not converge within ``max_picard_iterations``.
        """
        LOGGER.info("Starting the evolution solve.")
        self._require_session("evolution_solve")
        if self._first_absolute_iteration:
            self._clear_diagnostics()
        return self._run_picard_loop("Evolution", equilibrium=False)
