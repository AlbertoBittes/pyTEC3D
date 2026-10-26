"""Checks of a model, of the solver settings and of the time step, before they reach the C++ solver.

Each check collects every problem before raising one :class:`~pytec3d.errors.ValidationError` that lists them all.
Problems that do not prevent a solve (unused ids, no loads...) are logged as warnings.

The model data is checked against the registries compiled into the solver module (domains, element types, material
properties, real constants, laws): names, value shapes, ids against the mesh, laws and their curves. Which properties a
kernel needs depends on the element type and domain; the solver checks those when the model is loaded.
"""

from __future__ import annotations

import importlib
from collections.abc import Mapping
from functools import cache
from typing import TYPE_CHECKING, Any, TypeGuard

import numpy as np

from . import element_library
from . import mesh as meshes
from ._logging import get_logger
from ._typing import Location
from .coordinate_system import coordinate_system_from_dict
from .errors import ModelLockedError, ValidationError

if TYPE_CHECKING:
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

    from .model import Model

LOGGER = get_logger("pytec3d.validation")

#: KSP types accepted in the solver settings.
KSP_TYPES = frozenset({
    "preonly", "gmres", "fgmres", "bcgs", "bcgsl", "cgs", "cg", "minres", "symmlq", "lgmres", "ibcgs", "gcr",
    "tfqmr", "tcqmr", "cr", "lsqr", "bicg", "richardson", "chebyshev",
})
#: PC types accepted in the solver settings (lu and cholesky use MUMPS).
PC_TYPES = frozenset({
    "lu", "cholesky", "ilu", "icc", "asm", "gasm", "bjacobi", "jacobi", "sor", "eisenstat", "hypre", "ml", "gamg",
    "mg", "ksp", "none",
})

#: Shapes a parameter of each type may have: a number stands for a vector with equal components.
_SHAPES: dict[str, list[tuple[int, ...]]] = {"scalar": [()], "vector": [(), (3,)], "tensor": [(3, 3)]}
_SHAPE_NAMES = {"scalar": "a number", "vector": "a number or 3 numbers", "tensor": "a 3x3 array"}

#: Keys of a row of the element_types section.
_ELEMENT_TYPE_KEYS = {"element_type", "option1", "option2", "option3"}
_CP_ARRAY = meshes.COUPLING_ARRAY


def _is_real(value: Any) -> TypeGuard[float]:
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)


def _is_index(value: Any, size: int) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(value, bool) and 0 <= int(value) < size


# ---------------------------------------------------------------------------------------------------- arrays


def check_required_arrays(mesh: vtkUnstructuredGrid, element_types: list[str] | None = None) -> None:
    """Check that the mesh has the model arrays and those of the element types, with the right shapes.

    Args:
        mesh: The mesh.
        element_types: Element type names; all of the library when None.

    Raises:
        ValidationError: Listing every missing or misshapen array.
    """
    meshes.check_grid(mesh)
    problems: list[str] = []
    cell_data = mesh.GetCellData()
    for name in meshes.MODEL_CELL_ARRAYS:
        if cell_data.GetArray(name) is None:
            problems.append(f"cell array '{name}' is missing.")
    angles = mesh.GetPointData().GetArray(meshes.NODE_ANGLES_ARRAY)
    if angles is None:
        problems.append(f"point array '{meshes.NODE_ANGLES_ARRAY}' is missing.")
    elif angles.GetNumberOfComponents() != 3:
        problems.append(
            f"point array '{meshes.NODE_ANGLES_ARRAY}' has {angles.GetNumberOfComponents()} component(s), expected 3."
        )

    for grid_array in element_library.grid_arrays(element_types):
        where = "point" if grid_array.location == "node" else "cell"
        array = meshes.attributes(mesh, grid_array.location).GetArray(grid_array.name)
        if array is None:
            problems.append(f"{where} array '{grid_array.name}' is missing.")
            continue
        if array.GetNumberOfComponents() != grid_array.components:
            problems.append(
                f"{where} array '{grid_array.name}' has {array.GetNumberOfComponents()} component(s),"
                f" expected {grid_array.components}."
            )
        expected = meshes.size(mesh, grid_array.location)
        if array.GetNumberOfTuples() != expected:
            problems.append(
                f"{where} array '{grid_array.name}' has {array.GetNumberOfTuples()} tuples, expected {expected}."
            )
    if problems:
        raise ValidationError("The mesh misses arrays of the model", problems)


# ---------------------------------------------------------------------------------------------------- solver settings


def check_solver_settings(settings: Mapping[str, Any]) -> None:
    """Check the solver settings (:data:`pytec3d.solver.DEFAULT_SOLVER_SETTINGS` merged with the user's).

    * ``linear``: ``ksp_type`` and ``pc_type`` names of :data:`KSP_TYPES` and :data:`PC_TYPES`; tolerances finite
      and > 0, the absolute one not above the relative one; ``divergence_tolerance`` None or > 0;
      ``max_iterations`` >= 1; ``initial_guess_nonzero`` a bool; ``factor_solver_type`` None or a name;
    * ``nonlinear``: ``max_iterations`` >= 1; tolerances finite and >= 0; ``minimum_step_length`` in (0, 1];
      ``linear_tolerance_multiplier`` > 0; the switches bools;
    * ``petsc_options``: names starting with "-", values strings or numbers.

    Unknown keys are errors, so that a misspelled setting is not silently ignored.

    Raises:
        ValidationError: Listing every invalid setting.
    """
    from .solver import DEFAULT_SOLVER_SETTINGS

    problems: list[str] = []
    unknown_groups = sorted(set(settings) - set(DEFAULT_SOLVER_SETTINGS))
    if unknown_groups:
        problems.append(f"unknown groups {unknown_groups}; use {sorted(DEFAULT_SOLVER_SETTINGS)}.")

    def group(name: str) -> Mapping[str, Any]:
        values = settings.get(name, {})
        if not isinstance(values, Mapping):
            problems.append(f"{name} must be a dict, got {type(values).__name__}.")
            return {}
        unknown = sorted(set(values) - set(DEFAULT_SOLVER_SETTINGS[name]))
        if unknown:
            problems.append(f"{name}: unknown settings {unknown}; known: {sorted(DEFAULT_SOLVER_SETTINGS[name])}.")
        return values

    def positive(values: Mapping[str, Any], name: str, key: str, allow_zero: bool = False) -> None:
        value = values.get(key)
        if not _is_real(value) or not np.isfinite(value) or value < 0.0 or (value == 0.0 and not allow_zero):
            problems.append(f"{name}.{key} must be a finite number {'>=' if allow_zero else '>'} 0, got {value!r}.")

    def integer(values: Mapping[str, Any], name: str, key: str) -> None:
        value = values.get(key)
        if not isinstance(value, (int, np.integer)) or isinstance(value, bool) or value < 1:
            problems.append(f"{name}.{key} must be an integer >= 1, got {value!r}.")

    def flags(values: Mapping[str, Any], name: str, keys: tuple[str, ...]) -> None:
        for key in keys:
            if not isinstance(values.get(key), (bool, np.bool_)):
                problems.append(f"{name}.{key} must be a bool, got {values.get(key)!r}.")

    linear = group("linear")
    for key, allowed in (("ksp_type", KSP_TYPES), ("pc_type", PC_TYPES)):
        value = linear.get(key)
        if not isinstance(value, str) or value not in allowed:
            problems.append(f"linear.{key} {value!r} is not supported; use one of {sorted(allowed)}.")
    positive(linear, "linear", "relative_tolerance")
    positive(linear, "linear", "absolute_tolerance")
    rtol, atol = linear.get("relative_tolerance"), linear.get("absolute_tolerance")
    if _is_real(rtol) and _is_real(atol) and atol > rtol:
        problems.append("linear.absolute_tolerance must not exceed linear.relative_tolerance.")
    if linear.get("divergence_tolerance") is not None:
        positive(linear, "linear", "divergence_tolerance")
    integer(linear, "linear", "max_iterations")
    flags(linear, "linear", ("initial_guess_nonzero",))
    factor = linear.get("factor_solver_type")
    if factor is not None and (not isinstance(factor, str) or not factor):
        problems.append(f"linear.factor_solver_type must be a name (e.g. 'mumps') or None, got {factor!r}.")

    nonlinear = group("nonlinear")
    integer(nonlinear, "nonlinear", "max_iterations")
    for key in ("absolute_residual_tolerance", "relative_residual_tolerance", "absolute_step_tolerance",
                "relative_step_tolerance"):
        positive(nonlinear, "nonlinear", key, allow_zero=True)
    step = nonlinear.get("minimum_step_length")
    if not _is_real(step) or not 0.0 < step <= 1.0:
        problems.append(f"nonlinear.minimum_step_length must be in (0, 1], got {step!r}.")
    positive(nonlinear, "nonlinear", "linear_tolerance_multiplier")
    flags(nonlinear, "nonlinear", ("require_residual_reduction", "require_finite_residual", "brent_line_search",
                                   "continue_after_max_iterations", "continue_after_backtrack_failure", "verbose"))

    options = settings.get("petsc_options", {})
    if not isinstance(options, Mapping):
        problems.append(f"petsc_options must be a dict, got {type(options).__name__}.")
    else:
        for name, value in options.items():
            if not isinstance(name, str) or not name.startswith("-"):
                problems.append(f"petsc_options: {name!r} must be an option name starting with '-'.")
            if value is not None and not isinstance(value, (str, int, float)):
                problems.append(f"petsc_options[{name!r}] must be a string or a number, got {value!r}.")

    if problems:
        raise ValidationError("Invalid solver settings", problems)


# ---------------------------------------------------------------------------------------------------- time step


def check_time_step(time_step: Mapping[str, Any]) -> None:
    """Check the time step.

    * ``time``: end time of the step, finite and > 0;
    * ``theta``: theta method (first-order formulations), in (0, 1];
    * ``beta`` (default 0.25) and ``gamma`` (default 0.5): Newmark parameters (second-order formulations), in
      (0, 0.5] and [0, 1];
    * ``transient_effects``: a bool (or 0 / 1): False for a static solve.

    Raises:
        ValidationError: Listing every invalid entry.
    """
    if not isinstance(time_step, Mapping):
        raise ValidationError("Invalid time step", [f"time_step must be a dict, got {type(time_step).__name__}."])
    problems: list[str] = []

    def check(key: str, low: float, high: float, low_open: bool, default: Any = None) -> None:
        value = time_step.get(key, default)
        valid = _is_real(value) and np.isfinite(value) and (value > low if low_open else value >= low) and value <= high
        if not valid:
            interval = f"{'(' if low_open else '['}{low:g}, {high:g}]"
            problems.append(f"time_step['{key}'] must be in {interval}, got {value!r}.")

    check("time", 0.0, np.inf, low_open=True)
    check("theta", 0.0, 1.0, low_open=True)
    check("beta", 0.0, 0.5, low_open=True, default=0.25)
    check("gamma", 0.0, 1.0, low_open=False, default=0.5)
    transient = time_step.get("transient_effects")
    if not (isinstance(transient, (bool, np.bool_)) or (_is_real(transient) and transient in (0, 1))):
        problems.append(f"time_step['transient_effects'] must be a bool (or 0 / 1), got {transient!r}.")
    if problems:
        raise ValidationError("Invalid time step", problems)


#: The arrays that define the model in a solver session, read when the model is loaded: the cell ids (the element
#: matrices are built from them once), the element coordinate systems, the nodal frames and the DOF couplings
#: (``CP_k``, added by :func:`load_time_arrays`).
LOAD_TIME_ARRAYS: tuple[tuple[Location, str], ...] = (
    ("elem", "TYPE"), ("elem", "MATERIAL"), ("elem", "CONSTANT"), ("elem", "COORDINATE_SYSTEM"),
    ("node", meshes.NODE_ANGLES_ARRAY),
)


def load_time_arrays(mesh: vtkUnstructuredGrid) -> dict[tuple[Location, str], np.ndarray]:
    """Copies of the arrays of :data:`LOAD_TIME_ARRAYS` and of the ``CP_k`` arrays, taken when a model is loaded."""
    point_data = mesh.GetPointData()
    names = list(LOAD_TIME_ARRAYS) + [
        ("node", name) for name in (point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays()))
        if name and _CP_ARRAY.match(name)
    ]
    return {
        (location, name): np.array(meshes.array_values(mesh, name, location), copy=True)
        for location, name in names
        if meshes.attributes(mesh, location).GetArray(name) is not None
    }


def check_load_time_arrays(mesh: vtkUnstructuredGrid, snapshot: Mapping[tuple[Location, str], np.ndarray]) -> None:
    """Compare the arrays of a loaded model with their copies taken at the load (:func:`load_time_arrays`).

    A change would be ignored by the session, or mixed with the old values (``MATERIAL`` and ``CONSTANT`` are read
    again by the loads and the output fields). The cost is one comparison of the arrays, about 1 % of a solve on
    large models.

    Raises:
        ModelLockedError: If one of the arrays changed or was removed.
    """
    for (location, name), values in snapshot.items():
        present = meshes.attributes(mesh, location).GetArray(name) is not None
        if not present or not np.array_equal(meshes.array_values(mesh, name, location), values):
            raise ModelLockedError(
                f"The array '{name}' changed after the model was loaded. The session builds the element matrices"
                " (TYPE, MATERIAL, CONSTANT), the element coordinate systems, the nodal frames and the DOF couplings"
                " when the model is loaded: set them before load_model_into_solver(), or close the session and load"
                " the model again."
            )


def check_before_solve(
    model: Model,
    time_step: Mapping[str, Any],
    loaded_arrays: Mapping[tuple[Location, str], np.ndarray] | None = None,
) -> None:
    """Checks run before every solve: the time step, the arrays, and warnings for loads and constraints.

    Warns when every load array is zero (nothing drives the problem) or every Dirichlet array is NaN (no boundary
    condition).

    Args:
        model: The loaded model.
        time_step: The step to solve.
        loaded_arrays: The copies of :func:`load_time_arrays` taken at the load, to check they did not change.

    Raises:
        ValidationError: If the time step is invalid or arrays are missing.
        ModelLockedError: If an array of :data:`LOAD_TIME_ARRAYS` or a ``CP_k`` array changed since the load.
    """
    check_time_step(time_step)
    types = model.element_type_names() or None
    mesh = model.mesh
    check_required_arrays(mesh, types)
    if loaded_arrays is not None:
        check_load_time_arrays(mesh, loaded_arrays)

    loads = element_library.load_arrays(types)
    if loads and all(np.all(np.nan_to_num(meshes.array_values(mesh, a.name, a.location)) == 0.0) for a in loads):
        LOGGER.warning("No loads: the arrays %s are all zero.", ", ".join(a.name for a in loads))
    constraints = ["_" + name for name in element_library.dof_fields(types)]
    if constraints and all(np.all(np.isnan(meshes.array_values(mesh, name, "node"))) for name in constraints):
        LOGGER.warning("No boundary conditions: the arrays %s are all NaN.", ", ".join(constraints))


# ---------------------------------------------------------------------------------------------------- model data


@cache
def solver_registries() -> dict[str, Any]:
    """Domains, fields, element types, material properties, real constants and laws compiled into the solver.

    Raises:
        ImportError: If the solver module is not built.
    """
    try:
        module = importlib.import_module("pytec3d_solver.pytec3d_solver")
    except ImportError as exc:
        raise ImportError(
            "The model data is checked against the registries of the C++ solver module "
            "'pytec3d_solver.pytec3d_solver', which could not be imported: build it first (build.sh)."
        ) from exc
    registries: dict[str, Any] = module.registries()
    return registries


class _ModelChecker:
    """The checks of :func:`check_model`, one method per section of the model data."""

    def __init__(self, model: Model) -> None:
        registries = solver_registries()
        self.model = model
        self.data = model.data
        self.mesh = model.mesh
        self.problems: list[str] = []
        self.domains = list(registries["domains"])
        self.element_type_names = list(registries["element_types"])
        self.element_type_dimensions = [int(d) for d in registries["element_types"].values()]
        self.material_properties: dict[str, str] = dict(registries["material_properties"])
        self.real_constants: dict[str, str] = dict(registries["real_constants"])
        self.laws: dict[str, dict[str, Any]] = dict(registries.get("laws", {}))
        self.law_names = list(self.laws)

        def cell_ids(name: str) -> np.ndarray:
            return np.asarray(meshes.array_values(self.mesh, name, "elem"), dtype=np.int64).reshape(-1)

        self.type_ids = cell_ids("TYPE")
        self.material_ids = cell_ids("MATERIAL")
        self.constants_ids = cell_ids("CONSTANT")
        self.coordinate_system_ids = cell_ids("COORDINATE_SYSTEM")

    def run(self) -> None:
        self.coordinate_systems()
        materials = self.property_rows("materials", self.material_properties, "material property")
        self.mesh_ids("materials", materials, self.material_ids, missing_is_error=True)
        element_types = self.element_types()
        self.dof_couplings()
        constants = self.property_rows("constants", self.real_constants, "real constant")
        self.mesh_ids("constants", constants, self.constants_ids, missing_is_error=False)
        self.nonlinear_laws()
        if self.problems:
            raise ValidationError("The model is invalid", self.problems)
        LOGGER.info(
            "Model data checked: %d element types, %d materials, %d real constants.",
            len(element_types), len(materials), len(constants),
        )

    # -------------------------------------------------------------------------------------------- helpers

    def section(self, name: str) -> dict[int, Any]:
        """``data[name]`` keyed by int id; ids that are not positive integers are reported."""
        entries = self.data.get(name, {})
        if not isinstance(entries, Mapping):
            self.problems.append(f"data['{name}'] must be an object, got {type(entries).__name__}.")
            return {}
        out: dict[int, Any] = {}
        for key, value in entries.items():
            try:
                entry_id = int(key)
                valid = entry_id > 0 and str(entry_id) == str(key).strip()
            except (TypeError, ValueError):
                valid = False
            if not valid:
                self.problems.append(f"{name}: the id {key!r} must be a positive integer.")
            elif entry_id in out:
                self.problems.append(f"{name}: the id {entry_id} is given more than once.")
            else:
                out[entry_id] = value
        return out

    def value_shape(self, value: Any, value_type: str, label: str) -> None:
        """A number, a 3-vector (or a number for equal components) or a 3x3 tensor, finite."""
        try:
            array = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            self.problems.append(f"{label}: must be numeric.")
            return
        if array.shape not in _SHAPES[value_type]:
            self.problems.append(f"{label}: must be {_SHAPE_NAMES[value_type]}, got shape {array.shape}.")
        elif not np.all(np.isfinite(array)):
            self.problems.append(f"{label}: must be finite.")

    def mesh_ids(self, name: str, rows: dict[int, Any], ids: np.ndarray, missing_is_error: bool) -> None:
        """Ids of the cells against the rows of a section.

        Ids used by cells but not defined are errors (or warnings); defined but unused ones are warnings.
        """
        if ids.size == 0:
            return
        used = set(np.unique(ids).tolist())
        missing = sorted(i for i in used - set(rows) if missing_is_error or i > 0)
        if missing:
            message = f"{name}: cell ids absent from data['{name}']: {missing}."
            if missing_is_error:
                self.problems.append(message)
            else:
                LOGGER.warning(message)
        unused = sorted(set(rows) - used)
        if unused:
            LOGGER.warning("%s: ids not used by any cell: %s.", name, unused)

    # -------------------------------------------------------------------------------------------- sections

    def coordinate_systems(self) -> None:
        defined = self.section("coordinate_systems")
        for cs_id, spec in defined.items():
            try:
                coordinate_system_from_dict(spec)
            except (ValueError, TypeError) as exc:
                self.problems.append(f"coordinate_systems[{cs_id}]: {exc}")
        ids = self.coordinate_system_ids
        if ids.size == 0:
            return
        negative = sorted(np.unique(ids[ids < 0]).tolist())
        if negative:
            self.problems.append(f"COORDINATE_SYSTEM: ids must be >= 0, found {negative}.")
        missing = sorted(set(np.unique(ids[ids > 0]).tolist()) - set(defined))
        if missing:
            self.problems.append(
                f"COORDINATE_SYSTEM: ids used by the mesh but not in data['coordinate_systems']: {missing}."
            )

    def property_rows(self, name: str, known: dict[str, str], what: str) -> dict[int, Any]:
        """The materials or constants rows: known names and value shapes."""
        rows = self.section(name)
        for row_id, row in rows.items():
            if not isinstance(row, Mapping):
                self.problems.append(f"{name}[{row_id}]: must be an object of {what} values.")
                continue
            for key, value in row.items():
                label = f"{name}[{row_id}].{key}"
                if key not in known:
                    self.problems.append(f"{label}: unknown {what}; known: {sorted(known)}.")
                else:
                    self.value_shape(value, known[key], label)
        return rows

    def element_types(self) -> dict[int, Any]:
        rows = self.section("element_types")
        if not rows:
            self.problems.append("element_types cannot be empty.")

        resolved: dict[int, dict[str, int]] = {}
        for type_id, row in rows.items():
            label = f"element_types[{type_id}]"
            if not isinstance(row, Mapping):
                self.problems.append(f"{label}: must be an object.")
                continue
            unknown = sorted(set(row) - _ELEMENT_TYPE_KEYS)
            if unknown:
                self.problems.append(f"{label}: unknown keys {unknown}.")
            missing = [key for key in ("element_type", "option1") if key not in row]
            if missing:
                self.problems.append(f"{label}: missing keys {missing}.")
                continue

            element_type = self.registry_index(row["element_type"], self.element_type_names)
            domain = self.registry_index(row["option1"], self.domains)
            if element_type is None:
                self.problems.append(
                    f"{label}: element_type must be an index below {len(self.element_type_names)} or one of"
                    f" {self.element_type_names}; got {row['element_type']!r}."
                )
            if domain is None:
                self.problems.append(
                    f"{label}: option1 (domain) must be an index below {len(self.domains)} or one of"
                    f" {self.domains}; got {row['option1']!r}."
                )
            order = row.get("option2", 1)
            if isinstance(order, bool) or order not in (1, 2):
                self.problems.append(f"{label}: option2 (FE order) must be 1 or 2; got {order!r}.")
            axisymmetric = row.get("option3", 0)
            if axisymmetric not in (0, 1):
                self.problems.append(f"{label}: option3 (axisymmetric) must be 0 or 1; got {axisymmetric!r}.")
            if element_type is not None and domain is not None:
                type_name = self.element_type_names[element_type]
                supported = element_library.element_types().get(type_name, {}).get("domains")
                if supported is not None and self.domains[domain] not in supported:
                    self.problems.append(
                        f"{label}: the domain {self.domains[domain]} (option1) is not one of the domains of"
                        f" {type_name}: {supported}."
                    )
            if element_type is not None and domain is not None and order in (1, 2) and axisymmetric in (0, 1):
                resolved[type_id] = {"element_type": element_type, "axisymmetric": int(axisymmetric)}

        self.mesh_ids("element_types", rows, self.type_ids, missing_is_error=True)
        self.cell_dimensions(resolved)
        return rows

    @staticmethod
    def registry_index(value: Any, names: list[str]) -> int | None:
        """Index of a name, or an integer index, in a registry; None when invalid."""
        if isinstance(value, str):
            return names.index(value) if value in names else None
        return int(value) if _is_index(value, len(names)) else None

    def cell_dimensions(self, resolved: dict[int, dict[str, int]]) -> None:
        """One element type and one option3 for the cells of the mesh, of the dimension of the cells; 3D cells not
        axisymmetric. Only the rows used by cells are checked: the model data may define more (unused) rows."""
        active = sorted(set(np.unique(self.type_ids).tolist()) & set(resolved))
        if not active:
            return
        used = {resolved[t]["element_type"] for t in active}
        if len(used) > 1:
            self.problems.append(
                "element_types: all the cells must have the same element_type, found "
                f"{ {t: self.element_type_names[resolved[t]['element_type']] for t in active} }."
            )
        if len({resolved[t]["axisymmetric"] for t in active}) > 1:
            self.problems.append(
                "element_types: the cells must be all planar or all axisymmetric (option3), found "
                f"{ {t: resolved[t]['axisymmetric'] for t in active} }."
            )
        dimensions = meshes.cell_dimensions(self.mesh)
        for type_id in active:
            spec = resolved[type_id]
            cells = dimensions[self.type_ids == type_id]
            name = self.element_type_names[spec["element_type"]]
            expected = self.element_type_dimensions[spec["element_type"]]
            if np.any(cells == 2) and expected != 2:
                self.problems.append(f"element_types[{type_id}]: 2D cells need a 2D element type, got {name!r}.")
            if np.any(cells == 3):
                if expected != 3:
                    self.problems.append(f"element_types[{type_id}]: 3D cells need a 3D element type, got {name!r}.")
                if spec["axisymmetric"] != 0:
                    self.problems.append(f"element_types[{type_id}]: option3 (axisymmetric) must be 0 for 3D cells.")

    def dof_couplings(self) -> None:
        """Each coupling a non-empty list of DOF fields; every CP id of the mesh defined."""
        fields = set(element_library.dof_fields(self.model.element_type_names() or None))
        couplings = self.section("dof_couplings")
        for cp_id, variables in couplings.items():
            label = f"dof_couplings[{cp_id}]"
            if not isinstance(variables, (list, tuple)) or not variables:
                self.problems.append(f"{label}: must be a non-empty list of DOF fields.")
                continue
            unknown = [v for v in variables if v not in fields]
            if unknown:
                self.problems.append(f"{label}: the DOF fields must be in {sorted(fields)}; got {unknown}.")

        used: set[int] = set()
        point_data = self.mesh.GetPointData()
        for i in range(point_data.GetNumberOfArrays()):
            name = point_data.GetArrayName(i)
            if name and _CP_ARRAY.fullmatch(name):
                values = np.asarray(meshes.array_values(self.mesh, name, "node"), dtype=float).reshape(-1)
                used.update(np.unique(values[np.isfinite(values) & (values > 0.0)].astype(int)).tolist())
        missing = sorted(used - set(couplings))
        if missing:
            self.problems.append(f"CP_k arrays: coupling ids not in data['dof_couplings']: {missing}.")

    def nonlinear_laws(self) -> None:
        """Each entry a law of nonlinear_laws.json with its data (the curve arrays: same length >= 2, x increasing),
        at most once per id; ids not used by cells of the law's source are warned about."""
        if "bh_curves" in self.data:
            self.problems.append(
                "bh_curves: the curves are given with their law: data['nonlinear_laws'][<MATERIAL id>] = "
                "{'type': 'bh_curve', 'B': [...], 'H': [...]}."
            )
        cell_ids = {"materials": self.material_ids, "constants": self.constants_ids, "element_types": self.type_ids}
        for row_id, value in self.section("nonlinear_laws").items():
            entries = value if isinstance(value, list) else [value]
            seen: set[str] = set()
            for n, entry in enumerate(entries):
                label = f"nonlinear_laws[{row_id}]" + (f"[{n}]" if isinstance(value, list) else "")
                if not isinstance(entry, Mapping) or "type" not in entry:
                    self.problems.append(f"{label}: must be an object with the 'type' of the law.")
                    continue
                law = entry["type"]
                name = law if isinstance(law, str) else None
                if _is_index(law, len(self.law_names)):
                    name = self.law_names[int(law)]
                if name not in self.laws:
                    self.problems.append(
                        f"{label}: {law!r} is not a law of nonlinear_laws.json; known: {self.law_names}."
                    )
                    continue
                if name in seen:
                    self.problems.append(f"{label}: the law '{name}' is given twice for this id.")
                seen.add(name)
                curve = self.laws[name].get("curve")
                expected = {"type", *(curve or ())}
                unknown = sorted(set(entry) - expected)
                if unknown:
                    self.problems.append(
                        f"{label}: unknown keys {unknown} for the law '{name}'; expected {sorted(expected)}."
                    )
                if curve:
                    self.curve(entry, curve[0], curve[1], label)
                source = self.laws[name]["source"]
                if cell_ids[source].size and row_id not in set(np.unique(cell_ids[source]).tolist()):
                    LOGGER.warning("%s: no cell has the %s id %d.", label, source.rstrip("s").replace("_", " "), row_id)

    def curve(self, entry: Mapping[str, Any], x_key: str, y_key: str, label: str) -> None:
        """The x and y arrays of a law: numeric, finite, same length >= 2, x strictly increasing."""
        if not all(isinstance(entry.get(key), list) for key in (x_key, y_key)):
            self.problems.append(f"{label}: needs the arrays '{x_key}' and '{y_key}'.")
            return
        try:
            x = np.asarray(entry[x_key], dtype=float)
            y = np.asarray(entry[y_key], dtype=float)
        except (TypeError, ValueError):
            self.problems.append(f"{label}: '{x_key}' and '{y_key}' must be numeric arrays.")
            return
        if x.ndim != 1 or y.shape != x.shape or x.size < 2:
            self.problems.append(f"{label}: '{x_key}' and '{y_key}' must be 1D arrays of the same length >= 2.")
        elif not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
            self.problems.append(f"{label}: the values must be finite.")
        elif np.any(np.diff(x) <= 0.0):
            self.problems.append(f"{label}: '{x_key}' must be strictly increasing.")


def check_model(model: Model) -> None:
    """Check the model: the arrays of its mesh and its data against the solver registries and the mesh.

    Warns about ids defined but not used by any cell.

    Raises:
        ValidationError: Listing every problem found.
        ImportError: If the solver module is not built.
    """
    check_required_arrays(model.mesh, model.element_type_names() or None)
    _ModelChecker(model).run()
