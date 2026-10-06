from __future__ import annotations

import os
import re
import math
from os import PathLike
from typing import Any, Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd


def ensure_vtk_unstructured_grid(mesh: Any, arg_name: str = "mesh"):
    """Validate that the input is a vtk.vtkUnstructuredGrid instance.

    Args:
        mesh: Object expected to be a vtk.vtkUnstructuredGrid.
        arg_name: Parameter name used in error messages.

    Returns:
        Any: The validated mesh object.
    """
    if mesh is None:
        raise ValueError(f"{arg_name} must not be None.")

    try:
        import vtk
    except Exception as exc:
        raise RuntimeError(
            "VTK is required to validate mesh types, but it could not be imported."
        ) from exc

    if not isinstance(mesh, vtk.vtkUnstructuredGrid):
        raise TypeError(
            f"{arg_name} must be a vtk.vtkUnstructuredGrid, got {type(mesh).__name__}."
        )

    return mesh


def ensure_non_empty_string(value: Any, arg_name: str) -> str:
    """Validate a non-empty string argument.

    Args:
        value: Value to validate.
        arg_name: Parameter name used in error messages.

    Returns:
        str: Normalized string value.
    """
    if not isinstance(value, str):
        raise TypeError(f"{arg_name} must be a string, got {type(value).__name__}.")
    if not value.strip():
        raise ValueError(f"{arg_name} must be a non-empty string.")
    return value


def ensure_pathlike(value: Union[str, PathLike[str]], arg_name: str) -> str:
    """Validate and normalize a path-like argument.

    Args:
        value: Path value.
        arg_name: Parameter name used in error messages.

    Returns:
        str: Normalized non-empty path string.
    """
    try:
        path = os.fspath(value)
    except TypeError as exc:
        raise TypeError(f"{arg_name} must be a str or PathLike path.") from exc
    if not isinstance(path, str) or not path.strip():
        raise ValueError(f"{arg_name} must be a non-empty path.")
    return path


def ensure_dataframe(
    value: Any,
    arg_name: str,
    required_columns: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Validate that a value is a pandas DataFrame with optional required columns.

    Args:
        value: Candidate dataframe.
        arg_name: Parameter name used in error messages.
        required_columns: Optional iterable of required column names.

    Returns:
        pd.DataFrame: Validated dataframe.
    """
    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{arg_name} must be a pandas.DataFrame, got {type(value).__name__}.")
    if required_columns is not None:
        missing = set(required_columns).difference(value.columns)
        if missing:
            raise ValueError(
                f"{arg_name} is missing required columns: {sorted(missing)}."
            )
    return value


def ensure_numpy_2d(value: Any, arg_name: str, cols: Optional[int] = None) -> np.ndarray:
    """Validate that a value is a 2D NumPy-compatible array.

    Args:
        value: Candidate array-like value.
        arg_name: Parameter name used in error messages.
        cols: Optional expected number of columns.

    Returns:
        np.ndarray: Validated 2D array.
    """
    arr = np.asarray(value)
    if arr.ndim != 2:
        raise ValueError(f"{arg_name} must be a 2D array, got {arr.ndim}D.")
    if cols is not None and arr.shape[1] != cols:
        raise ValueError(
            f"{arg_name} must have shape (N, {cols}), got {arr.shape}."
        )
    return arr


def ensure_in_options(value: Any, arg_name: str, options: Sequence[Any]) -> Any:
    """Validate that a value belongs to a fixed set of options.

    Args:
        value: Candidate value.
        arg_name: Parameter name used in error messages.
        options: Allowed values.

    Returns:
        Any: Original validated value.
    """
    if value not in options:
        raise ValueError(f"{arg_name} must be one of {list(options)}, got {value!r}.")
    return value


def validate_required_arrays_present(mesh: Any) -> None:
    """Ensure required VTK point/cell arrays are present with expected shapes."""
    ensure_vtk_unstructured_grid(mesh, arg_name="mesh")

    from .constants import VTK_ARRAYS_TO_INIT

    errors: list[str] = []
    cell_data = mesh.GetCellData()
    point_data = mesh.GetPointData()

    n_cells = mesh.GetNumberOfCells()
    n_points = mesh.GetNumberOfPoints()

    for required in ("TYPE", "MATERIAL", "CONSTANT", "COORDINATE_SYSTEM"):
        arr = cell_data.GetArray(required) if cell_data is not None else None
        if arr is None:
            errors.append(f"CellData array '{required}' is required but was not found.")

    for name, info in VTK_ARRAYS_TO_INIT.items():
        loc = info["location"]
        data_obj = point_data if loc == "node" else cell_data
        arr = data_obj.GetArray(name) if data_obj is not None else None
        if arr is None:
            errors.append(f"{loc} data array '{name}' is required but was not found.")
            continue

        expected_components = int(info["components"])
        if arr.GetNumberOfComponents() != expected_components:
            errors.append(
                f"{loc} data array '{name}' must have {expected_components} component(s),"
                f" got {arr.GetNumberOfComponents()}."
            )

        expected_tuples = n_points if loc == "node" else n_cells
        if arr.GetNumberOfTuples() != expected_tuples:
            errors.append(
                f"{loc} data array '{name}' must have {expected_tuples} tuple(s),"
                f" got {arr.GetNumberOfTuples()}."
            )

    if errors:
        raise ValueError(
            "Loaded mesh is missing required arrays or has invalid array shapes:\n"
            + "\n".join(f"  - {e}" for e in errors)
        )


# ---------------------------------------------------------------------------
# Solver-settings validation
# ---------------------------------------------------------------------------

#: Vacuum permeability [H/m], used as the BH-curve saturation-slope reference.
_MU0: float = 4.0 * math.pi * 1e-7

#: KSP types supported by PETSc / pyTEC3D (direct and iterative).
_VALID_KSP_TYPES: frozenset[str] = frozenset({
    "preonly",       # direct solve (pair with lu / cholesky)
    "gmres",
    "fgmres",
    "bcgs",
    "bcgsl",
    "cgs",
    "cg",
    "minres",
    "symmlq",
    "lgmres",
    "ibcgs",
    "gcr",
    "tfqmr",
    "tcqmr",
    "cr",
    "lsqr",
    "bicg",
    "richardson",
    "chebyshev",
})

#: PC types supported by PETSc / pytec3d.
_VALID_PC_TYPES: frozenset[str] = frozenset({
    "lu",            # direct LU  → MUMPS in parallel
    "cholesky",      # direct Cholesky → MUMPS in parallel
    "ilu",
    "icc",
    "asm",
    "gasm",
    "bjacobi",
    "jacobi",
    "sor",
    "eisenstat",
    "hypre",
    "ml",
    "gamg",
    "mg",
    "ksp",
    "none",
})


def validate_solver_settings(settings: dict[str, Any]) -> None:
    """Validate the solver_settings dict against allowed types and value ranges.

    Checks:
    - ``ksp_type``: string, one of the recognised PETSc KSP type names.
    - ``pc_type``: string, one of the recognised PETSc PC type names.
    - ``initial_linear_tolerance``: positive finite float (relative residual target).
    - ``minimum_linear_tolerance``: positive finite float, ≤ ``initial_linear_tolerance``.
    - ``max_linear_iterations``: integer ≥ 1.

    Args:
        settings: Dictionary of solver settings (typically ``solver.solver_settings``).

    Returns:
        None

    Raises:
        TypeError: When a value has the wrong Python type.
        ValueError: When a value is outside the allowed range or not a recognised option.
    """
    errors: list[str] = []

    # --- ksp_type -----------------------------------------------------------
    ksp = settings.get("ksp_type")
    if not isinstance(ksp, str):
        errors.append(
            f"ksp_type must be a string, got {type(ksp).__name__}."
        )
    elif ksp not in _VALID_KSP_TYPES:
        errors.append(
            f"ksp_type {ksp!r} is not recognised. "
            f"Valid options: {sorted(_VALID_KSP_TYPES)}."
        )

    # --- pc_type ------------------------------------------------------------
    pc = settings.get("pc_type")
    if not isinstance(pc, str):
        errors.append(
            f"pc_type must be a string, got {type(pc).__name__}."
        )
    elif pc not in _VALID_PC_TYPES:
        errors.append(
            f"pc_type {pc!r} is not recognised. "
            f"Valid options: {sorted(_VALID_PC_TYPES)}."
        )

    # --- initial_linear_tolerance -------------------------------------------
    rtol = settings.get("initial_linear_tolerance")
    if not isinstance(rtol, (int, float)):
        errors.append(
            f"initial_linear_tolerance must be a real number, got {type(rtol).__name__}."
        )
    elif not (np.isfinite(rtol) and rtol > 0.0):
        errors.append(
            f"initial_linear_tolerance must be a finite positive number, got {rtol!r}."
        )

    # --- minimum_linear_tolerance -------------------------------------------
    abstol = settings.get("minimum_linear_tolerance")
    if not isinstance(abstol, (int, float)):
        errors.append(
            f"minimum_linear_tolerance must be a real number, got {type(abstol).__name__}."
        )
    elif not (np.isfinite(abstol) and abstol > 0.0):
        errors.append(
            f"minimum_linear_tolerance must be a finite positive number, got {abstol!r}."
        )
    elif isinstance(rtol, (int, float)) and np.isfinite(rtol) and abstol > rtol:
        errors.append(
            f"minimum_linear_tolerance ({abstol!r}) must be ≤ "
            f"initial_linear_tolerance ({rtol!r})."
        )

    # --- max_linear_iterations ----------------------------------------------
    maxits = settings.get("max_linear_iterations")
    if not isinstance(maxits, int) or isinstance(maxits, bool):
        errors.append(
            f"max_linear_iterations must be an integer, got {type(maxits).__name__}."
        )
    elif maxits < 1:
        errors.append(
            f"max_linear_iterations must be ≥ 1, got {maxits!r}."
        )

    if errors:
        raise ValueError(
            "Invalid solver_settings:\n" + "\n".join(f"  - {e}" for e in errors)
        )


# ---------------------------------------------------------------------------
# VTK cell-type classification arrays (sorted for np.isin, built at import)
# ---------------------------------------------------------------------------

#: 2-D surface cell types (VTK dimension = 2).
_VTK_2D_CELL_TYPES: np.ndarray = np.array(sorted([
    5, 6, 7, 8, 9, 22, 23, 28, 30, 34, 36,
    52, 53, 54, 61, 62, 63, 69, 70, 76, 77,
]), dtype=np.int32)

#: 3-D volumetric cell types (VTK dimension = 3).
_VTK_3D_CELL_TYPES: np.ndarray = np.array(sorted([
    10, 11, 12, 13, 14, 15, 16, 24, 25, 26, 27, 29,
    31, 32, 33, 37, 41, 42, 55, 56, 64, 65, 66, 67,
    71, 72, 73, 74, 78, 79,
]), dtype=np.int32)

#: Second-order (quadratic / higher-order) cell types.
_VTK_QUADRATIC_CELL_TYPES: np.ndarray = np.array(sorted([
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 61, 62, 63, 64, 65, 66, 67,
    69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
]), dtype=np.int32)


def _solver_registries() -> dict[str, Any]:
    """Domains, fields, element types, material properties and real constants compiled into the solver."""
    import importlib

    try:
        module = importlib.import_module("pytec3d_solver.pytec3d_solver")
    except ImportError as exc:
        raise ImportError(
            "Could not import the C++ solver module 'pytec3d_solver.pytec3d_solver', which holds the"
            " registries the model data is validated against. Build the extension first."
        ) from exc
    return module.registries()


def _section_entries(data: Any, section: str, errors: list[str]) -> dict[int, Any]:
    """Return ``data[section]`` keyed by int id; ids that are not positive integers are reported."""
    entries = data.get(section, {})
    if not isinstance(entries, dict):
        errors.append(f"data['{section}'] must be an object, got {type(entries).__name__}.")
        return {}
    out: dict[int, Any] = {}
    for key, value in entries.items():
        try:
            entry_id = int(key)
            valid = entry_id > 0 and str(entry_id) == str(key).strip()
        except (TypeError, ValueError):
            valid = False
        if not valid:
            errors.append(f"{section}: the id {key!r} must be a positive integer.")
            continue
        if entry_id in out:
            errors.append(f"{section}: the id {entry_id} is given more than once.")
            continue
        out[entry_id] = value
    return out


def _check_value_shape(value: Any, value_type: str, label: str, errors: list[str]) -> Optional[np.ndarray]:
    """Check a number, a 3-vector (a number stands for equal components) or a 3x3 tensor; return it as an array."""
    shapes = {"scalar": [()], "vector": [(), (3,)], "tensor": [(3, 3)]}[value_type]
    try:
        arr = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        errors.append(f"{label}: must be numeric.")
        return None
    if arr.shape not in shapes:
        expected = {"scalar": "a number", "vector": "a number or 3 numbers", "tensor": "a 3x3 array"}[value_type]
        errors.append(f"{label}: must be {expected}, got shape {arr.shape}.")
        return None
    if not np.all(np.isfinite(arr)):
        errors.append(f"{label}: must be finite.")
        return None
    return arr


def validate_pytec3dodel(model: Any) -> None:
    """Perform a comprehensive consistency check on a loaded pyTEC3D model.

    Validates ``model.data`` (``element_types``, ``materials``, ``constants``, ``dof_couplings``,
    ``bh_curves``, ``coordinate_systems``) against the registries compiled into the solver and
    against the VTK mesh stored in *model*. All violations are collected before raising so that a
    single call reports every problem found. Unused identifiers are emitted as log warnings.

    A material property or real constant may be left undefined. Which ones a kernel needs depends
    on the element type and domain; the solver checks that when the model is loaded. Here, the values
    that are defined are checked, and the real constants a stranded coil always needs.

    Args:
        model: ``pyTEC3D_Model`` instance whose mesh and data are to be validated.

    Returns:
        None

    Raises:
        ValueError: When one or more consistency constraints are violated.
    """
    import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]
    from .coordinate_system import coordinate_system_from_dict
    from .logging_utils import get_logger

    _log = get_logger("pytec3d.validation")
    errors: list[str] = []

    mesh = model.mesh
    if mesh is None:
        raise ValueError("No mesh loaded; call read_mesh() before validating.")

    # Keep this helper centralized so mesh-load and full model validation
    # enforce the exact same array contract.
    validate_required_arrays_present(mesh)

    registries = _solver_registries()
    domains: list[str] = list(registries["domains"])
    element_type_names: list[str] = list(registries["element_types"])
    element_type_dimensions: list[int] = [int(d) for d in registries["element_types"].values()]
    material_properties: dict[str, str] = dict(registries["material_properties"])
    real_constants: dict[str, str] = dict(registries["real_constants"])

    data = model.data
    point_data = mesh.GetPointData()
    cell_data = mesh.GetCellData()

    # ── per-cell VTK type classification ────────────────────────────────────
    vtk_ct = np.asarray(vtknp.vtk_to_numpy(mesh.GetCellTypesArray()), dtype=np.int32)
    is_2d = np.isin(vtk_ct, _VTK_2D_CELL_TYPES)
    is_3d = np.isin(vtk_ct, _VTK_3D_CELL_TYPES)

    def _get_int_array(name: str) -> np.ndarray:
        arr = cell_data.GetArray(name)
        return np.asarray(vtknp.vtk_to_numpy(arr), dtype=np.int64).reshape(-1)

    type_ids = _get_int_array("TYPE")
    mat_ids = _get_int_array("MATERIAL")
    const_ids = _get_int_array("CONSTANT")
    coordinate_system_ids = _get_int_array("COORDINATE_SYSTEM")

    # ========================================================================
    # 1. coordinate_systems
    # ========================================================================
    coordinate_systems = _section_entries(data, "coordinate_systems", errors)
    for cs_id, spec in coordinate_systems.items():
        try:
            coordinate_system_from_dict(spec)
        except (ValueError, TypeError) as exc:
            errors.append(f"coordinate_systems[{cs_id}]: {exc}")

    if coordinate_system_ids.size > 0:
        bad_cs = sorted(np.unique(coordinate_system_ids[coordinate_system_ids < 0]).tolist())
        if bad_cs:
            errors.append(f"COORDINATE_SYSTEM: values must be >= 0; found invalid values: {bad_cs}.")
        missing_cs = sorted(
            set(np.unique(coordinate_system_ids[coordinate_system_ids > 0]).tolist()) - set(coordinate_systems)
        )
        if missing_cs:
            errors.append(
                "COORDINATE_SYSTEM: positive IDs used by the mesh are not defined in"
                f" data['coordinate_systems']: {missing_cs}."
            )

    # ========================================================================
    # 2. materials
    # ========================================================================
    materials = _section_entries(data, "materials", errors)
    for material_id, material in materials.items():
        if not isinstance(material, dict):
            errors.append(f"materials[{material_id}]: must be an object of material properties.")
            continue
        for name, value in material.items():
            label = f"materials[{material_id}].{name}"
            if name not in material_properties:
                errors.append(f"{label}: unknown material property; known: {sorted(material_properties)}.")
                continue
            arr = _check_value_shape(value, material_properties[name], label, errors)
            if arr is None:
                continue
            # mur and rho are tensors; their diagonal is what has a sign.
            diagonal = np.diagonal(arr) if arr.ndim == 2 else arr.reshape(-1)
            if name == "mur" and not np.all(diagonal > 0.0):
                errors.append(f"{label}: the diagonal must be positive.")
            if name == "rho" and not np.all(diagonal >= 0.0):
                errors.append(f"{label}: the diagonal must be non-negative.")

    if mat_ids.size > 0:
        mesh_materials = set(np.unique(mat_ids).tolist())
        missing_m = sorted(mesh_materials - set(materials))
        if missing_m:
            errors.append(f"materials: cell MATERIAL values absent from data['materials']: {missing_m}.")
        unused_m = sorted(set(materials) - mesh_materials)
        if unused_m:
            _log.warning("materials: ids not used by any cell: %s.", unused_m)

    # ========================================================================
    # 3. bh_curves
    # ========================================================================
    bh_curves = _section_entries(data, "bh_curves", errors)
    mesh_material_set = set(np.unique(mat_ids).tolist()) if mat_ids.size > 0 else set()
    for material_id, curve in bh_curves.items():
        label = f"bh_curves[{material_id}]"
        if mesh_material_set and material_id not in mesh_material_set:
            errors.append(f"{label}: material ID {material_id} is not present in mesh MATERIAL cell values.")
        if not isinstance(curve, dict) or set(curve) != {"H", "B"}:
            errors.append(f"{label}: must be an object with the arrays 'H' and 'B'.")
            continue
        try:
            h_values = np.asarray(curve["H"], dtype=float)
            b_values = np.asarray(curve["B"], dtype=float)
        except (TypeError, ValueError):
            errors.append(f"{label}: 'H' and 'B' must be numeric arrays.")
            continue
        if h_values.ndim != 1 or b_values.shape != h_values.shape or h_values.size < 2:
            errors.append(f"{label}: 'H' and 'B' must be 1D arrays of the same length >= 2.")
            continue
        if not np.all(np.isfinite(h_values)) or not np.all(np.isfinite(b_values)):
            errors.append(f"{label}: H and B values must be finite.")
            continue

        if h_values[0] != 0.0:
            errors.append(f"{label}: H must start at 0.0, got {h_values[0]!r}.")
        if b_values[0] != 0.0:
            errors.append(f"{label}: B must start at 0.0, got {b_values[0]!r}.")
        if np.any(h_values < 0.0) or np.any(b_values < 0.0):
            errors.append(f"{label}: H and B values must be non-negative.")

        dh = np.diff(h_values)
        db = np.diff(b_values)
        if np.any(dh <= 0.0):
            errors.append(f"{label}: H values must be strictly monotonic increasing.")
        if np.any(db <= 0.0):
            errors.append(f"{label}: B values must be strictly monotonic increasing.")
            continue

        # The C++ solver floors the extrapolated slope beyond b_values[-1] at
        # 1/mu0 (vacuum permeability). If the curve's own last segment is far
        # softer than that, dH/dB will visibly kink right at the last point.
        last_slope = dh[-1] / db[-1]
        if last_slope < 0.5 / _MU0:
            _log.warning(
                "%s: last segment slope (dH/dB=%.3g) is well below "
                "vacuum reluctivity (1/mu0=%.3g); extrapolation beyond B=%.4g will "
                "show a slope discontinuity there. Consider extending the table so "
                "it saturates toward the vacuum slope.",
                label, last_slope, 1.0 / _MU0, b_values[-1],
            )

        # mur is ignored by the solver for materials with a bh_curve (nonlinear
        # elements always use the B-H table instead of a constant permeability),
        # so a mur != 1 on those materials likely indicates an unintended input.
        material = materials.get(material_id)
        if isinstance(material, dict) and "mur" in material:
            mur = np.asarray(material["mur"], dtype=float)
            if mur.shape != (3, 3) or np.any(np.abs(mur - np.eye(3)) > 1e-9):
                _log.warning(
                    "materials[%d]: a bh_curve is defined, so mur=%s is ignored by the solver"
                    " (nonlinear elements use the B-H table only).",
                    material_id, mur.tolist(),
                )

    # ========================================================================
    # 4. element_types
    # ========================================================================
    element_types = _section_entries(data, "element_types", errors)
    if not element_types:
        errors.append("element_types cannot be empty.")

    def _registry_index(value: Any, names: list[str]) -> Optional[int]:
        """Index in a registry of a name or of an integer index; None when invalid."""
        if isinstance(value, str):
            return names.index(value) if value in names else None
        if isinstance(value, (int, np.integer)) and not isinstance(value, bool) and 0 <= int(value) < len(names):
            return int(value)
        return None

    resolved: dict[int, dict[str, int]] = {}
    for type_id, entry in element_types.items():
        label = f"element_types[{type_id}]"
        if not isinstance(entry, dict):
            errors.append(f"{label}: must be an object.")
            continue
        unknown = sorted(set(entry) - {"element_type", "option1", "option2", "option3", "option4", "option5"})
        if unknown:
            errors.append(f"{label}: unknown keys {unknown}.")
        missing = [key for key in ("element_type", "option1") if key not in entry]
        if missing:
            errors.append(f"{label}: missing keys {missing}.")
            continue

        element_type = _registry_index(entry["element_type"], element_type_names)
        domain = _registry_index(entry["option1"], domains)
        if element_type is None:
            errors.append(
                f"{label}: element_type must be an index below {len(element_type_names)} or one of"
                f" {element_type_names}; got {entry['element_type']!r}."
            )
        if domain is None:
            errors.append(
                f"{label}: option1 (domain) must be an index below {len(domains)} or one of"
                f" {domains}; got {entry['option1']!r}."
            )
        order = entry.get("option2", 1)
        if order not in (1, 2) or isinstance(order, bool):
            errors.append(f"{label}: option2 (FE order) must be 1 or 2; got {order!r}.")
        axisymmetric = entry.get("option3", 0)
        if axisymmetric not in (0, 1):
            errors.append(f"{label}: option3 (axisymmetric) must be 0 or 1; got {axisymmetric!r}.")
        if entry.get("option4", 0) not in (0, 1):
            errors.append(f"{label}: option4 must be 0 or 1; got {entry.get('option4')!r}.")
        if element_type is not None and domain is not None and order in (1, 2) and axisymmetric in (0, 1):
            resolved[type_id] = {
                "element_type": element_type,
                "domain": domain,
                "order": int(order),
                "axisymmetric": int(axisymmetric),
            }

    if resolved:
        axisymmetric_values = sorted({spec["axisymmetric"] for spec in resolved.values()})
        if len(axisymmetric_values) > 1:
            errors.append(
                "element_types: option3 (axisymmetric) cannot contain distinct values; "
                "it must be globally all 0 or all 1."
                f" Found: { {i: spec['axisymmetric'] for i, spec in resolved.items()} }."
            )

    mesh_type_set = set(np.unique(type_ids).tolist()) if type_ids.size > 0 else set()
    missing_t = sorted(mesh_type_set - set(element_types))
    if missing_t:
        errors.append(f"element_types: cell TYPE values absent from data['element_types']: {missing_t}.")
    unused_t = sorted(set(element_types) - mesh_type_set)
    if unused_t:
        _log.warning("element_types: TYPE ids not used by any cell: %s.", unused_t)

    active_types = sorted(mesh_type_set & set(resolved))
    if active_types:
        used_element_types = {resolved[t]["element_type"] for t in active_types}
        if len(used_element_types) > 1:
            errors.append(
                "element_types: all cells in the mesh must share the same element_type,"
                f" but found: { {t: element_type_names[resolved[t]['element_type']] for t in active_types} }."
            )

    coil_domain = domains.index("STRANDED_COIL") if "STRANDED_COIL" in domains else -1
    cell_is_coil = np.zeros(type_ids.shape, dtype=bool)
    cell_coil_dimension = np.zeros(type_ids.shape, dtype=np.int8)
    cell_coil_axisymmetric = np.zeros(type_ids.shape, dtype=bool)
    for type_id in active_types:
        spec = resolved[type_id]
        mask = type_ids == type_id
        dimension = element_type_dimensions[spec["element_type"]]
        if bool(np.any(is_2d[mask])) and dimension != 2:
            errors.append(f"element_types[{type_id}]: 2D cells need a 2D element type, got {element_type_names[spec['element_type']]!r}.")
        if bool(np.any(is_3d[mask])):
            if dimension != 3:
                errors.append(f"element_types[{type_id}]: 3D cells need a 3D element type, got {element_type_names[spec['element_type']]!r}.")
            if spec["axisymmetric"] != 0:
                errors.append(f"element_types[{type_id}]: option3 (axisymmetric) must be 0 for 3D cells.")
        if spec["domain"] == coil_domain:
            cell_is_coil |= mask
            cell_coil_dimension[mask] = dimension
            cell_coil_axisymmetric[mask] = bool(spec["axisymmetric"])

    # ========================================================================
    # 5. dof_couplings
    # ========================================================================
    _valid_vars = {"A", "Az", "VOLT", "V", "EMF"}
    dof_couplings = _section_entries(data, "dof_couplings", errors)
    for cp_id, variables in dof_couplings.items():
        label = f"dof_couplings[{cp_id}]"
        if not isinstance(variables, (list, tuple)) or not variables:
            errors.append(f"{label}: must be a non-empty list of variable names.")
            continue
        bad_vars = [v for v in variables if v not in _valid_vars]
        if bad_vars:
            errors.append(f"{label}: variables must be in {sorted(_valid_vars)}; got {bad_vars}.")

    # CP_xx arrays: all values > 0 must exist in dof_couplings.
    cp_name = re.compile(r"^CP_(\d+)$")
    cp_ids_in_mesh: set[int] = set()
    if point_data is not None:
        for i in range(point_data.GetNumberOfArrays()):
            arr = point_data.GetArray(i)
            if arr is None:
                continue
            name = arr.GetName()
            if not isinstance(name, str) or cp_name.fullmatch(name) is None:
                continue
            vals = np.asarray(vtknp.vtk_to_numpy(arr), dtype=float).reshape(-1)
            pos = vals[np.isfinite(vals) & (vals > 0.0)]
            if pos.size > 0:
                cp_ids_in_mesh.update(np.unique(pos.astype(int)).tolist())
    if cp_ids_in_mesh:
        missing_cp_ids = sorted(cp_ids_in_mesh - set(dof_couplings))
        if missing_cp_ids:
            errors.append(
                "CP_xx arrays contain positive coupling IDs not found in data['dof_couplings']:"
                f" {missing_cp_ids}."
            )

    # ========================================================================
    # 6. constants
    # ========================================================================
    constants = _section_entries(data, "constants", errors)
    constant_arrays: dict[int, dict[str, np.ndarray]] = {}
    for constants_id, entry in constants.items():
        if not isinstance(entry, dict):
            errors.append(f"constants[{constants_id}]: must be an object of real constants.")
            continue
        arrays: dict[str, np.ndarray] = {}
        for name, value in entry.items():
            label = f"constants[{constants_id}].{name}"
            if name not in real_constants:
                errors.append(f"{label}: unknown real constant; known: {sorted(real_constants)}.")
                continue
            arr = _check_value_shape(value, real_constants[name], label, errors)
            if arr is None:
                continue
            arrays[name] = arr
            if name in ("turns", "cross_section", "volume", "radius", "resistance", "thickness") and not float(arr) > 0.0:
                errors.append(f"{label}: must be > 0; got {float(arr)}.")
            if name == "symmetry" and not float(arr) >= 1.0:
                errors.append(f"{label}: must be >= 1; got {float(arr)}.")
            if name == "T":
                t = np.broadcast_to(arr, (3,))
                if not np.all(np.abs(t) <= 1.0) or abs(float(np.linalg.norm(t)) - 1.0) > 1e-4:
                    errors.append(f"{label}: must be a unit vector; got {t.tolist()}.")
        constant_arrays[constants_id] = arrays

    if const_ids.size > 0:
        mesh_constants = set(np.unique(const_ids).tolist())
        unused_c = sorted(set(constants) - mesh_constants)
        if unused_c:
            _log.warning("constants: ids not used by any cell: %s.", unused_c)
        absent_c = sorted(c for c in mesh_constants - set(constants) if c > 0)
        if absent_c:
            _log.warning("constants: cell CONSTANT values absent from data['constants']: %s.", absent_c)

    # Real constants a stranded coil always needs, by cell group.
    coil_constants_ids = np.unique(const_ids[cell_is_coil]) if cell_is_coil.any() else np.empty((0,), dtype=np.int64)
    for constants_id in coil_constants_ids.tolist():
        group = cell_is_coil & (const_ids == constants_id)
        dimension = int(cell_coil_dimension[group][0])
        axisymmetric = bool(cell_coil_axisymmetric[group][0])
        needed = ["turns", "cross_section", "resistance", "symmetry", "T"]
        needed += ["volume"] if dimension == 3 else ["thickness"] + (["radius"] if axisymmetric else [])
        arrays = constant_arrays.get(int(constants_id), {})
        absent = [name for name in needed if name not in arrays]
        if absent:
            errors.append(
                f"constants[{int(constants_id)}]: a stranded coil ({dimension}D) needs {absent}."
            )
        elif dimension == 2 and abs(abs(float(np.broadcast_to(arrays["T"], (3,))[2])) - 1.0) > 1e-9:
            errors.append(f"constants[{int(constants_id)}]: T must be (0, 0, -1) or (0, 0, 1) for a 2D coil.")
        if axisymmetric and "thickness" in arrays and float(arrays["thickness"]) > 1.0:
            errors.append(f"constants[{int(constants_id)}]: thickness must be <= 1 for axisymmetric cells.")

    if errors:
        raise ValueError(
            "pyTEC3D model validation failed with {} error(s):\n  - {}".format(
                len(errors), "\n  - ".join(errors)
            )
        )

    _log.info("Validated model data: %d element types, %d materials, %d real constants.",
              len(element_types), len(materials), len(constants))


def validate_pre_solve_inputs(model: Any, time_step: Any) -> None:
    """Validate runtime inputs right before solve and emit actionable warnings.

    This check is intended to run on every solve call so users get immediate
    feedback if loads/BCs/time-step settings would lead to a no-op or unstable
    run setup.

    Args:
        model: ``pyTEC3D_Model`` with the mesh to solve.
        time_step: Time-step dict of the solver (``time``, ``theta``, ``transient_effects``).
    """
    from .logging_utils import get_logger

    import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

    _log = get_logger("pytec3d.validation")
    mesh = model.mesh

    if mesh is None:
        raise ValueError("No mesh loaded; call read_mesh() before solving.")

    point_data = mesh.GetPointData()
    cell_data = mesh.GetCellData()

    errors: list[str] = []

    def _get_array_values(data_obj: Any, name: str) -> np.ndarray:
        arr = data_obj.GetArray(name) if data_obj is not None else None
        if arr is None:
            errors.append(f"Required array '{name}' was not found before solve.")
            return np.empty((0,), dtype=float)
        return np.asarray(vtknp.vtk_to_numpy(arr), dtype=float)

    # 1) Loads check: Js (elem), AMPS/CURR (node)
    js = _get_array_values(cell_data, "Js")
    amps = _get_array_values(point_data, "AMPS")
    curr = _get_array_values(point_data, "CURR")
    if js.size > 0 and amps.size > 0 and curr.size > 0:
        if np.all(np.nan_to_num(js, nan=0.0) == 0.0) and np.all(np.nan_to_num(amps, nan=0.0) == 0.0) and np.all(np.nan_to_num(curr, nan=0.0) == 0.0):
            _log.warning(
                "No loads detected: arrays Js, AMPS, and CURR are all zero."
            )

    # 2) Boundary conditions check: all BC arrays all-NaN
    bc_names = ["_A", "_Az", "_VOLT", "_V", "_EMF"]
    bc_arrays = [_get_array_values(point_data, name) for name in bc_names]
    if all(arr.size > 0 for arr in bc_arrays):
        if all(np.all(np.isnan(arr)) for arr in bc_arrays):
            _log.warning(
                "No boundary conditions detected: arrays _A, _Az, _VOLT, _V, and _EMF are all NaN."
            )

    # 3) Time-step constraints
    if not isinstance(time_step, dict):
        errors.append("time_step must be a dictionary with keys 'time' and 'theta'.")
    else:
        t = time_step.get("time")
        th = time_step.get("theta")

        if not isinstance(t, (int, float)) or not np.isfinite(t) or float(t) <= 0.0:
            errors.append(f"time_step['time'] must be a finite value > 0, got {t!r}.")

        if not isinstance(th, (int, float)) or not np.isfinite(th) or not (0.0 <= float(th) <= 1.0):
            errors.append(
                f"time_step['theta'] must be a finite value in [0, 1], got {th!r}."
            )

    if errors:
        raise ValueError(
            "Invalid pre-solve inputs:\n" + "\n".join(f"  - {e}" for e in errors)
        )

    validate_required_arrays_present(mesh)

