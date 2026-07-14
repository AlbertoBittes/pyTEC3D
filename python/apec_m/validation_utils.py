from __future__ import annotations

import os
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


# ---------------------------------------------------------------------------
# Solver-settings validation
# ---------------------------------------------------------------------------

#: KSP types supported by PETSc / APEC-M (direct and iterative).
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

#: PC types supported by PETSc / APEC-M.
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
        settings: Dictionary of solver settings (typically ``apec.solver_settings``).

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


def validate_apec_model(apec: Any) -> None:
    """Perform a comprehensive consistency check on a loaded APEC model.

    Validates the ``element_types``, ``dof_couplings``, ``materials``, and
    ``constants`` DataFrames against each other and against the VTK mesh stored
    in *apec*.  All constraint violations are collected before raising so that
    a single call reports every problem found.  Unused identifiers are emitted
    as log warnings rather than errors.

    Args:
        apec: ``APEC_M`` instance whose mesh and metadata DataFrames are to
            be validated.

    Returns:
        None

    Raises:
        ValueError: When one or more consistency constraints are violated.
    """
    import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]
    from .logging_utils import get_logger
    from .constants import (
        ELEMENT_TYPES_COLUMNS_TYPES,
        MATERIALS_COLUMNS_TYPES,
        CONSTANTS_COLUMNS_TYPES,
        DOF_COUPLINGS_COLUMNS_TYPES,
    )

    _log = get_logger("apec_m.validation")
    errors: list[str] = []

    mesh = apec._mesh
    if mesh is None:
        raise ValueError("No mesh loaded; call read_mesh() before validating.")

    et = apec.element_types
    mat = apec.materials
    const = apec.constants
    dof = apec.dof_couplings

    # ── solver settings ─────────────────────────────────────────────────────
    validate_solver_settings(apec.solver_settings)

    n_cells: int = mesh.GetNumberOfCells()
    cell_data = mesh.GetCellData()

    # ── per-cell VTK type classification ────────────────────────────────────
    vtk_ct = np.asarray(vtknp.vtk_to_numpy(mesh.GetCellTypesArray()), dtype=np.int32)
    is_2d = np.isin(vtk_ct, _VTK_2D_CELL_TYPES)
    is_3d = np.isin(vtk_ct, _VTK_3D_CELL_TYPES)
    is_quad = np.isin(vtk_ct, _VTK_QUADRATIC_CELL_TYPES)

    def _get_int_array(name: str) -> Optional[np.ndarray]:
        arr = cell_data.GetArray(name)
        if arr is None:
            _log.warning(
                "CellData array '%s' not found; related checks skipped.", name
            )
            return None
        return np.asarray(vtknp.vtk_to_numpy(arr), dtype=np.int32)

    type_ids = _get_int_array("TYPE")
    mat_ids = _get_int_array("MATERIAL")
    const_ids = _get_int_array("CONSTANT")

    # ========================================================================
    # 0. bh_curves dictionary
    # ========================================================================
    bh_curves = getattr(apec, "bh_curves", {})
    if bh_curves is None:
        bh_curves = {}

    if not isinstance(bh_curves, dict):
        errors.append(f"bh_curves must be a dict, got {type(bh_curves).__name__}.")
    else:
        material_ids_in_grid: set[int] = set()
        if mat_ids is not None:
            material_ids_in_grid = set(np.unique(mat_ids).astype(int).tolist())

        for material_id_raw, component_map in bh_curves.items():
            if isinstance(material_id_raw, bool) or not isinstance(material_id_raw, (int, np.integer)):
                errors.append(
                    "bh_curves: material key must be an integer material ID; "
                    f"got key={material_id_raw!r} ({type(material_id_raw).__name__})."
                )
                continue

            material_id = int(material_id_raw)
            if material_id <= 0:
                errors.append(f"bh_curves: material key must be > 0, got {material_id}.")
                continue

            if material_ids_in_grid and material_id not in material_ids_in_grid:
                errors.append(
                    f"bh_curves: material ID {material_id} is not present in mesh MATERIAL cell values."
                )

            if not isinstance(component_map, dict):
                errors.append(
                    f"bh_curves[{material_id}]: value must be a dict indexed by component 0/1/2."
                )
                continue

            missing_components = {0, 1, 2} - {int(k) for k in component_map.keys() if isinstance(k, (int, np.integer)) and not isinstance(k, bool)}
            if missing_components:
                errors.append(
                    f"bh_curves[{material_id}]: missing component curves for {sorted(missing_components)}."
                )

            for component_raw, curve_values in component_map.items():
                if isinstance(component_raw, bool) or not isinstance(component_raw, (int, np.integer)):
                    errors.append(
                        f"bh_curves[{material_id}]: component key must be integer 0/1/2, got {component_raw!r}."
                    )
                    continue

                component = int(component_raw)
                if component not in (0, 1, 2):
                    errors.append(
                        f"bh_curves[{material_id}]: component key must be 0, 1, or 2, got {component}."
                    )
                    continue

                try:
                    curve = np.asarray(curve_values, dtype=float)
                except Exception:
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: could not be converted to a float NumPy array."
                    )
                    continue

                if curve.ndim != 2 or curve.shape[0] != 2 or curve.shape[1] < 2:
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: expected shape (2, N) with N >= 2, got {curve.shape}."
                    )
                    continue

                h_values = curve[0, :]
                b_values = curve[1, :]

                if not np.all(np.isfinite(h_values)) or not np.all(np.isfinite(b_values)):
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: H and B values must be finite."
                    )
                    continue

                if h_values[0] != 0.0:
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: H must start at 0.0, got {h_values[0]!r}."
                    )
                if b_values[0] != 0.0:
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: B must start at 0.0, got {b_values[0]!r}."
                    )

                if np.any(h_values < 0.0) or np.any(b_values < 0.0):
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: H and B values must be non-negative."
                    )

                dh = np.diff(h_values)
                db = np.diff(b_values)
                if np.any(dh <= 0.0):
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: H values must be strictly monotonic increasing."
                    )
                if np.any(db <= 0.0):
                    errors.append(
                        f"bh_curves[{material_id}][{component}]: B values must be strictly monotonic increasing."
                    )

    # ========================================================================
    # 1. element_types DataFrame
    # ========================================================================
    req_et = list(ELEMENT_TYPES_COLUMNS_TYPES.keys())
    miss_et = set(req_et) - set(et.columns)
    et_ok = not bool(miss_et)
    if miss_et:
        errors.append(f"element_types missing columns: {sorted(miss_et)}.")

    if et_ok and not et.empty:
        # TYPE_ID must be a positive integer
        bad = et.loc[~(et["TYPE_ID"] > 0), "TYPE_ID"].tolist()
        if bad:
            errors.append(f"element_types: TYPE_ID must be > 0; bad values: {bad}.")

        # No duplicate TYPE_IDs
        dups = et.loc[et["TYPE_ID"].duplicated(keep=False), "TYPE_ID"].unique().tolist()
        if dups:
            errors.append(f"element_types: duplicate TYPE_ID values: {dups}.")

        def _et_pairs(mask: Any, col: str) -> list:
            return list(zip(
                et.loc[mask, "TYPE_ID"].astype(int).tolist(),
                et.loc[mask, col].tolist(),
            ))

        # ELEMENT_TYPE ∈ {0, 1, 2}
        mask = ~et["ELEMENT_TYPE"].isin((0, 1, 2))
        if mask.any():
            errors.append(
                f"element_types: ELEMENT_TYPE must be 0, 1 or 2;"
                f" bad (TYPE_ID, value) pairs: {_et_pairs(mask, 'ELEMENT_TYPE')}."
            )

        # OPTION1 ∈ [0, 6]
        mask = ~et["OPTION1"].between(0, 6)
        if mask.any():
            errors.append(
                f"element_types: OPTION1 must be in [0, 6];"
                f" bad (TYPE_ID, value) pairs: {_et_pairs(mask, 'OPTION1')}."
            )

        # OPTION2 ∈ {1, 2}
        mask = ~et["OPTION2"].isin((1, 2))
        if mask.any():
            errors.append(
                f"element_types: OPTION2 must be 1 or 2;"
                f" bad (TYPE_ID, value) pairs: {_et_pairs(mask, 'OPTION2')}."
            )

        # OPTION3 ∈ {0, 1}
        mask = ~et["OPTION3"].isin((0, 1))
        if mask.any():
            errors.append(
                f"element_types: OPTION3 must be 0 or 1;"
                f" bad (TYPE_ID, value) pairs: {_et_pairs(mask, 'OPTION3')}."
            )

        # OPTION4 ∈ {0, 1}
        mask = ~et["OPTION4"].isin((0, 1))
        if mask.any():
            errors.append(
                f"element_types: OPTION4 must be 0 or 1;"
                f" bad (TYPE_ID, value) pairs: {_et_pairs(mask, 'OPTION4')}."
            )

        # ── cross-validate against mesh ──────────────────────────────────────
        if type_ids is not None:
            et_id_set: set[int] = set(et["TYPE_ID"].astype(int))
            mesh_type_set: set[int] = set(np.unique(type_ids).tolist())

            miss_t = mesh_type_set - et_id_set
            if miss_t:
                errors.append(
                    f"element_types: cell TYPE values absent from TYPE_ID: {sorted(miss_t)}."
                )

            unused_t = et_id_set - mesh_type_set
            if unused_t:
                _log.warning(
                    "element_types: TYPE_IDs not used by any cell: %s.", sorted(unused_t)
                )

            # All TYPE_IDs used by cells must share the same ELEMENT_TYPE value.
            # Only the TYPE_IDs that are actually defined in element_types matter here.
            active_tids = mesh_type_set & et_id_set
            if active_tids:
                active_et_vals = (
                    et.loc[et["TYPE_ID"].isin(active_tids), ["TYPE_ID", "ELEMENT_TYPE"]]
                    .set_index("TYPE_ID")["ELEMENT_TYPE"]
                    .astype(int)
                )
                unique_et_vals = active_et_vals.unique().tolist()
                if len(unique_et_vals) > 1:
                    detail = {
                        int(tid): int(active_et_vals[tid]) for tid in active_et_vals.index
                    }
                    errors.append(
                        f"element_types: all cells in the mesh must share the same"
                        f" ELEMENT_TYPE value, but found multiple: {sorted(unique_et_vals)}."
                        f" TYPE_ID → ELEMENT_TYPE mapping: {detail}."
                    )

            # Build column lookup: type_id → {ELEMENT_TYPE, OPTION2, OPTION3}
            et_lut: dict[int, dict[str, int]] = (
                et.set_index("TYPE_ID")[["ELEMENT_TYPE", "OPTION2", "OPTION3"]]
                .astype(int)
                .to_dict(orient="index")
            )

            # Sort cells by type_id for slice-based groupby – O(n log n) total
            et_sort = np.argsort(type_ids, kind="stable")
            s_tid = type_ids[et_sort]
            s_is2d = is_2d[et_sort]
            s_is3d = is_3d[et_sort]
            s_isq = is_quad[et_sort]

            u_tids, et_g_start = np.unique(s_tid, return_index=True)
            et_g_end = np.append(et_g_start[1:], len(s_tid))

            for idx in range(len(u_tids)):
                tid = int(u_tids[idx])
                if tid not in et_id_set:
                    continue
                sl = slice(int(et_g_start[idx]), int(et_g_end[idx]))
                row = et_lut[tid]
                elem_type = row["ELEMENT_TYPE"]
                opt2 = row["OPTION2"]
                opt3 = row["OPTION3"]
                any_2d = bool(np.any(s_is2d[sl]))
                any_3d = bool(np.any(s_is3d[sl]))
                any_quad = bool(np.any(s_isq[sl]))
                any_lin = bool(np.any(~s_isq[sl]))

                if any_2d and elem_type != 0:
                    errors.append(
                        f"element_types TYPE_ID={tid}: 2D cells require ELEMENT_TYPE=0,"
                        f" got {elem_type}."
                    )
                if any_3d and elem_type not in (1, 2):
                    errors.append(
                        f"element_types TYPE_ID={tid}: 3D cells require ELEMENT_TYPE 1 or 2,"
                        f" got {elem_type}."
                    )
                if any_3d and opt3 != 0:
                    errors.append(
                        f"element_types TYPE_ID={tid}: OPTION3 must be 0 for 3D cells,"
                        f" got {opt3}."
                    )
                if elem_type == 2 and any_lin:
                    errors.append(
                        f"element_types TYPE_ID={tid}: ELEMENT_TYPE=2 requires all cells"
                        " to be quadratic; linear VTK cells found."
                    )
                if elem_type in (0, 1):
                    if any_quad and any_lin:
                        errors.append(
                            f"element_types TYPE_ID={tid}: cells mix linear and quadratic"
                            " VTK types (OPTION2 must be uniform)."
                        )
                    elif any_quad and opt2 != 2:
                        errors.append(
                            f"element_types TYPE_ID={tid}: quadratic cells require OPTION2=2,"
                            f" got {opt2}."
                        )
                    elif any_lin and opt2 != 1:
                        errors.append(
                            f"element_types TYPE_ID={tid}: linear cells require OPTION2=1,"
                            f" got {opt2}."
                        )

    # ========================================================================
    # 2. dof_couplings DataFrame
    # ========================================================================
    if not dof.empty:
        req_dof = list(DOF_COUPLINGS_COLUMNS_TYPES.keys())
        miss_dof = set(req_dof) - set(dof.columns)
        if miss_dof:
            errors.append(f"dof_couplings missing columns: {sorted(miss_dof)}.")
        else:
            _valid_vars = {"A", "Az", "VOLT", "EMF"}
            mask_vars = ~dof["VARIABLE"].isin(_valid_vars)
            if mask_vars.any():
                bad_pairs_var = list(zip(
                    dof.loc[mask_vars, "CP_ID"].tolist(),
                    dof.loc[mask_vars, "VARIABLE"].tolist(),
                ))
                errors.append(
                    f"dof_couplings: VARIABLE must be one of {sorted(_valid_vars)};"
                    f" bad (CP_ID, value) pairs: {bad_pairs_var}."
                )
            bad_cp = dof.loc[~(dof["CP_ID"] > 0), "CP_ID"].tolist()
            if bad_cp:
                errors.append(
                    f"dof_couplings: CP_ID must be integer > 0; bad values: {bad_cp}."
                )

    # ========================================================================
    # 3. materials DataFrame
    # ========================================================================
    if not mat.empty:
        req_mat = list(MATERIALS_COLUMNS_TYPES.keys())
        miss_mat = set(req_mat) - set(mat.columns)
        if miss_mat:
            errors.append(f"materials missing columns: {sorted(miss_mat)}.")
        else:
            bad_mid = mat.loc[~(mat["MATERIAL_ID"] > 0), "MATERIAL_ID"].tolist()
            if bad_mid:
                errors.append(
                    f"materials: MATERIAL_ID must be > 0; bad values: {bad_mid}."
                )

            dups_m = (
                mat.loc[mat["MATERIAL_ID"].duplicated(keep=False), "MATERIAL_ID"]
                .unique()
                .tolist()
            )
            if dups_m:
                errors.append(f"materials: duplicate MATERIAL_ID values: {dups_m}.")

            # All non-ID columns must be positive real numbers
            for col in req_mat:
                if col == "MATERIAL_ID" or col not in mat.columns:
                    continue
                col_f = mat[col].astype(float)
                bad_mask = ~(col_f > 0.0)
                if bad_mask.any():
                    bad_pairs = list(
                        zip(
                            mat.loc[bad_mask, "MATERIAL_ID"].astype(int).tolist(),
                            mat.loc[bad_mask, col].tolist(),
                        )
                    )
                    preview = bad_pairs[:5]
                    suffix = "..." if len(bad_pairs) > 5 else ""
                    errors.append(
                        f"materials: '{col}' must be a positive real number;"
                        f" bad (MATERIAL_ID, value) pairs: {preview}{suffix}."
                    )

            if mat_ids is not None:
                mat_id_set: set[int] = set(mat["MATERIAL_ID"].astype(int))
                mesh_mat_set: set[int] = set(np.unique(mat_ids).tolist())
                miss_m = mesh_mat_set - mat_id_set
                if miss_m:
                    errors.append(
                        f"materials: cell MATERIAL values absent from MATERIAL_ID:"
                        f" {sorted(miss_m)}."
                    )
                unused_m = mat_id_set - mesh_mat_set
                if unused_m:
                    _log.warning(
                        "materials: MATERIAL_IDs not used by any cell: %s.", sorted(unused_m)
                    )

    # ========================================================================
    # 4. constants DataFrame
    # ========================================================================
    if not const.empty:
        req_const = list(CONSTANTS_COLUMNS_TYPES.keys())
        miss_const = set(req_const) - set(const.columns)
        if miss_const:
            errors.append(f"constants missing columns: {sorted(miss_const)}.")
        else:
            bad_cid = const.loc[~(const["CONSTANTS_ID"] > 0), "CONSTANTS_ID"].tolist()
            if bad_cid:
                errors.append(
                    f"constants: CONSTANTS_ID must be > 0; bad values: {bad_cid}."
                )

            dups_c = (
                const.loc[const["CONSTANTS_ID"].duplicated(keep=False), "CONSTANTS_ID"]
                .unique()
                .tolist()
            )
            if dups_c:
                errors.append(f"constants: duplicate CONSTANTS_ID values: {dups_c}.")

            # All float columns must contain finite values
            float_cols_c = [c for c in req_const if c != "CONSTANTS_ID"]
            for col in float_cols_c:
                if col not in const.columns:
                    continue
                col_f = const[col].astype(float)
                bad_mask_f = ~np.isfinite(col_f)
                if bad_mask_f.any():
                    bad_pairs_f = list(zip(
                        const.loc[bad_mask_f, "CONSTANTS_ID"].astype(int).tolist(),
                        const.loc[bad_mask_f, col].tolist(),
                    ))
                    preview = bad_pairs_f[:5]
                    suffix = "..." if len(bad_pairs_f) > 5 else ""
                    errors.append(
                        f"constants: '{col}' contains non-finite values;"
                        f" bad (CONSTANTS_ID, value) pairs: {preview}{suffix}."
                    )

            # ── cross-validate constants against mesh ────────────────────────
            if (
                const_ids is not None
                and type_ids is not None
                and et_ok
                and not et.empty
            ):
                const_id_set: set[int] = set(const["CONSTANTS_ID"].astype(int))
                mesh_const_set: set[int] = set(np.unique(const_ids).tolist())

                miss_c = mesh_const_set - const_id_set
                if miss_c:
                    errors.append(
                        f"constants: cell CONSTANT values absent from CONSTANTS_ID:"
                        f" {sorted(miss_c)}."
                    )

                unused_c = const_id_set - mesh_const_set
                if unused_c:
                    _log.warning(
                        "constants: CONSTANTS_IDs not used by any cell: %s.",
                        sorted(unused_c),
                    )

                # Build cell-level attribute arrays from element_types lookup.
                # Cells whose TYPE is not in element_types keep sentinel value -1
                # and are excluded from context-specific checks.
                et_dict: dict[int, tuple[int, int, int]] = {
                    int(r["TYPE_ID"]): (
                        int(r["ELEMENT_TYPE"]),
                        int(r["OPTION1"]),
                        int(r["OPTION3"]),
                    )
                    for _, r in et.iterrows()
                }
                cell_et = np.full(n_cells, -1, dtype=np.int8)
                cell_opt1 = np.full(n_cells, -1, dtype=np.int8)
                cell_opt3 = np.full(n_cells, -1, dtype=np.int8)
                for tid_i, (etype_i, opt1_i, opt3_i) in et_dict.items():
                    m = type_ids == tid_i
                    cell_et[m] = etype_i
                    cell_opt1[m] = opt1_i
                    cell_opt3[m] = opt3_i

                # Index constants table for O(1) row access
                const_idx_df = const.set_index("CONSTANTS_ID")

                # Sort cells by const_id for slice-based groupby – O(n log n) total
                c_sort = np.argsort(const_ids, kind="stable")
                s_cid = const_ids[c_sort]
                s_is2d_c = is_2d[c_sort]
                s_is3d_c = is_3d[c_sort]
                s_et_c = cell_et[c_sort]
                s_opt1_c = cell_opt1[c_sort]
                s_opt3_c = cell_opt3[c_sort]

                u_cids, c_g_start = np.unique(s_cid, return_index=True)
                c_g_end = np.append(c_g_start[1:], len(s_cid))

                for ci in range(len(u_cids)):
                    cid = int(u_cids[ci])
                    if cid not in const_id_set:
                        continue
                    sl = slice(int(c_g_start[ci]), int(c_g_end[ci]))
                    g_is2d = s_is2d_c[sl]
                    g_is3d = s_is3d_c[sl]
                    g_et = s_et_c[sl]
                    g_opt1 = s_opt1_c[sl]
                    g_opt3 = s_opt3_c[sl]

                    crow = const_idx_df.loc[cid]
                    is_2d_t0 = g_is2d & (g_et == 0)
                    is_3d_t12 = g_is3d & ((g_et == 1) | (g_et == 2))
                    any_2d_t0 = bool(np.any(is_2d_t0))
                    any_3d_t12 = bool(np.any(is_3d_t12))

                    # THK – only relevant for 2D type-0 cells
                    if any_2d_t0:
                        thk = float(crow["THK"])
                        if not (thk > 0.0):
                            errors.append(
                                f"constants CONSTANTS_ID={cid}: THK must be > 0"
                                f" for 2D type-0 cells; got {thk}."
                            )
                        elif bool(np.any(is_2d_t0 & (g_opt3 == 1))) and thk > 1.0:
                            errors.append(
                                f"constants CONSTANTS_ID={cid}: THK must be <= 1"
                                f" for OPTION3=1 cells; got {thk}."
                            )

                    # CS / TURNS / RAD / T_Z / R / S  (2D type-0, OPTION1=6)
                    if any_2d_t0 and bool(np.any(is_2d_t0 & (g_opt1 == 6))):
                        cs = float(crow["CS"])
                        turns = float(crow["TURNS"])
                        rad = float(crow["RAD"])
                        t_z = float(crow["T_Z"])
                        r = float(crow["R"])
                        s = float(crow["S"])
                        pfx = f"constants CONSTANTS_ID={cid} (2D OPTION1=6):"
                        if not (cs > 0.0):
                            errors.append(f"{pfx} CS must be > 0; got {cs}.")
                        if not (turns > 0.0):
                            errors.append(f"{pfx} TURNS must be > 0; got {turns}.")
                        if not (rad > 0.0):
                            errors.append(f"{pfx} RAD must be > 0; got {rad}.")
                        if t_z not in (-1.0, 1.0):
                            errors.append(f"{pfx} T_Z must be -1 or 1; got {t_z}.")
                        if not (r > 0.0):
                            errors.append(f"{pfx} R must be > 0; got {r}.")
                        if not (s >= 1.0):
                            errors.append(f"{pfx} S must be >= 1; got {s}.")

                    # CS / TURNS / VOL / T_X / T_Y / T_Z / R / S (3D type-1/2, OPTION1=6)
                    if any_3d_t12 and bool(np.any(is_3d_t12 & (g_opt1 == 6))):
                        cs = float(crow["CS"])
                        turns = float(crow["TURNS"])
                        vol = float(crow["VOL"])
                        t_x = float(crow["T_X"])
                        t_y = float(crow["T_Y"])
                        t_z = float(crow["T_Z"])
                        r = float(crow["R"])
                        s = float(crow["S"])
                        pfx = f"constants CONSTANTS_ID={cid} (3D OPTION1=6):"
                        if not (cs > 0.0):
                            errors.append(f"{pfx} CS must be > 0; got {cs}.")
                        if not (turns > 0.0):
                            errors.append(f"{pfx} TURNS must be > 0; got {turns}.")
                        if not (vol > 0.0):
                            errors.append(f"{pfx} VOL must be > 0; got {vol}.")
                        if not (r > 0.0):
                            errors.append(f"{pfx} R must be > 0; got {r}.")
                        if not (s >= 1.0):
                            errors.append(f"{pfx} S must be >= 1; got {s}.")
                        tx_ok = -1.0 <= t_x <= 1.0
                        ty_ok = -1.0 <= t_y <= 1.0
                        tz_ok = -1.0 <= t_z <= 1.0
                        if not (tx_ok and ty_ok and tz_ok):
                            errors.append(
                                f"{pfx} T_X/T_Y/T_Z must each be in [-1, 1];"
                                f" got ({t_x}, {t_y}, {t_z})."
                            )
                        else:
                            norm = (t_x ** 2 + t_y ** 2 + t_z ** 2) ** 0.5
                            if abs(norm - 1.0) > 1e-4:
                                errors.append(
                                    f"{pfx} (T_X, T_Y, T_Z) must form a unit vector"
                                    f" (norm={norm:.6f}); got ({t_x}, {t_y}, {t_z})."
                                )

    if errors:
        raise ValueError(
            "APEC model validation failed with {} error(s):\n  - {}".format(
                len(errors), "\n  - ".join(errors)
            )
        )

    def _format_df_for_log(df: pd.DataFrame) -> str:
        if df.empty:
            return "<empty>"
        return df.to_string(index=False)

    _log.info("Validated element types:\n%s", _format_df_for_log(et))
    _log.info("Validated materials:\n%s", _format_df_for_log(mat))
    _log.info("Validated real constants:\n%s", _format_df_for_log(const))
