"""The pyTEC3D model: a mesh and its model data."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import vtkStringArray
from vtkmodules.vtkCommonDataModel import vtkFieldData, vtkPointLocator, vtkUnstructuredGrid

from . import element_library, model_edit, selection, tree_gauge
from . import mesh as meshes
from ._logging import get_logger
from ._typing import PathType
from .coordinate_system import CoordinateSystem, coordinate_system_from_dict, coordinate_system_to_dict
from .errors import ModelLockedError, SolverStateError

LOGGER = get_logger("pytec3d.model")

#: Prefix of the field-data arrays written by the model.
FIELD_PREFIX = "pyTEC3D_"
#: Field-data array holding the model data as JSON text.
DATA_ARRAY_NAME = f"{FIELD_PREFIX}MODEL_DATA"
#: Field-data array holding the model name.
NAME_ARRAY_NAME = f"{FIELD_PREFIX}META_MESH_NAME"

#: Sections of the model data, each an object keyed by id.
DATA_SECTIONS = ("element_types", "materials", "constants", "dof_couplings", "nonlinear_laws", "coordinate_systems")


def empty_model_data() -> dict[str, Any]:
    """Model data with every section empty."""
    return {section: {} for section in DATA_SECTIONS}


_LOCKED_DATA_MESSAGE = (
    "The model data cannot change while the model is loaded into a solver session: the session uses the copy made "
    "at the load. Close the session (Solver.close_session()), change the data and load the model again."
)


def _refuse(*_: Any, **__: Any) -> Any:
    raise ModelLockedError(_LOCKED_DATA_MESSAGE)


class _FrozenDict(dict):  # type: ignore[type-arg]
    """A dict of the model data while the model is loaded: reading works as usual, any change raises."""

    __setitem__ = __delitem__ = __ior__ = _refuse
    clear = pop = popitem = setdefault = update = _refuse

    def __copy__(self) -> dict[str, Any]:
        return _thaw(self)

    def __deepcopy__(self, memo: dict[int, Any]) -> dict[str, Any]:
        return _thaw(self)

    def __reduce__(self) -> Any:
        return (dict, (_thaw(self),))


class _FrozenList(list):  # type: ignore[type-arg]
    """A list of the model data while the model is loaded: reading works as usual, any change raises."""

    __setitem__ = __delitem__ = __iadd__ = __imul__ = _refuse
    append = clear = extend = insert = pop = remove = reverse = sort = _refuse

    def __copy__(self) -> list[Any]:
        return _thaw(self)

    def __deepcopy__(self, memo: dict[int, Any]) -> list[Any]:
        return _thaw(self)

    def __reduce__(self) -> Any:
        return (list, (_thaw(self),))


def _freeze(value: Any) -> Any:
    """A read-only deep copy of model data (dicts and lists nested at any depth)."""
    if isinstance(value, Mapping):
        return _FrozenDict((key, _freeze(item)) for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return _FrozenList(_freeze(item) for item in value)
    if isinstance(value, np.ndarray):
        return _FrozenList(_freeze(item) for item in value.tolist())
    return value


def _thaw(value: Any) -> Any:
    """A plain, editable deep copy of (possibly frozen) model data."""
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_thaw(item) for item in value]
    return value


def _to_jsonable(value: Any, path: str = "data") -> Any:
    """Plain JSON types of dicts, sequences and NumPy values; keys become strings.

    Raises:
        ValueError: If two keys become the same string (``1`` and ``"1"``).
        TypeError: If a value is not a number, string, None, sequence or mapping.
    """
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for key, item in value.items():
            text = str(key)
            if text in out:
                raise ValueError(f"{path}: the key {text!r} is given more than once (for example as 1 and '1').")
            out[text] = _to_jsonable(item, f"{path}[{text}]")
        return out
    if isinstance(value, np.ndarray):
        return _to_jsonable(value.tolist(), path)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(item, path) for item in value]
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    raise TypeError(f"{path}: a {type(value).__name__} cannot be stored in the model data; use numbers, strings, "
                    "lists and dicts.")


class Model:
    """A pyTEC3D model: the mesh (a VTK unstructured grid) and the model data.

    The model data, :attr:`data`, is a JSON-like dict with these sections (:data:`DATA_SECTIONS`):

    * ``element_types``: ``{"<TYPE id>": {"element_type": 0, "option1": 6, "option2": 1, "option3": 0}}``;
    * ``materials``: ``{"<MATERIAL id>": {"mur": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "rho": [[...]], "M": [0, 0, 0]}}``;
    * ``constants``: ``{"<CONSTANT id>": {"thickness": 1.0, "T": [0, 0, 1], ...}}``;
    * ``dof_couplings``: ``{"<CP id>": ["A", "VOLT"]}``;
    * ``nonlinear_laws``: ``{"<id>": {"type": "bh_curve", "B": [0.0, ...], "H": [0.0, ...]}}`` (or a list of such
      entries): the laws of ``nonlinear_laws.json`` selected on the cells of the id, a MATERIAL, CONSTANT or TYPE id
      as the source of the law says, with the data of the law (here its curve);
    * ``coordinate_systems``: ``{"<COORDINATE_SYSTEM id>": {"type": "cylindrical", "origin": [...], "matrix": [...]}}``.

    The names in a materials or constants row are the ``mat.`` and ``rc.`` parameters of ``element_types.json``; a
    property no kernel of the element's domain needs can be left out. The cells select their rows through the cell
    arrays TYPE, MATERIAL, CONSTANT and COORDINATE_SYSTEM.

    The model data is stored in the mesh (VTK field data, as JSON text) when the model is saved, and handed to the
    solver with the views of the mesh arrays when it is loaded (:mod:`pytec3d.solver_input`).

    While the model is loaded into a solver session it is :attr:`locked`: the mesh and the model data cannot change,
    only the values of the mesh arrays (loads, boundary conditions). Close the session
    (:meth:`pytec3d.Solver.close_session`) to change them, then load the model again.
    """

    def __init__(self, mesh: vtkUnstructuredGrid | None = None, mesh_name: str | None = None) -> None:
        """Create a model, empty or around a mesh (whose arrays are not touched; see :meth:`from_grid`).

        Args:
            mesh: The mesh.
            mesh_name: Name of the model, used in output file names; ``"model"`` by default.

        Raises:
            TypeError: If ``mesh`` is not a ``vtkUnstructuredGrid`` or ``mesh_name`` not a str.
        """
        if mesh_name is not None and not isinstance(mesh_name, str):
            raise TypeError(f"mesh_name must be a str, got {type(mesh_name).__name__}.")
        self._mesh: vtkUnstructuredGrid | None = None
        self.mesh_name: str = mesh_name or "model"
        LOGGER.debug("Model '%s' created.", self.mesh_name)
        self._data: dict[str, Any] = empty_model_data()
        #: Absolute path of the JSON file the model data was read from (:meth:`read_json`); None once it is replaced.
        self.data_file: str | None = None
        self._locked = False
        self._frozen_data: dict[str, Any] | None = None
        if mesh is not None:
            self.mesh = mesh

    def __repr__(self) -> str:
        size = (
            f"{self._mesh.GetNumberOfPoints()} points, {self._mesh.GetNumberOfCells()} cells"
            if self._mesh is not None else "no mesh"
        )
        return f"Model({self.mesh_name!r}, {size}{', locked' if self._locked else ''})"

    # ------------------------------------------------------------------------------------------- creation

    @classmethod
    def from_grid(cls, mesh: vtkUnstructuredGrid, mesh_name: str | None = None) -> Model:
        """Create a model around a mesh, with the model data and the name stored in it, if any.

        The arrays the model needs are added to the mesh (:meth:`update_grid`).

        Args:
            mesh: The mesh (kept, not copied).
            mesh_name: Name of the model; the one stored in the mesh, or ``"model"``, when None.
        """
        model = cls(mesh=mesh)
        model._read_field_data()
        if mesh_name is not None:
            model.mesh_name = mesh_name
        model.update_grid()
        return model

    @classmethod
    def read_model(cls, path: PathType, reset_arrays: bool = False) -> Model:
        """Read a model from a ``.vtu`` / ``.pvtu`` file: the mesh and the model data stored in it.

        The model is named after the file.

        Args:
            path: File to read.
            reset_arrays: Set the arrays of the element types back to their default values, clear the selection and
                compute the midside-node mask.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file is not a valid mesh.
        """
        file_path = os.path.abspath(os.fspath(path))
        mesh = meshes.read_mesh(file_path, reset_arrays=reset_arrays)
        if reset_arrays:
            selection.ensure_selection_arrays(mesh)
            tree_gauge.add_midside_nodes_mask(mesh)
        return cls.from_grid(mesh, mesh_name=os.path.splitext(os.path.basename(file_path))[0])

    def save_model(self, path: PathType, binary: bool = True) -> None:
        """Write the mesh, with the name and the model data in its field data, to a ``.vtu`` file.

        Raises:
            SolverStateError: If the model has no mesh.
            ValueError, TypeError: If the model data cannot be written as JSON (see :meth:`data_json`).
            OSError: If the file cannot be written.
        """
        mesh = self.mesh
        field_data = mesh.GetFieldData()
        self._remove_field_arrays(field_data)
        self._add_string(field_data, NAME_ARRAY_NAME, self.mesh_name)
        self._add_string(field_data, DATA_ARRAY_NAME, self.data_json())
        meshes.save_mesh(mesh, path, binary=binary)

    # ------------------------------------------------------------------------------------------- lock

    @property
    def locked(self) -> bool:
        """True while the model is loaded into a solver session."""
        return self._locked

    def lock(self) -> None:
        """Fix the mesh and the model data; called by the solver when it loads the model.

        From now on :attr:`data` is a read-only copy: an edit raises :class:`~pytec3d.errors.ModelLockedError`.
        """
        self._frozen_data = _freeze(self._data)
        self._locked = True

    def unlock(self) -> None:
        """Allow changes to the mesh and the model data again; called by the solver when it closes the session."""
        self._locked = False
        self._frozen_data = None

    def _require_unlocked(self, what: str) -> None:
        if self._locked:
            raise ModelLockedError(
                f"Cannot change {what} while the model is loaded into a solver session: only array values can "
                "change. Close the session (Solver.close_session()) and load the model again."
            )

    # ------------------------------------------------------------------------------------------- mesh

    @property
    def has_mesh(self) -> bool:
        """True when the model has a mesh."""
        return self._mesh is not None

    @property
    def mesh(self) -> vtkUnstructuredGrid:
        """The mesh.

        Raises:
            SolverStateError: When reading it before a mesh is set.
            TypeError: When setting something else than a ``vtkUnstructuredGrid``.
            ModelLockedError: When setting it while the model is locked.
        """
        if self._mesh is None:
            raise SolverStateError("The model has no mesh: read one first (Solver.read_model).")
        return self._mesh

    @mesh.setter
    def mesh(self, value: vtkUnstructuredGrid) -> None:
        meshes.check_grid(value)
        self._require_unlocked("the mesh")
        self._mesh = value

    def update_grid(self, reset: bool = False) -> None:
        """Add to the mesh the arrays the model needs (:func:`pytec3d.mesh.update_grid`).

        Those of the element types of the model data, or of every element type of the library while it names none.

        Args:
            reset: Also set the arrays of the element types back to their default values.
        """
        meshes.update_grid(self.mesh, self.element_type_names() or None, reset=reset)

    # ------------------------------------------------------------------------------------------- data

    @property
    def data(self) -> dict[str, Any]:
        """The model data: edit it in place, or assign a new dict (missing sections are added, empty).

        The solver copies it when the model is loaded. While the model is loaded the data is read-only: reading works
        as usual (the values are ``dict`` and ``list`` objects), any edit raises.

        Raises:
            ModelLockedError: When assigning it, or editing it in place, while the model is locked.
            TypeError: When assigning something else than a dict of dicts.
        """
        if self._locked and self._frozen_data is not None:
            return self._frozen_data
        return self._data

    @data.setter
    def data(self, value: Mapping[str, Any]) -> None:
        self._require_unlocked("the model data")
        if not isinstance(value, Mapping):
            raise TypeError(f"data must be a dict, got {type(value).__name__}.")
        data = empty_model_data()
        data.update(value)
        for section in DATA_SECTIONS:
            if not isinstance(data[section], Mapping):
                raise TypeError(f"data['{section}'] must be a dict, got {type(data[section]).__name__}.")
        self._data = data
        self.data_file = None

    def data_json(self) -> str:
        """The model data as compact JSON text: what the solver and the saved mesh receive.

        Raises:
            ValueError: If a value is not finite, or a key is given twice (as 1 and "1").
            TypeError: If a value is not a number, string, list or dict.
        """
        return json.dumps(_to_jsonable(self._data), separators=(",", ":"), allow_nan=False)

    def read_json(self, path: PathType) -> None:
        """Replace the model data by the content of a JSON file, and add the arrays its element types need.

        Raises:
            ModelLockedError: If the model is locked.
            ValueError: If the file does not hold a JSON object.
            OSError: If the file cannot be read.
        """
        self._require_unlocked("the model data")
        file_path = os.fspath(path)
        LOGGER.info("Reading the model data from '%s'.", file_path)
        with open(file_path, encoding="utf-8") as handle:
            loaded = json.load(handle)
        if not isinstance(loaded, dict):
            raise ValueError(f"'{file_path}' must hold a JSON object, got {type(loaded).__name__}.")
        unknown = sorted(set(loaded) - set(DATA_SECTIONS))
        if unknown:
            LOGGER.warning("Model data sections not used by the solver: %s.", unknown)
        self.data = loaded
        self.data_file = os.path.abspath(file_path)
        laws = sum(len(v) if isinstance(v, list) else 1 for v in self._data["nonlinear_laws"].values())
        LOGGER.info(
            "Model data: element types %s; %d material(s), %d constants row(s), %d DOF coupling(s), %d nonlinear "
            "law(s), %d coordinate system(s).",
            sorted(set(self.element_type_names())) or "none", len(self._data["materials"]),
            len(self._data["constants"]), len(self._data["dof_couplings"]), laws,
            len(self._data["coordinate_systems"]),
        )
        if self._mesh is not None:
            self.update_grid()

    def write_json(self, path: PathType) -> None:
        """Write the model data to a JSON file (indented).

        Raises:
            ValueError, TypeError: As :meth:`data_json`.
        """
        file_path = os.fspath(path)
        with open(file_path, "w", encoding="utf-8") as handle:
            json.dump(_to_jsonable(self._data), handle, indent=2, allow_nan=False)
        LOGGER.info("Model data written to '%s'.", file_path)

    def element_type_names(self) -> list[str]:
        """Names of the element types of the ``element_types`` rows."""
        return [
            element_library.element_type_name(row["element_type"])
            for row in self._data.get("element_types", {}).values()
            if isinstance(row, Mapping) and "element_type" in row
        ]

    def _row(self, section: str, row_id: int) -> Any:
        """``data[section]`` row ``row_id`` (keyed as int or str), or None."""
        rows = self._data[section]
        for key in (row_id, str(row_id)):
            if key in rows:
                return rows[key]
        return None

    def real_constant(self, constants_id: int, name: str, default: Any = None) -> Any:
        """The real constant ``name`` of the constants row ``constants_id``, or ``default`` when undefined."""
        row = self._row("constants", int(constants_id))
        return default if row is None else row.get(name, default)

    def coordinate_system(self, cs_id: int) -> CoordinateSystem:
        """The coordinate system ``cs_id`` of ``data["coordinate_systems"]``, as an object.

        Used for the rotations of the elements of a COORDINATE_SYSTEM (:mod:`pytec3d.solver_input`).

        Raises:
            KeyError: If the model data does not define it.
            ValueError: If its attributes are invalid.
        """
        spec = self._row("coordinate_systems", int(cs_id))
        if spec is None:
            raise KeyError(f"data['coordinate_systems'] has no entry for the id {int(cs_id)}.")
        try:
            return coordinate_system_from_dict(spec)
        except ValueError as exc:
            raise ValueError(f"coordinate_systems[{int(cs_id)}]: {exc}") from exc

    def set_coordinate_system(self, cs_id: int, coordinate_system: CoordinateSystem) -> None:
        """Store a coordinate system as ``data["coordinate_systems"][cs_id]`` (type, origin, matrix).

        Raises:
            ModelLockedError: If the model is locked.
            ValueError: If ``cs_id`` is not > 0 (0 is the global system).
        """
        self._require_unlocked("the model data")
        if int(cs_id) <= 0:
            raise ValueError(f"cs_id must be > 0 (0 is the global system), got {cs_id}.")
        self._data["coordinate_systems"][str(int(cs_id))] = coordinate_system_to_dict(coordinate_system)
        LOGGER.debug("Coordinate system %d: %r.", int(cs_id), coordinate_system)

    def add_dof_couplings(self, couplings: Mapping[int, Sequence[str]] | None) -> None:
        """Add coupled DOF fields to CP ids, keeping those already there.

        Args:
            couplings: ``{cp_id: [field, ...]}``.

        Raises:
            ModelLockedError: If the model is locked.
        """
        self._require_unlocked("the model data")
        for cp_id, fields in (couplings or {}).items():
            existing = self._row("dof_couplings", int(cp_id))
            merged = list(existing) if existing is not None else []
            merged.extend(str(f) for f in fields if str(f) not in merged)
            self._data["dof_couplings"].pop(int(cp_id), None)
            self._data["dof_couplings"][str(int(cp_id))] = merged

    def periodic_boundaries(
        self,
        coordinate_system: CoordinateSystem,
        variables: Sequence[str] | str,
        dx: float = 0.0,
        dy: float = 0.0,
        dz: float = 0.0,
        array_name: str = "CP_0",
        tolerance: float = 1e-6,
    ) -> dict[int, list[str]]:
        """Couple the selected points with their periodic images.

        Each selected point is moved by ``(dx, dy, dz)`` in the local coordinates of ``coordinate_system`` (e.g.
        ``dy=30`` turns by 30 degrees in a cylindrical system) and matched with the nearest mesh point. Both points of a
        pair get the same new CP id in ``array_name``, after the ids already used by any CP_k array, and the new ids
        couple ``variables`` in ``data["dof_couplings"]``. Points that already have an id in ``array_name``, or whose
        image is farther than ``tolerance``, are skipped.

        Args:
            coordinate_system: System of the shift.
            variables: The DOF fields the pairs share, e.g. ``("A", "VOLT")``.
            dx: Shift of the first local coordinate.
            dy: Shift of the second local coordinate (degrees for an angle).
            dz: Shift of the third local coordinate.
            array_name: CP point array that receives the ids.
            tolerance: Largest distance between an image and the point matched with it.

        Returns:
            ``{cp_id: variables}`` of the new pairs (also added to the model data).

        Raises:
            ModelLockedError: If the model is locked.
            MissingArrayError: If ``array_name`` is missing.
            ValueError: If ``tolerance`` is not positive.
        """
        self._require_unlocked("the model data")
        if tolerance <= 0:
            raise ValueError(f"tolerance must be positive, got {tolerance}.")
        fields = [variables] if isinstance(variables, str) else [str(v) for v in variables]
        mesh = self.mesh

        local, sources = model_edit.to_local_coordinates(mesh, coordinate_system, "node", selected_only=True)
        if local.size == 0:
            LOGGER.warning("periodic_boundaries: no point is selected; nothing coupled.")
            return {}
        shifted = np.asarray(local, dtype=float).copy()
        shifted[:, 0] += float(dx)
        shifted[:, 1] += float(dy)
        shifted[:, 2] += float(dz)
        images, _ = model_edit.to_global_coordinates(mesh, shifted, coordinate_system, "node", selected_only=True)

        coordinates = meshes.points(mesh)
        locator = vtkPointLocator()
        locator.SetDataSet(mesh)
        locator.BuildLocator()
        cp_array = meshes.get_array(mesh, array_name, "node")
        cp_values = vtknp.vtk_to_numpy(cp_array)
        used = self._used_coupling_ids()
        next_id = 1
        while next_id in used:
            next_id += 1

        new_ids: list[int] = []
        paired: set[int] = set()
        unmatched = 0
        for source, image in zip(sources.tolist(), images, strict=True):
            if source in paired or cp_values[source] != 0:
                continue
            target = int(locator.FindClosestPoint(image))
            if target < 0 or target == source or target in paired or cp_values[target] != 0:
                continue
            if np.linalg.norm(coordinates[target] - image) > tolerance:
                unmatched += 1
                continue
            cp_values[source] = cp_values[target] = next_id
            paired.update((source, target))
            used.add(next_id)
            new_ids.append(next_id)
            while next_id in used:
                next_id += 1
        cp_array.Modified()

        couplings = {cp_id: list(fields) for cp_id in new_ids}
        self.add_dof_couplings(couplings)
        LOGGER.info(
            "Periodic boundaries: %d pair(s) of %d selected point(s) coupled in '%s' (ids %s), fields %s; "
            "%d image(s) farther than %g.",
            len(new_ids), sources.size, array_name,
            f"{new_ids[0]}..{new_ids[-1]}" if new_ids else "none", fields, unmatched, tolerance,
        )
        return couplings

    def _used_coupling_ids(self) -> set[int]:
        """CP ids already used by the contiguous family of point arrays CP_0, CP_1, ..."""
        used: set[int] = set()
        point_data = self.mesh.GetPointData()
        index = 0
        while (array := point_data.GetArray(f"CP_{index}")) is not None:
            values = vtknp.vtk_to_numpy(array)
            used.update(values[np.isfinite(values) & (values != 0)].astype(int).tolist())
            index += 1
        return used

    def set_nonlinear_law(self, row_id: int, law: str, **data: Any) -> None:
        """Select a nonlinear law on the cells of an id, with its data: ``data["nonlinear_laws"][row_id]``.

        The id is a MATERIAL, CONSTANT or TYPE id, as the source of the law says (``nonlinear_laws.json``). An entry of
        the same law for that id is replaced; entries of other laws are kept.

        Args:
            row_id: The MATERIAL, CONSTANT or TYPE id.
            law: Name of the law, e.g. ``"bh_curve"``.
            **data: The data of the law, e.g. its curve arrays ``B=[...], H=[...]``.

        Raises:
            ModelLockedError: If the model is locked.
        """
        self._require_unlocked("the model data")
        laws = self._data["nonlinear_laws"]
        key = str(int(row_id))
        existing = laws.pop(int(row_id), None) or laws.pop(key, None) or []
        entries = [e for e in (existing if isinstance(existing, list) else [existing]) if e.get("type") != law]
        entries.append({"type": law, **{name: _to_jsonable(value) for name, value in data.items()}})
        laws[key] = entries[0] if len(entries) == 1 else entries
        LOGGER.info("Nonlinear law '%s' set on id %s (data: %s).", law, key, ", ".join(sorted(data)) or "none")

    def set_bh_curve(self, material_id: int, table: Any) -> None:
        """Select the ``bh_curve`` law on a material, with its BH curve given as rows ``[H, B]``.

        Args:
            material_id: The MATERIAL id.
            table: Shape (N, 2): H and B, B increasing.

        Raises:
            ModelLockedError: If the model is locked.
            ValueError: If ``table`` does not have shape (N, 2).
        """
        rows = np.asarray(table, dtype=float)
        if rows.ndim != 2 or rows.shape[1] != 2:
            raise ValueError(f"table must have shape (N, 2) with columns [H, B], got {rows.shape}.")
        self.set_nonlinear_law(material_id, "bh_curve", B=rows[:, 1].tolist(), H=rows[:, 0].tolist())

    # ------------------------------------------------------------------------------------------- field data

    @staticmethod
    def _remove_field_arrays(field_data: vtkFieldData) -> None:
        names = [field_data.GetArrayName(i) for i in range(field_data.GetNumberOfArrays())]
        for name in names:
            if name and name.startswith(FIELD_PREFIX):
                field_data.RemoveArray(name)

    @staticmethod
    def _add_string(field_data: vtkFieldData, name: str, text: str) -> None:
        array = vtkStringArray()
        array.SetName(name)
        array.SetNumberOfValues(1)
        array.SetValue(0, text)
        field_data.AddArray(array)

    @staticmethod
    def _read_string(field_data: vtkFieldData, name: str) -> str | None:
        array = vtkStringArray.SafeDownCast(field_data.GetAbstractArray(name))
        return None if array is None or array.GetNumberOfValues() == 0 else array.GetValue(0)

    def _read_field_data(self) -> None:
        """Take the name and the model data stored in the field data of the mesh, if any."""
        field_data = self.mesh.GetFieldData()
        name = self._read_string(field_data, NAME_ARRAY_NAME)
        if name:
            self.mesh_name = name
        text = self._read_string(field_data, DATA_ARRAY_NAME)
        if text is None:
            self._data = empty_model_data()
        else:
            self.data = json.loads(text)
