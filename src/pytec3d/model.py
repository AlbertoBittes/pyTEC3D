import json
import os
from os import PathLike
from typing import Any, Dict, Mapping, Optional, Sequence, Union

import numpy as np
import vtk

from .coordinate_system import (
    CoordinateSystem,
    coordinate_system_from_dict,
    coordinate_system_to_dict,
)
from .logging_utils import get_logger
from .model_edit import ModelEdit


LOGGER = get_logger("pytec3d.pyTEC3D_Model")

FIELD_PREFIX = "pyTEC3D_"

#: Name of the VTK field-data array holding the model data as JSON text.
DATA_ARRAY_NAME = f"{FIELD_PREFIX}MODEL_DATA"

#: Sections of the model data. Every one is an object keyed by id (or a name for the coupling ids).
DATA_SECTIONS = (
    "element_types",
    "materials",
    "constants",
    "dof_couplings",
    "bh_curves",
    "coordinate_systems",
)


def empty_model_data() -> Dict[str, Any]:
    """Return model data with every section empty."""
    return {section: {} for section in DATA_SECTIONS}


def _to_jsonable(value: Any, path: str = "data") -> Any:
    """Convert dicts, sequences and NumPy values to plain JSON types; dict keys become strings."""
    if isinstance(value, Mapping):
        out: Dict[str, Any] = {}
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
    raise TypeError(f"{path}: {type(value).__name__} cannot be stored in the model data; use numbers, strings, lists and dicts.")


class pyTEC3D_Model:
    """pyTEC3D model: the mesh and the model data.

    The model data is a JSON-like dict, ``model.data``, with these sections (see ``DATA_SECTIONS``):

    - ``element_types``: ``{"<TYPE id>": {"element_type": 0, "option1": 6, "option2": 1, "option3": 0}}``
    - ``materials``: ``{"<MATERIAL id>": {"mur": [[1,0,0],[0,1,0],[0,0,1]], "rho": [[...]], "M": [0, 0, 0]}}``
    - ``constants``: ``{"<CONSTANT id>": {"thickness": 1.0, "T": [0, 0, 1], ...}}``
    - ``dof_couplings``: ``{"<CP id>": ["A", "VOLT"]}``
    - ``bh_curves``: ``{"<MATERIAL id>": {"H": [0.0, ...], "B": [0.0, ...]}}``
    - ``coordinate_systems``: ``{"<COORDINATE_SYSTEM id>": {"type": "cylindrical", "origin": [...], "matrix": [[...]]}}``

    The names inside a material or constants entry are the ``mat.`` and ``rc.`` parameters of
    ``element_types.json``; a property no kernel of the element's domain needs can be left out.
    The data is stored in the mesh (VTK field data) as JSON text when the model is saved, and
    handed to the solver as JSON text when it is loaded.
    """

    def __init__(
        self,
        mesh: Optional[vtk.vtkUnstructuredGrid] = None,
        mesh_name: Optional[str] = None,
    ) -> None:
        """Create a model, optionally around an existing mesh.

        Args:
            mesh: Optional unstructured grid.
            mesh_name: Optional logical name of the model, used in output file names.

        Raises:
            TypeError: If mesh or mesh_name has the wrong type.
        """
        if mesh is not None and not isinstance(mesh, vtk.vtkUnstructuredGrid):
            raise TypeError("mesh must be a vtk.vtkUnstructuredGrid or None.")
        if mesh_name is not None and not isinstance(mesh_name, str):
            raise TypeError("mesh_name must be a str or None.")

        self._mesh: Optional[vtk.vtkUnstructuredGrid] = None
        self._mesh_name: Optional[str] = mesh_name if mesh_name is not None else "model"
        self._data: Dict[str, Any] = empty_model_data()
        if mesh is not None:
            self.mesh = mesh

    @property
    def mesh(self) -> vtk.vtkUnstructuredGrid:
        """The mesh of the model.

        Raises:
            RuntimeError: If no mesh is loaded.
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        return self._mesh

    @mesh.setter
    def mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        if not isinstance(value, vtk.vtkUnstructuredGrid):
            raise TypeError("mesh must be a vtk.vtkUnstructuredGrid.")
        self._mesh = value

    @property
    def data(self) -> Dict[str, Any]:
        """The model data dict; edit it in place or assign a new one."""
        return self._data

    @data.setter
    def data(self, value: Mapping[str, Any]) -> None:
        if not isinstance(value, Mapping):
            raise TypeError("data must be a dict.")
        data = empty_model_data()
        data.update(value)
        for section in DATA_SECTIONS:
            if not isinstance(data[section], Mapping):
                raise TypeError(f"data['{section}'] must be a dict.")
        self._data = data

    # ------------------------------------------------------------------ JSON

    def data_json(self) -> str:
        """Return the model data as compact JSON text (what the solver and the VTK field data get).

        Raises:
            ValueError: If a value is not finite or a key is repeated.
            TypeError: If a value is not a number, string, list or dict.
        """
        return json.dumps(_to_jsonable(self._data), separators=(",", ":"), allow_nan=False)

    def read_json(self, json_file_path: Union[str, PathLike[str]]) -> None:
        """Replace the model data with the content of a JSON file.

        Args:
            json_file_path: Path of the model data JSON.

        Raises:
            ValueError: If the file does not contain a JSON object.
        """
        path = os.fspath(json_file_path)
        LOGGER.info("Reading model data from '%s'.", path)
        with open(path, encoding="utf-8") as handle:
            loaded = json.load(handle)
        if not isinstance(loaded, dict):
            raise ValueError(f"'{path}' must contain a JSON object.")
        unknown = sorted(set(loaded) - set(DATA_SECTIONS))
        if unknown:
            LOGGER.warning("Model data sections not used by the solver: %s.", unknown)
        self.data = loaded

    def write_json(self, json_file_path: Union[str, PathLike[str]]) -> None:
        """Write the model data to a JSON file.

        Args:
            json_file_path: Path of the file to write.
        """
        path = os.fspath(json_file_path)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(_to_jsonable(self._data), handle, indent=2, allow_nan=False)
        LOGGER.info("Model data written to '%s'.", path)

    # ------------------------------------------------- access to the sections

    def _entry(self, section: str, entry_id: int) -> Any:
        """Return ``data[section]`` for an id given as int or str, or None when absent."""
        entries = self._data[section]
        for key in (entry_id, str(entry_id)):
            if key in entries:
                return entries[key]
        return None

    def coordinate_system(self, cs_id: int) -> CoordinateSystem:
        """Build the coordinate system ``cs_id`` from its attributes in ``data['coordinate_systems']``.

        Args:
            cs_id: COORDINATE_SYSTEM id used by the mesh.

        Returns:
            CoordinateSystem: A new object for the stored attributes.

        Raises:
            KeyError: If the model data does not define it.
            ValueError: If its attributes are invalid.
        """
        spec = self._entry("coordinate_systems", int(cs_id))
        if spec is None:
            raise KeyError(f"data['coordinate_systems'] has no entry for id {int(cs_id)}.")
        try:
            return coordinate_system_from_dict(spec)
        except ValueError as exc:
            raise ValueError(f"coordinate_systems[{int(cs_id)}]: {exc}") from exc

    def set_coordinate_system(self, cs_id: int, coordinate_system: CoordinateSystem) -> None:
        """Store the attributes (type, origin, matrix) of a coordinate system object under ``cs_id``.

        Args:
            cs_id: COORDINATE_SYSTEM id used by the mesh (> 0).
            coordinate_system: Object to describe in the model data.
        """
        if int(cs_id) <= 0:
            raise ValueError("cs_id must be > 0.")
        self._data["coordinate_systems"][str(int(cs_id))] = coordinate_system_to_dict(coordinate_system)

    def add_dof_couplings(self, couplings: Mapping[int, Sequence[str]]) -> None:
        """Add coupled variables to the CP ids, keeping those already defined.

        Args:
            couplings: ``{cp_id: [variable, ...]}``, as returned by ``ModelEdit.periodic_boundaries``.
        """
        for cp_id, variables in couplings.items():
            existing = self._entry("dof_couplings", int(cp_id))
            merged = list(existing) if existing is not None else []
            merged.extend(str(v) for v in variables if str(v) not in merged)
            self._data["dof_couplings"].pop(int(cp_id), None)
            self._data["dof_couplings"][str(int(cp_id))] = merged

    def set_bh_curve(self, material_id: int, table: Any) -> None:
        """Store a BH curve given as rows of ``[H, B]`` under the material id.

        Args:
            material_id: MATERIAL id used by the mesh (> 0).
            table: Array-like of shape (N, 2) with columns H and B.
        """
        arr = np.asarray(table, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError(f"table must have shape (N, 2) with columns [H, B], got {arr.shape}.")
        self._data["bh_curves"][str(int(material_id))] = {"H": arr[:, 0].tolist(), "B": arr[:, 1].tolist()}

    def real_constant(self, constants_id: int, name: str, default: Any = None) -> Any:
        """Return the real constant ``name`` of the constants ``constants_id``, or ``default`` if undefined."""
        entry = self._entry("constants", int(constants_id))
        if entry is None:
            return default
        return entry.get(name, default)

    # ------------------------------------------------------ VTK serialization

    @staticmethod
    def _meta_array_name(key: str) -> str:
        return f"{FIELD_PREFIX}META_{key}"

    @staticmethod
    def _remove_prefixed_field_arrays(field_data: vtk.vtkFieldData) -> None:
        """Remove every array of the field data written by the model."""
        names = []
        for i in range(field_data.GetNumberOfArrays()):
            arr = field_data.GetAbstractArray(i)
            if arr is None:
                continue
            name = arr.GetName()
            if isinstance(name, str) and name.startswith(FIELD_PREFIX):
                names.append(name)
        for name in names:
            field_data.RemoveArray(name)

    @staticmethod
    def _add_string_array(field_data: vtk.vtkFieldData, name: str, text: str) -> None:
        arr = vtk.vtkStringArray()
        arr.SetName(name)
        arr.SetNumberOfValues(1)
        arr.SetValue(0, text)
        field_data.AddArray(arr)

    @staticmethod
    def _read_string(field_data: vtk.vtkFieldData, name: str) -> Optional[str]:
        arr = vtk.vtkStringArray.SafeDownCast(field_data.GetAbstractArray(name))
        if arr is None or arr.GetNumberOfValues() == 0:
            return None
        return arr.GetValue(0)

    def _write_model_field_data(self) -> None:
        """Store the mesh name and the model data (as JSON text) in the field data of the mesh."""
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        field_data = self._mesh.GetFieldData()
        self._remove_prefixed_field_arrays(field_data)
        self._add_string_array(field_data, self._meta_array_name("MESH_NAME"), self._mesh_name or "")
        self._add_string_array(field_data, DATA_ARRAY_NAME, self.data_json())

    def _read_model_field_data(self) -> None:
        """Load the mesh name and the model data from the field data of the mesh."""
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        field_data = self._mesh.GetFieldData()
        mesh_name = self._read_string(field_data, self._meta_array_name("MESH_NAME"))
        if mesh_name:
            self._mesh_name = mesh_name

        text = self._read_string(field_data, DATA_ARRAY_NAME)
        if text is None:
            self._data = empty_model_data()
            return
        self.data = json.loads(text)

    @classmethod
    def read_model(
        cls,
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = False,
    ) -> "pyTEC3D_Model":
        """Read a model (mesh and the model data stored in it) from a VTK file.

        Args:
            vtk_file_path: Path of the VTK file.
            reset_arrays: Whether to reinitialize the required VTK arrays on load.

        Returns:
            pyTEC3D_Model: The loaded model.
        """
        absolute_path = os.path.abspath(os.fspath(vtk_file_path))
        loaded_mesh = ModelEdit.read_mesh(absolute_path, reset_arrays=reset_arrays)
        resolved_name = os.path.splitext(os.path.basename(absolute_path))[0]
        return cls.from_grid(loaded_mesh, mesh_name=resolved_name)

    @classmethod
    def from_grid(
        cls,
        mesh: vtk.vtkUnstructuredGrid,
        mesh_name: Optional[str] = None,
    ) -> "pyTEC3D_Model":
        """Create a model around a grid, taking the model data stored in it if any.

        Args:
            mesh: Unstructured grid.
            mesh_name: Optional logical name of the model.

        Returns:
            pyTEC3D_Model: The model.
        """
        model = cls(mesh=mesh, mesh_name=mesh_name)
        model._read_model_field_data()
        return model

    def save_model(self, output_file: Union[str, PathLike[str]], binary: bool = True) -> None:
        """Save the mesh with the model data stored in its field data.

        Args:
            output_file: Path of the VTK file.
            binary: Whether to write binary data.

        Raises:
            RuntimeError: If no mesh is loaded.
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        self._write_model_field_data()
        ModelEdit.save_mesh(self.mesh, output_file, binary=binary)
