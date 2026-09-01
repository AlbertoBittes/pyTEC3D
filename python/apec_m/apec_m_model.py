import copy
import os
from os import PathLike
from typing import Any, Dict, Optional, Sequence, Union

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

from .constants import (
    CONSTANTS_COLUMNS_TYPES,
    DEFAULT_COLUMNS_ATTRIBUTES_TYPES,
    DEFAULT_SOLVER_SETTINGS,
    DEFAULT_SOLVER_SETTINGS_TYPES,
    DEFAULT_TIME_PARAMETERS,
    DEFAULT_TIME_PARAMETERS_TYPES,
    DOF_COUPLINGS_COLUMNS_TYPES,
    ELEMENT_TYPES_COLUMNS_TYPES,
    MATERIALS_COLUMNS_TYPES,
    SHEET_PROCESSING_MAP,
)
from .excel import Excel
from .logging_utils import get_logger
from .model_edit import ModelEdit


LOGGER = get_logger("apec_m.APEC_M_Model")

FIELD_PREFIX = "APEC_M_"

TABLE_SPECS = {
    "element_types": ELEMENT_TYPES_COLUMNS_TYPES,
    "materials": MATERIALS_COLUMNS_TYPES,
    "constants": CONSTANTS_COLUMNS_TYPES,
    "dof_couplings": DOF_COUPLINGS_COLUMNS_TYPES,
}


class APEC_M_Model:
    """APEC model container for mesh and model metadata."""

    def __init__(
        self,
        mesh: Optional[vtk.vtkUnstructuredGrid] = None,
        mesh_name: Optional[str] = None,
    ) -> None:
        if mesh is not None and not isinstance(mesh, vtk.vtkUnstructuredGrid):
            raise TypeError("mesh must be a vtk.vtkUnstructuredGrid or None.")
        if mesh_name is not None and not isinstance(mesh_name, str):
            raise TypeError("mesh_name must be a str or None.")

        self._mesh: Optional[vtk.vtkUnstructuredGrid] = mesh
        self._mesh_name: Optional[str] = mesh_name

        self.element_types: pd.DataFrame = pd.DataFrame(
            columns=list(ELEMENT_TYPES_COLUMNS_TYPES.keys())
        )
        self.materials: pd.DataFrame = pd.DataFrame(
            columns=list(MATERIALS_COLUMNS_TYPES.keys())
        )
        self.constants: pd.DataFrame = pd.DataFrame(
            columns=list(CONSTANTS_COLUMNS_TYPES.keys())
        )
        self.dof_couplings: pd.DataFrame = pd.DataFrame(
            columns=list(DOF_COUPLINGS_COLUMNS_TYPES.keys())
        )
        self.time_step: Dict[str, Any] = copy.deepcopy(DEFAULT_TIME_PARAMETERS)
        self.solver_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self.coordinate_systems: Dict[int, Any] = {}
        self.bh_curves: Dict[int, np.ndarray] = {}

    @property
    def mesh(self) -> vtk.vtkUnstructuredGrid:
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        return self._mesh

    @mesh.setter
    def mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        raise AttributeError("Direct mesh assignment is not allowed. Use read_model().")

    @staticmethod
    def _table_rows_array_name(table_name: str) -> str:
        return f"{FIELD_PREFIX}TABLE_{table_name.upper()}_ROWS"

    @staticmethod
    def _table_column_array_name(table_name: str, column_name: str) -> str:
        return f"{FIELD_PREFIX}TABLE_{table_name.upper()}_{column_name}"

    @staticmethod
    def _meta_array_name(key: str) -> str:
        return f"{FIELD_PREFIX}META_{key}"

    @staticmethod
    def _solver_setting_array_name(key: str) -> str:
        return f"{FIELD_PREFIX}SOLVER_SETTING_{key}"

    @staticmethod
    def _remove_prefixed_field_arrays(field_data: vtk.vtkFieldData) -> None:
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
    def _add_string_array(field_data: vtk.vtkFieldData, name: str, values: Sequence[str]) -> None:
        arr = vtk.vtkStringArray()
        arr.SetName(name)
        arr.SetNumberOfValues(len(values))
        for i, value in enumerate(values):
            arr.SetValue(i, str(value))
        field_data.AddArray(arr)

    @staticmethod
    def _add_int_array(field_data: vtk.vtkFieldData, name: str, values: np.ndarray) -> None:
        int_values = np.asarray(values, dtype=np.int64)
        vtk_arr = vtknp.numpy_to_vtk(int_values, deep=True, array_type=vtk.VTK_LONG_LONG)
        vtk_arr.SetName(name)
        vtk_arr.SetNumberOfComponents(1)
        field_data.AddArray(vtk_arr)

    @staticmethod
    def _add_double_array(
        field_data: vtk.vtkFieldData,
        name: str,
        values: np.ndarray,
        n_components: int = 1,
    ) -> None:
        np_values = np.asarray(values, dtype=float)
        if n_components > 1:
            np_values = np_values.reshape(-1, n_components)
        vtk_arr = vtknp.numpy_to_vtk(np_values, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(name)
        vtk_arr.SetNumberOfComponents(n_components)
        field_data.AddArray(vtk_arr)

    @staticmethod
    def _read_string_array(field_data: vtk.vtkFieldData, name: str) -> Optional[list[str]]:
        arr = vtk.vtkStringArray.SafeDownCast(field_data.GetAbstractArray(name))
        if arr is None:
            return None
        return [arr.GetValue(i) for i in range(arr.GetNumberOfValues())]

    @staticmethod
    def _read_numeric_array(field_data: vtk.vtkFieldData, name: str) -> Optional[np.ndarray]:
        arr = field_data.GetArray(name)
        if arr is None:
            return None
        out = np.asarray(vtknp.vtk_to_numpy(arr))
        if arr.GetNumberOfComponents() > 1:
            return out.reshape(-1, arr.GetNumberOfComponents())
        return out

    @staticmethod
    def _cast_dataframe(df: pd.DataFrame, type_map: Dict[str, Any]) -> pd.DataFrame:
        out = df.copy()
        for column, column_type in type_map.items():
            if column not in out.columns:
                continue
            if column_type is str:
                out[column] = out[column].astype(str)
            elif column_type is int:
                out[column] = pd.to_numeric(out[column], errors="coerce").fillna(0).astype(int)
            elif column_type is float:
                out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)
            elif column_type is bool:
                out[column] = out[column].astype(bool)
        return out

    def _write_dataframe_table(
        self,
        field_data: vtk.vtkFieldData,
        table_name: str,
        df: pd.DataFrame,
        type_map: Dict[str, Any],
    ) -> None:
        work_df = df.copy()
        for col in type_map.keys():
            if col not in work_df.columns:
                work_df[col] = np.nan
        work_df = work_df.loc[:, list(type_map.keys())]

        self._add_int_array(
            field_data,
            self._table_rows_array_name(table_name),
            np.array([len(work_df)], dtype=np.int64),
        )

        for col_name, col_type in type_map.items():
            array_name = self._table_column_array_name(table_name, col_name)
            series = work_df[col_name]
            if col_type is str:
                self._add_string_array(field_data, array_name, ["" if pd.isna(v) else str(v) for v in series.tolist()])
            elif col_type is int:
                values = pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=np.int64)
                self._add_int_array(field_data, array_name, values)
            else:
                values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
                self._add_double_array(field_data, array_name, values, n_components=1)

    def _read_dataframe_table(
        self,
        field_data: vtk.vtkFieldData,
        table_name: str,
        type_map: Dict[str, Any],
    ) -> pd.DataFrame:
        rows_arr = self._read_numeric_array(field_data, self._table_rows_array_name(table_name))
        if rows_arr is None or rows_arr.size == 0:
            return pd.DataFrame(columns=list(type_map.keys()))

        n_rows = int(rows_arr.reshape(-1)[0])
        data: Dict[str, Any] = {}
        for col_name, col_type in type_map.items():
            array_name = self._table_column_array_name(table_name, col_name)
            if col_type is str:
                values = self._read_string_array(field_data, array_name)
                if values is None:
                    values = [""] * n_rows
                data[col_name] = values[:n_rows]
            else:
                values = self._read_numeric_array(field_data, array_name)
                if values is None:
                    if col_type is int:
                        data[col_name] = np.zeros((n_rows,), dtype=np.int64)
                    else:
                        data[col_name] = np.full((n_rows,), np.nan, dtype=float)
                else:
                    data[col_name] = values.reshape(-1)[:n_rows]

        return self._cast_dataframe(pd.DataFrame(data), type_map)

    def _write_coordinate_systems(self, field_data: vtk.vtkFieldData) -> None:
        try:
            from .coordinate_system import COORDINATE_CLASS_TO_ID
        except Exception:
            COORDINATE_CLASS_TO_ID = {}

        ids = []
        class_ids = []
        origins = []

        for cs_id, cs_obj in self.coordinate_systems.items():
            class_name = cs_obj.__class__.__name__
            class_id = COORDINATE_CLASS_TO_ID.get(class_name)
            if class_id is None:
                continue
            origin = getattr(cs_obj, "origin", None)
            if origin is None:
                continue

            ids.append(int(cs_id))
            class_ids.append(int(class_id))
            origins.append(np.asarray(origin, dtype=float).reshape(3))

        self._add_int_array(field_data, f"{FIELD_PREFIX}COORD_IDS", np.asarray(ids, dtype=np.int64))
        self._add_int_array(field_data, f"{FIELD_PREFIX}COORD_CLASS_IDS", np.asarray(class_ids, dtype=np.int64))
        if len(origins) == 0:
            self._add_double_array(field_data, f"{FIELD_PREFIX}COORD_ORIGIN", np.empty((0, 3), dtype=float), n_components=3)
        else:
            self._add_double_array(field_data, f"{FIELD_PREFIX}COORD_ORIGIN", np.asarray(origins, dtype=float), n_components=3)

    def _read_coordinate_systems(self, field_data: vtk.vtkFieldData) -> Dict[int, Any]:
        try:
            from .coordinate_system import coordinate_class_from_id
        except Exception:
            return {}

        ids = self._read_numeric_array(field_data, f"{FIELD_PREFIX}COORD_IDS")
        class_ids = self._read_numeric_array(field_data, f"{FIELD_PREFIX}COORD_CLASS_IDS")
        origins = self._read_numeric_array(field_data, f"{FIELD_PREFIX}COORD_ORIGIN")

        if ids is None or class_ids is None or origins is None:
            return {}

        ids = ids.reshape(-1)
        class_ids = class_ids.reshape(-1)
        origins = origins.reshape(-1, 3)
        n_rows = min(len(ids), len(class_ids), len(origins))

        out: Dict[int, Any] = {}
        for i in range(n_rows):
            cls = coordinate_class_from_id(int(class_ids[i]))
            if cls is None:
                continue
            try:
                out[int(ids[i])] = cls(origin=origins[i])
            except Exception:
                continue

        return out

    def _write_bh_curves(self, field_data: vtk.vtkFieldData) -> None:
        material_ids = []
        point_idx = []
        hb_values = []

        for material_id, curve in self.bh_curves.items():
            curve_np = np.asarray(curve, dtype=float)
            if curve_np.ndim != 2 or curve_np.shape[1] != 2:
                continue
            for idx in range(curve_np.shape[0]):
                material_ids.append(int(material_id))
                point_idx.append(int(idx))
                hb_values.append(curve_np[idx])

        self._add_int_array(field_data, f"{FIELD_PREFIX}BH_MATERIAL_ID", np.asarray(material_ids, dtype=np.int64))
        self._add_int_array(field_data, f"{FIELD_PREFIX}BH_POINT_INDEX", np.asarray(point_idx, dtype=np.int64))
        if len(hb_values) == 0:
            self._add_double_array(field_data, f"{FIELD_PREFIX}BH_HB", np.empty((0, 2), dtype=float), n_components=2)
        else:
            self._add_double_array(field_data, f"{FIELD_PREFIX}BH_HB", np.asarray(hb_values, dtype=float), n_components=2)

    def _read_bh_curves(self, field_data: vtk.vtkFieldData) -> Dict[int, np.ndarray]:
        material_ids = self._read_numeric_array(field_data, f"{FIELD_PREFIX}BH_MATERIAL_ID")
        point_idx = self._read_numeric_array(field_data, f"{FIELD_PREFIX}BH_POINT_INDEX")
        hb_values = self._read_numeric_array(field_data, f"{FIELD_PREFIX}BH_HB")

        if material_ids is None or point_idx is None or hb_values is None:
            return {}

        material_ids = material_ids.reshape(-1)
        point_idx = point_idx.reshape(-1)
        hb_values = hb_values.reshape(-1, 2)

        n_rows = min(len(material_ids), len(point_idx), len(hb_values))
        grouped: Dict[int, list[tuple[int, np.ndarray]]] = {}
        for i in range(n_rows):
            mat_id = int(material_ids[i])
            grouped.setdefault(mat_id, []).append((int(point_idx[i]), hb_values[i]))

        out: Dict[int, np.ndarray] = {}
        for mat_id, pairs in grouped.items():
            pairs.sort(key=lambda p: p[0])
            out[mat_id] = np.asarray([p[1] for p in pairs], dtype=float)

        return out

    def _write_model_field_data(self) -> None:
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        field_data = self._mesh.GetFieldData()
        self._remove_prefixed_field_arrays(field_data)

        self._add_string_array(field_data, self._meta_array_name("MESH_NAME"), [self._mesh_name or ""])

        for table_name, type_map in TABLE_SPECS.items():
            df = getattr(self, table_name)
            self._write_dataframe_table(field_data, table_name, df, type_map)

        time_step_keys = list(DEFAULT_TIME_PARAMETERS_TYPES.keys())
        time_step_values = np.array(
            [float(self.time_step.get(key, DEFAULT_TIME_PARAMETERS[key])) for key in time_step_keys],
            dtype=float,
        ).reshape(1, -1)
        self._add_double_array(
            field_data,
            f"{FIELD_PREFIX}TIME_STEP",
            time_step_values,
            n_components=len(time_step_keys),
        )

        for key, key_type in DEFAULT_SOLVER_SETTINGS_TYPES.items():
            name = self._solver_setting_array_name(key)
            value = self.solver_settings.get(key, DEFAULT_SOLVER_SETTINGS[key])
            if key_type is str:
                self._add_string_array(field_data, name, [str(value)])
            elif key_type is int:
                self._add_int_array(field_data, name, np.array([int(value)], dtype=np.int64))
            elif key_type is bool:
                self._add_int_array(field_data, name, np.array([1 if bool(value) else 0], dtype=np.int64))
            else:
                self._add_double_array(field_data, name, np.array([float(value)], dtype=float))

        self._write_coordinate_systems(field_data)
        self._write_bh_curves(field_data)

    def _read_model_field_data(self) -> None:
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        field_data = self._mesh.GetFieldData()

        mesh_name_values = self._read_string_array(field_data, self._meta_array_name("MESH_NAME"))
        if mesh_name_values and mesh_name_values[0]:
            self._mesh_name = mesh_name_values[0]

        for table_name, type_map in TABLE_SPECS.items():
            setattr(self, table_name, self._read_dataframe_table(field_data, table_name, type_map))

        time_step_array = self._read_numeric_array(field_data, f"{FIELD_PREFIX}TIME_STEP")
        if time_step_array is not None and time_step_array.size >= len(DEFAULT_TIME_PARAMETERS_TYPES):
            values = time_step_array.reshape(-1)
            self.time_step = {
                "time": float(values[0]),
                "transient_effects": bool(round(values[1])),
                "theta": float(values[2]),
            }

        for key, key_type in DEFAULT_SOLVER_SETTINGS_TYPES.items():
            name = self._solver_setting_array_name(key)
            if key_type is str:
                value = self._read_string_array(field_data, name)
                if value and len(value) > 0:
                    self.solver_settings[key] = value[0]
            else:
                value = self._read_numeric_array(field_data, name)
                if value is None or value.size == 0:
                    continue
                raw = value.reshape(-1)[0]
                if key_type is int:
                    self.solver_settings[key] = int(raw)
                elif key_type is bool:
                    self.solver_settings[key] = bool(int(raw))
                else:
                    self.solver_settings[key] = float(raw)

        self.coordinate_systems = self._read_coordinate_systems(field_data)
        self.bh_curves = self._read_bh_curves(field_data)

    @classmethod
    def read_model(
        cls,
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = False,
    ) -> "APEC_M_Model":
        absolute_path = os.path.abspath(os.fspath(vtk_file_path))
        loaded_mesh = ModelEdit.read_mesh(absolute_path, reset_arrays=reset_arrays)
        resolved_name = os.path.splitext(os.path.basename(absolute_path))[0]
        return cls.from_grid(loaded_mesh, mesh_name=resolved_name)

    @classmethod
    def from_grid(
        cls,
        mesh: vtk.vtkUnstructuredGrid,
        mesh_name: Optional[str] = None,
    ) -> "APEC_M_Model":
        model = cls(mesh=mesh, mesh_name=mesh_name)
        model._read_model_field_data()
        return model

    def read_excel(
        self,
        excel_file_path: Union[str, PathLike[str]],
        data_to_read: Sequence[str] = tuple(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys()),
    ) -> None:
        if not isinstance(data_to_read, Sequence):
            raise TypeError("data_to_read must be a sequence of sheet names.")

        LOGGER.info("Reading Excel workbook '%s'.", excel_file_path)
        for label in data_to_read:
            try:
                self.read_excel_sheet(excel_file_path, label)
            except ValueError as exc:
                LOGGER.warning("%s", exc)

    def read_excel_sheet(
        self,
        excel_file: Union[str, PathLike[str]],
        sheet_name: str,
    ) -> None:
        if sheet_name not in SHEET_PROCESSING_MAP:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(SHEET_PROCESSING_MAP.keys())}."
            )

        processing_info = SHEET_PROCESSING_MAP[sheet_name]
        LOGGER.info("Reading sheet '%s' from '%s'.", sheet_name, excel_file)
        if isinstance(processing_info, str):
            func = getattr(Excel, processing_info)
            result = func(excel_file, sheet_name)
            if result is not None:
                setattr(self, sheet_name, result)
            return

        if self._mesh is None:
            raise RuntimeError(
                "No mesh is loaded. Load a mesh before reading mesh-dependent Excel sheets."
            )

        method_name, data_type = processing_info
        func = getattr(Excel, method_name)
        func(self.mesh, excel_file, sheet_name, data_type)

    def _save_model(self, output_file: Union[str, PathLike[str]], binary: bool = True) -> None:
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        self._write_model_field_data()
        ModelEdit.save_mesh(self.mesh, output_file, binary=binary)
