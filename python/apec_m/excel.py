import os
from os import PathLike
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Union

import numpy as np
import pandas as pd

from .constants import DEFAULT_COLUMNS_ATTRIBUTES_TYPES, VARIABLES_GLOBAL_LIST
from .model_edit import ModelEdit
from .logging_utils import get_logger
from .validation_utils import (
    ensure_in_options,
    ensure_non_empty_string,
    ensure_pathlike,
    ensure_vtk_unstructured_grid,
)


LOGGER = get_logger("apec_m.excel")

if TYPE_CHECKING:  # pragma: no cover - used only for type hints
    import vtk


class Excel:
    """Excel helpers for reading and writing APEC-M tabular data."""

    @staticmethod
    def write_excel_cell_point_data(
        mesh: "vtk.vtkUnstructuredGrid",
        output_excel_file: Union[str, PathLike[str]],
        data_type: str = "elem",
    ) -> None:
        """
        Writes the loads or boundary conditions to an Excel file.

        Args:
            mesh: The VTK mesh object.
            output_excel_file: Path to the output Excel file.
            data_type: Data to write ("elem" or 'node'). Defaults to "elem".

        Returns:
            None

        Raises:
            ValueError: If the data_type is invalid.
        """
        ensure_vtk_unstructured_grid(mesh)
        output_excel_file = ensure_pathlike(output_excel_file, "output_excel_file")
        data_type = ensure_in_options(data_type, "data_type", ("elem", "node"))

        output_dir = os.path.dirname(output_excel_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Determine the VTK data and mapping based on the data_type
        vtk_data = mesh.GetCellData() if data_type == "elem" else mesh.GetPointData()

        # Retrieve the final_data as a structured DataFrame
        final_data = Excel._vtk_data_to_records(vtk_data)

        # Convert the NumPy array to a pandas DataFrame
        df = pd.DataFrame(final_data)

        # Replace the "load_bc_number" column with the corresponding "label"
        reverse_mapping = {i: v for i, v in enumerate(VARIABLES_GLOBAL_LIST)}
        df["label"] = df["load_bc_number"].map(reverse_mapping)

        # Reorder columns to make "label" the second column and remove "load_bc_number"
        df = df[["vtk_id", "label", "component", "value"]]

        # Check if the Excel file already exists
        # Overwrite the sheet if it exists
        if os.path.exists(output_excel_file):
            with pd.ExcelWriter(
                output_excel_file,
                engine="openpyxl",
                mode="a",
                if_sheet_exists="replace",
            ) as writer:
                df.to_excel(writer, sheet_name=data_type, index=False)
        else:
            with pd.ExcelWriter(output_excel_file, engine="openpyxl") as writer:
                df.to_excel(writer, sheet_name=data_type, index=False)
        LOGGER.info("%s written to '%s'.", data_type.capitalize(), output_excel_file)

    @staticmethod
    def _vtk_data_to_records(vtk_data: "vtk.vtkDataSetAttributes") -> pd.DataFrame:
        """Extract non-NaN values from named VTK arrays into a structured records table.

        Iterates over VARIABLES_GLOBAL_LIST and collects every non-NaN entry as a
        row with columns [vtk_id, load_bc_number, component, value].

        Args:
            vtk_data: VTK data attributes object (CellData or PointData).

        Returns:
            pd.DataFrame: Records with columns [vtk_id, load_bc_number, component, value].
        """
        import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

        records: List[Dict[str, Any]] = []
        for load_bc_number, var_name in enumerate(VARIABLES_GLOBAL_LIST):
            vtk_array = vtk_data.GetArray(var_name)
            if vtk_array is None:
                continue
            arr = np.asarray(vtknp.vtk_to_numpy(vtk_array), dtype=float)
            if arr.ndim == 1:
                valid_ids = np.flatnonzero(~np.isnan(arr))
                for vtk_id in valid_ids:
                    records.append(
                        {
                            "vtk_id": int(vtk_id),
                            "load_bc_number": load_bc_number,
                            "component": 0,
                            "value": float(arr[vtk_id]),
                        }
                    )
            else:
                valid_rows, valid_cols = np.where(~np.isnan(arr))
                for vtk_id, comp in zip(valid_rows.tolist(), valid_cols.tolist()):
                    records.append(
                        {
                            "vtk_id": int(vtk_id),
                            "load_bc_number": load_bc_number,
                            "component": int(comp),
                            "value": float(arr[vtk_id, comp]),
                        }
                    )
        return pd.DataFrame(
            records,
            columns=["vtk_id", "load_bc_number", "component", "value"],
        )

    @staticmethod
    def read_excel_dict(
        excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> Optional[Dict[str, Any]]:
        """
        Reads a specific sheet from an Excel file into a dictionary and returns it.
        Ensures that the keys in the 'SETTING' column match the default keys defined in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to read.

        Returns:
            Optional[Dict[str, Any]]: Loaded dictionary, or None if the sheet is missing.

        Raises:
            ValueError: If the sheet name is not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES or if the keys do not match.
            FileNotFoundError: If the Excel file does not exist.
        """

        excel_file = ensure_pathlike(excel_file, "excel_file")
        sheet_name = ensure_non_empty_string(sheet_name, "sheet_name")

        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            LOGGER.warning("Sheet '%s' not found in '%s'.", sheet_name, excel_file)
            return None
        required_columns = ["SETTING", "VALUE"]
        if not all(column in df.columns for column in required_columns):
            raise ValueError(
                f"Sheet '{sheet_name}' must contain the following columns: {required_columns}."
            )
        expected_keys = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[sheet_name].keys()
        invalid_keys = [key for key in df["SETTING"] if key not in expected_keys]
        if invalid_keys:
            raise ValueError(
                f"Invalid keys found in 'SETTING' column of sheet '{sheet_name}': {invalid_keys}. "
                f"Expected keys: {list(expected_keys)}."
            )
        type_map = DEFAULT_COLUMNS_ATTRIBUTES_TYPES.get(sheet_name, {})
        result_dict = {}
        for _, row in df.iterrows():
            setting = row["SETTING"]
            value = row["VALUE"]
            if setting in type_map:
                typ = type_map[setting]
                if typ == int:
                    value = int(value)
                elif typ == float:
                    value = float(value)
                elif typ == str:
                    value = str(value)
            result_dict[setting] = value
        LOGGER.info("Data from sheet '%s' successfully loaded into dictionary.", sheet_name)
        return result_dict

    @staticmethod
    def read_excel_df(
        excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> Optional[pd.DataFrame]:
        """
        Reads a specific sheet from an Excel file into a pandas DataFrame and returns it.
        Ensures that the columns in the sheet match the default columns defined in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to read.

        Returns:
            Optional[pd.DataFrame]: Loaded DataFrame, or None if the sheet is missing.

        Raises:
            ValueError: If the sheet name is not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES or if the columns do not match.
            FileNotFoundError: If the Excel file does not exist.
        """
        excel_file = ensure_pathlike(excel_file, "excel_file")
        sheet_name = ensure_non_empty_string(sheet_name, "sheet_name")

        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            LOGGER.warning("Sheet '%s' not found in '%s'.", sheet_name, excel_file)
            return None
        expected_columns = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[sheet_name]
        if not all(column in expected_columns for column in df.columns):
            raise ValueError(
                f"Columns in sheet '{sheet_name}' do not match the expected columns. "
                f"Expected: {expected_columns}, Found: {list(df.columns)}"
            )
        type_map = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[sheet_name]
        for col, typ in type_map.items():
            if col in df.columns:
                if typ == int:
                    df[col] = df[col].astype(int)
                elif typ == float:
                    df[col] = df[col].astype(float)
                elif typ == str:
                    df[col] = df[col].astype(str)
        LOGGER.info("Data from sheet '%s' successfully loaded into DataFrame.", sheet_name)
        return df

    @staticmethod
    def update_vtk_grid_from_excel(
        mesh: "vtk.vtkUnstructuredGrid",
        excel_file: Union[str, PathLike[str]],
        sheet_name: str,
        data_type: Optional[str] = None,
    ) -> None:
        """
        Reads a sheet from an Excel file and updates the VTK grid using the data.

        Args:
            mesh: The VTK mesh object.
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to read.
            data_type: Type of VTK data to update ('elem' or 'node'). Defaults to `sheet_name`.

        Returns:
            None

        Raises:
            ValueError: If the data_type is invalid or the sheet does not exist.
            FileNotFoundError: If the Excel file is missing.
        """
        ensure_vtk_unstructured_grid(mesh)
        excel_file = ensure_pathlike(excel_file, "excel_file")
        sheet_name = ensure_non_empty_string(sheet_name, "sheet_name")
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            LOGGER.warning("Sheet '%s' not found in '%s'.", sheet_name, excel_file)
            return

        required_columns = {"vtk_id", "label", "component", "value"}
        missing_columns = required_columns.difference(df.columns)
        if missing_columns:
            raise ValueError(
                f"Sheet '{sheet_name}' is missing required columns: {sorted(missing_columns)}."
            )

        if data_type is None:
            data_type = sheet_name
        data_type = ensure_in_options(data_type, "data_type", ("elem", "node"))
        # Determine the VTK data type to update
        if data_type == "elem":
            vtk_data = mesh.GetCellData()
        elif data_type == "node":
            vtk_data = mesh.GetPointData()

        ModelEdit.update_vtk_grid(vtk_data, df)
        LOGGER.info("VTK grid updated with data from sheet '%s'.", sheet_name)
