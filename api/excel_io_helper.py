import os
from os import PathLike
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Union

import numpy as np
import pandas as pd

from binary_helper import BinaryHelper
from constants import DEFAULT_COLUMNS_ATTRIBUTES_TYPES, VARIABLES_GLOBAL_LIST
from vtk_grid_helper import VTKGridHelper

if TYPE_CHECKING:  # pragma: no cover - used only for type hints
    import vtk


class ExcelIOHelper:
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

        # Determine the VTK data and mapping based on the data_type
        if data_type == "elem":
            vtk_data = mesh.GetCellData()
        elif data_type == 'node':
            vtk_data = mesh.GetPointData()
        else:
            raise ValueError("Invalid data_type. Must be 'elem' or 'node'.")

        # Retrieve the final_data NumPy array
        final_data = BinaryHelper.create_numpy_array_for_binary(vtk_data)

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
        print(f"{data_type.capitalize()} written to '{output_excel_file}'.")

    @staticmethod
    def read_excel_dict(excel_file: Union[str, PathLike[str]], sheet_name: str) -> Optional[Dict[str, Any]]:
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

        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
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
        print(f"Data from sheet '{sheet_name}' successfully loaded into dictionary.")
        return result_dict

    @staticmethod
    def read_excel_df(excel_file: Union[str, PathLike[str]], sheet_name: str) -> Optional[pd.DataFrame]:
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
        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
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
        print(f"Data from sheet '{sheet_name}' successfully loaded into DataFrame.")
        return df

    @staticmethod
    def _evaluate_periodic_matrix(raw_value: Any) -> np.ndarray:
        """
        Evaluate a periodic matrix expression into a NumPy array.

        Args:
            raw_value: Cell contents from the periodic matrix column.

        Returns:
            np.ndarray: Evaluated NumPy array.

        Raises:
            ValueError: If the expression cannot be evaluated or converted.
        """
        try:
            evaluated = eval(
                str(raw_value),
                {"np": np, "numpy": np, "__builtins__": {}},
            )
        except Exception as exc:
            raise ValueError(
                f"Failed to evaluate periodic matrix expression '{raw_value}': {exc}"
            ) from exc

        try:
            return np.asarray(evaluated, dtype=float)
        except Exception as exc:
            raise ValueError(
                f"Evaluated periodic matrix '{raw_value}' could not be converted to a NumPy array: {exc}"
            ) from exc

    @staticmethod
    def read_periodic_boundaries(
        excel_file: Union[str, PathLike[str]],
        sheet_name: str,
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Read periodic boundary definitions from an Excel sheet.

        The sheet must contain the columns defined in DEFAULT_COLUMNS_ATTRIBUTES_TYPES
        under 'periodic_boundaries': main_boundary, paired_boundary, periodic_matrix,
        and variable. The periodic_matrix column is evaluated via `eval` with NumPy
        available as `np`/`numpy`.

        Args:
            excel_file: Path to the Excel workbook.
            sheet_name: Name of the sheet containing periodic boundaries.

        Returns:
            Optional[List[Dict[str, Any]]]: Parsed periodic boundary entries, or None if the sheet is missing.

        Raises:
            FileNotFoundError: If the Excel file does not exist.
            ValueError: If required columns are missing or values cannot be coerced.
        """
        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")

        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
            return None

        df = df.dropna(how="all")
        expected_columns = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[sheet_name].keys()
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"Sheet '{sheet_name}' is missing required columns: {missing_columns}. "
                f"Expected columns: {list(expected_columns)}."
            )

        type_map = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[sheet_name]
        periodic_boundaries: List[Dict[str, Any]] = []
        for row_idx, row in df.iterrows():
            boundary: Dict[str, Any] = {}
            for column, typ in type_map.items():
                value = row[column]
                if pd.isna(value):
                    raise ValueError(
                        f"Value for column '{column}' is missing at row {row_idx + 2} in sheet '{sheet_name}'."
                    )

                if column == "periodic_matrix":
                    boundary[column] = ExcelIOHelper._evaluate_periodic_matrix(value)
                elif typ == int:
                    boundary[column] = int(value)
                elif typ == float:
                    boundary[column] = float(value)
                elif typ == str:
                    boundary[column] = str(value)
                else:
                    boundary[column] = value

            periodic_boundaries.append(boundary)

        print(f"Periodic boundaries loaded from sheet '{sheet_name}'.")
        return periodic_boundaries

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
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
            return

        if data_type is None:
            data_type = sheet_name
        # Determine the VTK data type to update
        if data_type == "elem":
            vtk_data = mesh.GetCellData()
        elif data_type == 'node':
            vtk_data = mesh.GetPointData()
        else:
            raise ValueError(
                f"Invalid data_type '{data_type}'. Must be 'elem' or 'node'."
            )

        VTKGridHelper.update_vtk_grid(vtk_data, df)
        print(f"VTK grid updated with data from sheet '{sheet_name}'.")
