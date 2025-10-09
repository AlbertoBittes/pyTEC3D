import struct
import numpy as np
import vtk.util.numpy_support as vtknp
import pandas as pd
from constants import DEFAULT_COLUMNS_ATTRIBUTES_TYPES, VARIABLES_GLOBAL_LIST


class BinaryHelper:
    @staticmethod
    def _write_single_value(binary_file, value):
        if isinstance(value, int):
            binary_file.write(struct.pack("i", 0))  # type code for int
            binary_file.write(struct.pack("i", value))
        elif isinstance(value, float):
            binary_file.write(struct.pack("i", 1))  # type code for float
            binary_file.write(struct.pack("d", value))
        elif isinstance(value, str):
            binary_file.write(struct.pack("i", 2))  # type code for str
            encoded = value.encode("utf-8")
            binary_file.write(struct.pack("I", len(encoded)))
            binary_file.write(encoded)
        else:
            raise ValueError(f"Unsupported single value type: {type(value)}")

    @staticmethod
    def _write_numpy_array(binary_file, arr):
        # Write type code for numpy array
        binary_file.write(struct.pack("i", 3))
        # Write dtype as string
        dtype_str = str(arr.dtype)
        dtype_bytes = dtype_str.encode("utf-8")
        binary_file.write(struct.pack("I", len(dtype_bytes)))
        binary_file.write(dtype_bytes)
        # Write number of dimensions
        binary_file.write(struct.pack("I", arr.ndim))
        # Write shape
        for dim in arr.shape:
            binary_file.write(struct.pack("I", dim))
        # Write array data
        if arr.dtype.kind in {"i", "u"}:  # int, uint
            binary_file.write(struct.pack("i", 0))
            binary_file.write(arr.astype(np.int64).tobytes())
        elif arr.dtype.kind == "f":  # float
            binary_file.write(struct.pack("i", 1))
            binary_file.write(arr.astype(np.float64).tobytes())
        elif arr.dtype.kind == "U" or arr.dtype.kind == "S":  # unicode or bytes
            binary_file.write(struct.pack("i", 2))
            # For string arrays, write each string's length and bytes
            flat_arr = arr.flatten()
            binary_file.write(struct.pack("I", len(flat_arr)))
            for s in flat_arr:
                s_bytes = str(s).encode("utf-8")
                binary_file.write(struct.pack("I", len(s_bytes)))
                binary_file.write(s_bytes)
        else:
            raise ValueError(f"Unsupported numpy array dtype: {arr.dtype}")

    @staticmethod
    def write_settings_to_binary(settings, binary_file, label):
        """
        Writes a settings dictionary to an open binary file stream.
        Supports int, float, str, and numpy arrays as values.
        """
        if label not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES:
            raise ValueError(
                f"Invalid label '{label}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        expected_types = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[label]
        binary_file.write(struct.pack("I", len(settings)))
        for key, value in settings.items():
            binary_file.write(struct.pack("I", len(key)))
            binary_file.write(key.encode("utf-8"))
            expected_type = expected_types[key]
            # Validate type
            if isinstance(value, np.ndarray):
                # Optionally check dtype matches expected_type
                BinaryHelper._write_numpy_array(binary_file, value)
            elif expected_type == int and isinstance(value, int):
                BinaryHelper._write_single_value(binary_file, value)
            elif expected_type == float and isinstance(value, float):
                BinaryHelper._write_single_value(binary_file, value)
            elif expected_type == str and isinstance(value, str):
                BinaryHelper._write_single_value(binary_file, value)
            else:
                raise TypeError(
                    f"Key '{key}' is expected to be of type {expected_type.__name__} or numpy.ndarray, but got {type(value).__name__}."
                )

    @staticmethod
    def write_dataframe_to_binary(df, binary_file, key):
        """
        Writes a pandas DataFrame to an open binary file stream using predefined types for each column.

        Parameters:
            df (pd.DataFrame): The DataFrame to write.
            binary_file (file object): An open binary file stream.
            key (str): The key to retrieve the column types from DEFAULT_COLUMNS_ATTRIBUTES_TYPES.

        Raises:
            ValueError: If the key is invalid or if a column in the DataFrame does not have a predefined type.
        """
        if key not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES:
            raise ValueError(
                f"Invalid key '{key}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        column_types = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[key]
        for column in df.columns:
            if column not in column_types:
                raise ValueError(
                    f"Column '{column}' does not have a predefined type in column_types."
                )
            expected_type = column_types[column]
            if expected_type == int and not np.issubdtype(df[column].dtype, np.integer):
                raise TypeError(
                    f"Column '{column}' is expected to be of type int, but got {df[column].dtype}."
                )
            elif expected_type == float and not np.issubdtype(
                df[column].dtype, np.floating
            ):
                raise TypeError(
                    f"Column '{column}' is expected to be of type float, but got {df[column].dtype}."
                )
            elif expected_type == str and not np.issubdtype(
                df[column].dtype, np.object_
            ):
                raise TypeError(
                    f"Column '{column}' is expected to be of type str, but got {df[column].dtype}."
                )
        num_rows = len(df)
        num_columns = len(df.columns)
        binary_file.write(struct.pack("I", num_rows))
        binary_file.write(struct.pack("I", num_columns))
        for column in df.columns:
            column_name = column
            column_type = column_types[column]
            binary_file.write(struct.pack("I", len(column_name)))
            binary_file.write(column_name.encode("utf-8"))
            if column_type == int:
                binary_file.write(struct.pack("i", 0))
            elif column_type == float:
                binary_file.write(struct.pack("i", 1))
            elif column_type == str:
                binary_file.write(struct.pack("i", 2))
            else:
                raise ValueError(
                    f"Unsupported column type '{column_type}' for column '{column}'."
                )
        for _, row in df.iterrows():
            for column in df.columns:
                value = row[column]
                column_type = column_types[column]
                if column_type == int:
                    binary_file.write(struct.pack("i", int(value)))
                elif column_type == float:
                    binary_file.write(struct.pack("d", float(value)))
                elif column_type == str:
                    encoded_value = value.encode("utf-8")
                    binary_file.write(struct.pack("I", len(encoded_value)))
                    binary_file.write(encoded_value)
                else:
                    raise ValueError(
                        f"Unsupported column type '{column_type}' for column '{column}'."
                    )

    def write_any_to_binary(obj, binary_file, label):
        """
        Writes a dict, DataFrame, or an iterable of dicts/DataFrames to a binary file.

        Parameters:
            obj: dict, pd.DataFrame, or iterable of dicts/DataFrames
            binary_file: open binary file stream
            label: str, used for type validation (dicts) or column types (DataFrames)
        """
        if isinstance(obj, dict):
            BinaryHelper.write_settings_to_binary(obj, binary_file, label)
        elif isinstance(obj, pd.DataFrame):
            BinaryHelper.write_dataframe_to_binary(obj, binary_file, label)
        elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes)):
            obj_list = list(obj)
            binary_file.write(struct.pack("I", len(obj_list)))
            if not obj_list:
                return
            if all(isinstance(item, dict) for item in obj_list):
                for settings in obj_list:
                    BinaryHelper.write_settings_to_binary(settings, binary_file, label)
            elif all(isinstance(item, pd.DataFrame) for item in obj_list):
                for df in obj_list:
                    BinaryHelper.write_dataframe_to_binary(df, binary_file, label)
            else:
                raise TypeError("Iterable must contain only dicts or only DataFrames.")
        else:
            raise TypeError("Object must be a dict, DataFrame, or iterable of those.")

    @staticmethod
    def create_numpy_array_for_binary(cell_point_data):
        """
        Creates a NumPy structured array for loads/nodal DOFs from a VTK cell/point data array.

        Parameters:
            cell_point_data (vtk.vtkDataSetAttributes): The VTK cell or point data.

        Returns:
            np.ndarray: A structured NumPy array containing the loads or boundary conditions.
        """
        total_rows = 0
        all_stacked_data = []
        dtype = np.dtype(
            [
                ("vtk_id", np.uint32),
                ("load_bc_number", np.uint32),
                ("component", np.uint32),
                ("value", np.float64),
            ]
        )
        for idx, load_bc_label in enumerate(VARIABLES_GLOBAL_LIST):
            load_bc_array = cell_point_data.GetArray(load_bc_label)
            if load_bc_array is None:
                # print(f"Warning: Array '{load_bc_label}' not found in cell/point data.")
                continue
            load_bc_array_np = vtknp.vtk_to_numpy(load_bc_array)
            num_components = load_bc_array.GetNumberOfComponents()
            num_cells = cell_point_data.GetArray(0).GetNumberOfTuples()
            for component in range(num_components):
                if load_bc_array_np.ndim == 1:
                    component_values = load_bc_array_np
                else:
                    component_values = load_bc_array_np[:, component]
                mask = ~np.isnan(component_values)
                filtered_cell_points_ids = np.arange(num_cells, dtype=np.uint32)[mask]
                filtered_values = component_values[mask]
                load_bc_number_array = np.full(
                    len(filtered_cell_points_ids),
                    idx,
                    dtype=np.uint32,
                )
                component_array = np.full(
                    len(filtered_cell_points_ids), component, dtype=np.uint32
                )
                structured_data = np.zeros(len(filtered_cell_points_ids), dtype=dtype)
                structured_data["vtk_id"] = filtered_cell_points_ids
                structured_data["load_bc_number"] = load_bc_number_array
                structured_data["component"] = component_array
                structured_data["value"] = filtered_values
                all_stacked_data.append(structured_data)
                total_rows += len(filtered_cell_points_ids)
        if all_stacked_data:
            final_data = np.concatenate(all_stacked_data)
            return final_data
        return np.zeros(0, dtype=dtype)

    @staticmethod
    def write_numpy_array_to_binary(mesh, binary_file, label):
        """
        Extracts the correct numpy array from mesh (cell or point data) and writes it to an open binary file stream.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to extract data from.
            binary_file (file object): An open binary file stream.
            label (str): Either 'cell_data' or 'point_data'.

        Raises:
            ValueError: If the label is not recognized or the extracted array is not a NumPy array.
        """
        if label == "cell_data":
            numpy_array = BinaryHelper.create_numpy_array_for_binary(mesh.GetCellData())
        elif label == "point_data":
            numpy_array = BinaryHelper.create_numpy_array_for_binary(
                mesh.GetPointData()
            )
        else:
            raise ValueError(f"Unsupported label for numpy array binary write: {label}")
        if not isinstance(numpy_array, np.ndarray):
            raise ValueError("Extracted data is not a NumPy array.")
        binary_file.write(struct.pack("I", len(numpy_array)))
        binary_file.write(numpy_array.tobytes())
