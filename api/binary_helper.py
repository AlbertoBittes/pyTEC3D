import struct
import numpy as np
import vtk.util.numpy_support as vtknp
import pandas as pd
from constants import DEFAULT_COLUMNS_ATTRIBUTES_TYPES, VARIABLES_GLOBAL_LIST


class BinaryHelper:

    @staticmethod
    def _ensure_settings_list(settings):
        if isinstance(settings, dict):
            return [settings]
        if hasattr(settings, "__iter__") and not isinstance(settings, (str, bytes)):
            settings_list = list(settings)
            for idx, item in enumerate(settings_list):
                if not isinstance(item, dict):
                    raise TypeError(
                        f"Expected dict at index {idx} in settings iterable, got {type(item).__name__}."
                    )
            return settings_list
        raise TypeError("Settings must be a dict or an iterable of dicts.")

    @staticmethod
    def _validate_settings_keys(settings_list):
        first_keys = list(settings_list[0].keys())
        reference = set(first_keys)
        for idx, settings in enumerate(settings_list[1:], start=1):
            keys_set = set(settings.keys())
            if keys_set != reference:
                raise ValueError(
                    f"Settings dict at index {idx} does not match keys of the first dict."
                )
        return first_keys

    @staticmethod
    def _validate_series_type(values, expected_type, key):
        if expected_type == int:
            ok = all(isinstance(v, (int, np.integer)) for v in values)
        elif expected_type == float:
            ok = all(
                isinstance(v, (float, int, np.floating, np.integer)) for v in values
            )
        elif expected_type == str:
            ok = all(isinstance(v, str) for v in values)
        else:
            ok = True
        if not ok:
            raise TypeError(
                f"Values for key '{key}' do not match expected type {expected_type.__name__}."
            )

    @staticmethod
    def _write_flattened_series(binary_file, values, key):
        if all(isinstance(v, (int, np.integer)) for v in values):
            data = np.asarray(values, dtype=np.int64)
            BinaryHelper._write_ndarray_payload(
                binary_file, "int64", [len(values)], 0, data
            )
            return

        if all(isinstance(v, (int, float, np.integer, np.floating)) for v in values):
            data = np.asarray(values, dtype=np.float64)
            BinaryHelper._write_ndarray_payload(
                binary_file, "float64", [len(values)], 1, data
            )
            return

        if all(isinstance(v, str) for v in values):
            BinaryHelper._write_ndarray_payload(
                binary_file, "str", [len(values)], 2, list(values)
            )
            return

        arrays = []
        for idx, value in enumerate(values):
            if isinstance(value, np.ndarray):
                arrays.append(value)
            elif isinstance(value, (list, tuple)):
                arrays.append(np.asarray(value))
            else:
                raise TypeError(
                    f"Key '{key}' contains unsupported value type '{type(value).__name__}' at index {idx}."
                )

        dtype_kind = arrays[0].dtype.kind
        if dtype_kind in {"i", "u"}:
            element_kind = 0
            dtype_str = "int64"
        elif dtype_kind == "f":
            element_kind = 1
            dtype_str = "float64"
        elif dtype_kind in {"U", "S", "O"}:
            element_kind = 2
            dtype_str = "str"
        else:
            raise TypeError(
                f"Key '{key}' contains numpy arrays with unsupported dtype kind '{dtype_kind}'."
            )

        base_shape = arrays[0].shape
        for idx, arr in enumerate(arrays[1:], start=1):
            if arr.shape != base_shape:
                raise ValueError(
                    f"Key '{key}' expects consistent array shapes; got {base_shape} and {arr.shape}."
                )

        full_shape = [len(arrays), *base_shape]
        if element_kind in {0, 1}:
            if arrays and arrays[0].size:
                flattened = np.concatenate([arr.reshape(-1) for arr in arrays])
            else:
                flattened = np.empty(0, dtype=arrays[0].dtype)
            cast_dtype = np.int64 if element_kind == 0 else np.float64
            BinaryHelper._write_ndarray_payload(
                binary_file, dtype_str, full_shape, element_kind, flattened.astype(cast_dtype)
            )
        else:
            str_values = []
            for arr in arrays:
                flattened = arr.reshape(-1)
                str_values.extend(str(item) for item in flattened)
            BinaryHelper._write_ndarray_payload(
                binary_file, dtype_str, full_shape, element_kind, str_values
            )

    @staticmethod
    def _write_ndarray_payload(binary_file, dtype_str, shape, element_kind, data):
        binary_file.write(struct.pack("i", 3))

        dtype_bytes = dtype_str.encode("utf-8")
        binary_file.write(struct.pack("I", len(dtype_bytes)))
        binary_file.write(dtype_bytes)

        binary_file.write(struct.pack("I", len(shape)))
        for dim in shape:
            binary_file.write(struct.pack("I", int(dim)))

        binary_file.write(struct.pack("i", element_kind))

        total_elements = 1
        for dim in shape:
            total_elements *= int(dim)

        if element_kind in {0, 1}:
            if hasattr(data, "size") and data.size != total_elements:
                raise ValueError("NDArray payload size mismatch for numeric data.")
            if total_elements:
                binary_file.write(data.tobytes())
            return

        if len(data) != total_elements:
            raise ValueError("NDArray payload size mismatch for string data.")
        binary_file.write(struct.pack("I", total_elements))
        for value in data:
            encoded = str(value).encode("utf-8")
            binary_file.write(struct.pack("I", len(encoded)))
            binary_file.write(encoded)

    @staticmethod
    def write_settings_to_binary(settings, binary_file, label):
        """
        Writes one or many settings dictionaries to an open binary file stream.
        For labels without predefined column types, values are flattened across the
        list into NDArray payloads so they can be reconstructed later.
        """
        settings_list = BinaryHelper._ensure_settings_list(settings)
        if not settings_list:
            binary_file.write(struct.pack("I", 0))
            return

        keys = BinaryHelper._validate_settings_keys(settings_list)
        binary_file.write(struct.pack("I", len(keys)))
        expected_types = DEFAULT_COLUMNS_ATTRIBUTES_TYPES.get(label)

        for key in keys:
            binary_file.write(struct.pack("I", len(key)))
            binary_file.write(key.encode("utf-8"))
            values = [settings[key] for settings in settings_list]
            if expected_types is not None and key in expected_types:
                BinaryHelper._validate_series_type(values, expected_types[key], key)
            BinaryHelper._write_flattened_series(binary_file, values, key)

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

    @staticmethod
    def write_to_binary(obj, binary_file, label):
        """
        Writes a dict, DataFrame, or an iterable of dicts/DataFrames to a binary file.
        A single dict/DataFrame is treated as an iterable containing that single element.

        Parameters:
            obj: dict, pd.DataFrame, or iterable of dicts/DataFrames
            binary_file: open binary file stream
            label: str, used for type validation (dicts) or column types (DataFrames)
        """
        if isinstance(obj, dict):
            obj_list = [obj]
        elif isinstance(obj, pd.DataFrame):
            obj_list = [obj]
        elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes)):
            obj_list = list(obj)
        else:
            raise TypeError("Object must be a dict, DataFrame, or iterable of those.")

        binary_file.write(struct.pack("I", len(obj_list)))
        if not obj_list:
            return

        if all(isinstance(item, dict) for item in obj_list):
            BinaryHelper.write_settings_to_binary(obj_list, binary_file, label)
        elif all(isinstance(item, pd.DataFrame) for item in obj_list):
            for df in obj_list:
                BinaryHelper.write_dataframe_to_binary(df, binary_file, label)
        else:
            raise TypeError("Iterable must contain only dicts or only DataFrames.")

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
