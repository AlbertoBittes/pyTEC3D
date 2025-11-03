import struct
import numpy as np
import vtk.util.numpy_support as vtknp
import pandas as pd
from constants import VARIABLES_GLOBAL_LIST


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
    def _dataframe_to_column_dict(df):
        column_dict = {}
        for column_name in df.columns:
            series = df[column_name]
            if pd.api.types.is_integer_dtype(series):
                column_dict[column_name] = series.to_numpy(dtype=np.int64, copy=True)
            elif pd.api.types.is_float_dtype(series):
                column_dict[column_name] = series.to_numpy(dtype=np.float64, copy=True)
            elif pd.api.types.is_bool_dtype(series):
                column_dict[column_name] = series.astype(np.int64).to_numpy(copy=True)
            else:
                column_dict[column_name] = series.astype(str).tolist()
        return column_dict

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
    def write_data_to_binary(settings, binary_file, label=None):
        """
        Writes one or many settings dictionaries to an open binary file stream.
        Values are flattened into NDArray payloads so they can be reconstructed later.
        """
        settings_list = BinaryHelper._ensure_settings_list(settings)
        if not settings_list:
            binary_file.write(struct.pack("I", 0))
            return

        keys = BinaryHelper._validate_settings_keys(settings_list)
        binary_file.write(struct.pack("I", len(keys)))

        for key in keys:
            binary_file.write(struct.pack("I", len(key)))
            binary_file.write(key.encode("utf-8"))
            values = [settings[key] for settings in settings_list]
            BinaryHelper._write_flattened_series(binary_file, values, key)

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
            obj_list = [BinaryHelper._dataframe_to_column_dict(obj)]
        elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes)):
            obj_list = []
            for item in obj:
                if isinstance(item, dict):
                    obj_list.append(item)
                elif isinstance(item, pd.DataFrame):
                    obj_list.append(BinaryHelper._dataframe_to_column_dict(item))
                else:
                    raise TypeError("Iterable must contain only dicts or DataFrames.")
        else:
            raise TypeError("Object must be a dict, DataFrame, or iterable of those.")

        binary_file.write(struct.pack("I", len(obj_list)))
        if not obj_list:
            return

        BinaryHelper.write_data_to_binary(obj_list, binary_file, label)

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
