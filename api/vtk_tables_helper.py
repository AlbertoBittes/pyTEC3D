import vtk
import numpy as np
import vtk.util.numpy_support as vtknp
import pandas as pd
from collections.abc import Iterable
import warnings


class VTKTablesHelper:
    @staticmethod
    def dataframe_to_vtk_field_arrays(mesh, df, overwrite=True):
        """
        Converts a Pandas DataFrame into vtkArrays and adds them to the mesh's field data.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to add field data to.
            df (pd.DataFrame): The input DataFrame.
            overwrite (bool): If True, overwrite arrays with the same name.

        Returns:
            None
        """
        arrays = []
        for column_name in df.columns:
            col = df[column_name]
            if pd.api.types.is_integer_dtype(col):
                vtk_array = vtk.vtkIntArray()
            elif pd.api.types.is_float_dtype(col):
                vtk_array = vtk.vtkDoubleArray()
            else:
                vtk_array = vtk.vtkStringArray()
            vtk_array.SetName(column_name)
            for value in col:
                if isinstance(vtk_array, vtk.vtkIntArray):
                    vtk_array.InsertNextValue(int(value))
                elif isinstance(vtk_array, vtk.vtkDoubleArray):
                    vtk_array.InsertNextValue(float(value))
                else:
                    vtk_array.InsertNextValue(str(value))
            arrays.append(vtk_array)
        VTKTablesHelper.add_field_arrays_to_mesh(mesh, arrays, overwrite=overwrite)

    @staticmethod
    def dicts_to_vtk_field_arrays(mesh, dicts, overwrite=True):
        """
        Converts a list (or single) of dictionaries into vtkArrays and adds them to the mesh's field data.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to add field data to.
            dicts (dict or list of dict): The input data.
            overwrite (bool): If True, overwrite arrays with the same name.

        Returns:
            None
        """
        if isinstance(dicts, dict):
            dicts = [dicts]

        arrays = []
        keys = dicts[0].keys()
        for key in keys:
            first_value = dicts[0][key]

            # ---- SCALAR FIELDS ----
            if isinstance(first_value, (int, float, str)):
                if isinstance(first_value, int):
                    vtk_array = vtk.vtkIntArray()
                elif isinstance(first_value, float):
                    vtk_array = vtk.vtkDoubleArray()
                else:
                    vtk_array = vtk.vtkStringArray()

                vtk_array.SetName(str(key))
                for d in dicts:
                    vtk_array.InsertNextValue(d[key])
                arrays.append(vtk_array)
                continue

            # ---- ITERABLE FIELDS ----
            if isinstance(first_value, Iterable) and not isinstance(
                first_value, (str, bytes)
            ):
                flat_values = []
                shape_values = []
                offset_values = []

                cursor_shape = 0
                for d in dicts:
                    arr = np.asarray(d[key], dtype=object)
                    ndim = arr.ndim
                    shape = list(arr.shape)
                    elem0 = arr.flat[0] if arr.size > 0 else ""

                    # Determine type code
                    if isinstance(elem0, (int, np.integer)):
                        type_code = 0
                        arr_flat = arr.flatten().astype(int).tolist()
                    elif isinstance(elem0, (float, np.floating)):
                        type_code = 1
                        arr_flat = arr.flatten().astype(float).tolist()
                    elif isinstance(elem0, str):
                        type_code = 2
                        arr_flat = arr.flatten().tolist()
                    else:
                        warnings.warn(
                            f"Unsupported element type for key '{key}': {type(elem0)}"
                        )
                        continue

                    flat_values.extend(arr_flat)

                    # Store offset for this row
                    offset_values.append(cursor_shape)

                    # Shape encoding: [type_code, ndim, dim1, dim2, ..., optional string lengths]
                    shape_entry = [type_code, ndim] + shape
                    if type_code == 2:
                        shape_entry += [len(s) for s in arr_flat]
                    shape_values.extend(shape_entry)
                    cursor_shape += len(shape_entry)

                # ---- CREATE VTK ARRAYS ----
                if not flat_values:
                    continue

                if type_code == 0:
                    vtk_flat = vtk.vtkIntArray()
                    for v in flat_values:
                        vtk_flat.InsertNextValue(int(v))
                elif type_code == 1:
                    vtk_flat = vtk.vtkDoubleArray()
                    for v in flat_values:
                        vtk_flat.InsertNextValue(float(v))
                elif type_code == 2:
                    vtk_flat = vtk.vtkStringArray()
                    for v in flat_values:
                        vtk_flat.InsertNextValue(v)
                else:
                    warnings.warn(f"Unsupported flattened data type for key '{key}'")
                    continue

                vtk_flat.SetName(str(key))
                vtk_shape = vtk.vtkIntArray()
                vtk_shape.SetName(f"{key}_shape")
                for s in shape_values:
                    vtk_shape.InsertNextValue(int(s))
                vtk_offset = vtk.vtkIntArray()
                vtk_offset.SetName(f"{key}_offset")
                for o in offset_values:
                    vtk_offset.InsertNextValue(int(o))

                arrays.extend([vtk_flat, vtk_shape, vtk_offset])
                continue

            warnings.warn(f"Unsupported type for key '{key}': {type(first_value)}")

        VTKTablesHelper.add_field_arrays_to_mesh(mesh, arrays, overwrite=overwrite)

    @staticmethod
    def add_field_arrays_to_mesh(mesh, vtk_arrays, overwrite=True):
        """
        Adds a list of vtkArrays to the mesh's field data.
        If an array with the same name exists, it is overwritten (remove + add) if overwrite=True,
        otherwise it is skipped.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to add field data to.
            vtk_arrays (list): List of vtkArrays to add.
            overwrite (bool): If True, overwrite arrays with the same name.

        Returns:
            None
        """
        field_data = mesh.GetFieldData()
        if not field_data:
            field_data = vtk.vtkFieldData()
            mesh.SetFieldData(field_data)

        for arr in vtk_arrays:
            name = arr.GetName()
            if field_data.HasArray(name):
                if overwrite:
                    field_data.RemoveArray(name)
                    field_data.AddArray(arr)
                else:
                    continue
            else:
                field_data.AddArray(arr)

    @staticmethod
    def get_dataframe_from_field_data(mesh):
        """
        Extracts a pandas DataFrame from the field data of the given mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to extract field data from.

        Returns:
            pd.DataFrame: A DataFrame containing the field data.

        Raises:
            ValueError: If the field data is missing in the mesh.
        """
        field_data = mesh.GetFieldData()
        if not field_data:
            raise ValueError("Field data is missing in the mesh.")
        data = {}
        for i in range(field_data.GetNumberOfArrays()):
            vtk_array = field_data.GetArray(i)
            column_name = vtk_array.GetName()
            data[column_name] = vtknp.vtk_to_numpy(vtk_array)
        df = pd.DataFrame(data)
        return df

    @staticmethod
    def update_mesh_field_data(mesh, data, overwrite=True):
        """
        Updates the mesh field data using input data (dict(s) or DataFrame).

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to update.
            data (dict, list of dict, or pd.DataFrame): The input data.
            overwrite (bool): If True, overwrite arrays with the same name.

        Returns:
            None
        """
        if isinstance(data, pd.DataFrame):
            VTKTablesHelper.dataframe_to_vtk_field_arrays(
                mesh, data, overwrite=overwrite
            )
        elif isinstance(data, dict) or (
            isinstance(data, list) and all(isinstance(d, dict) for d in data)
        ):
            VTKTablesHelper.dicts_to_vtk_field_arrays(mesh, data, overwrite=overwrite)
        else:
            raise TypeError(
                "Input data must be a pandas DataFrame, a dict, or a list of dicts."
            )
