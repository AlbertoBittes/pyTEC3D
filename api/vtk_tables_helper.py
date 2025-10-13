import vtk
import numpy as np
import vtk.util.numpy_support as vtknp
import pandas as pd
from collections.abc import Iterable
import warnings


class VTKTablesHelper:
    @staticmethod
    def dataframe_to_vtk_table(df):
        """
        Converts a Pandas DataFrame into a vtkTable.

        Parameters:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            vtk.vtkTable: The converted vtkTable.
        """
        vtk_table = vtk.vtkTable()
        for column_name in df.columns:
            if pd.api.types.is_integer_dtype(df[column_name]):
                vtk_array = vtk.vtkIntArray()
            elif pd.api.types.is_float_dtype(df[column_name]):
                vtk_array = vtk.vtkDoubleArray()
            else:
                vtk_array = vtk.vtkStringArray()
            vtk_array.SetName(column_name)
            for value in df[column_name]:
                if isinstance(vtk_array, vtk.vtkIntArray):
                    vtk_array.InsertNextValue(int(value))
                elif isinstance(vtk_array, vtk.vtkDoubleArray):
                    vtk_array.InsertNextValue(float(value))
                else:
                    vtk_array.InsertNextValue(str(value))
            vtk_table.AddColumn(vtk_array)
        return vtk_table

    @staticmethod
    def dicts_to_vtk_table(dicts):
        """
        Converts a list (or single) of dictionaries into a vtkTable.
        Scalars are stored directly.
        Iterables (lists, tuples, numpy arrays) are flattened and stored with companion
        <key>_shape and <key>_offset arrays for reconstruction.

        Shape encoding (per iterable):
          [type_code, ndim, dim1, dim2, ..., optional element sizes for strings]

        Type codes:
            0 = int
            1 = float
            2 = str

        Parameters
        ----------
        dicts : Iterable[dict] or dict
            Input data dictionaries (all must share same keys).

        Returns
        -------
        vtk.vtkTable
            Resulting VTK table.
        """
        if isinstance(dicts, dict):
            dicts = [dicts]

        vtk_table = vtk.vtkTable()
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
                vtk_table.AddColumn(vtk_array)
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
                # Values
                if not flat_values:
                    continue

                if type_code in (0, 1):
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

                # Shape array
                vtk_shape = vtk.vtkIntArray()
                vtk_shape.SetName(f"{key}_shape")
                for s in shape_values:
                    vtk_shape.InsertNextValue(int(s))

                # Offset array
                vtk_offset = vtk.vtkIntArray()
                vtk_offset.SetName(f"{key}_offset")
                for o in offset_values:
                    vtk_offset.InsertNextValue(int(o))

                vtk_table.AddColumn(vtk_flat)
                vtk_table.AddColumn(vtk_shape)
                vtk_table.AddColumn(vtk_offset)
                continue

            warnings.warn(f"Unsupported type for key '{key}': {type(first_value)}")

        return vtk_table

    @staticmethod
    def add_table_to_mesh(mesh, vtk_table, overwrite=True):
        """
        Adds a vtkTable as field data to the given mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to add field data to.
            vtk_table (vtk.vtkTable): The vtkTable whose columns will be added as field data arrays.
            overwrite (bool, optional):
                - If True (default), removes any existing field data arrays with names matching
                the vtkTable columns before adding new ones (preserves previous behavior).
                - If False, retains existing arrays and only adds new ones if they
                do not already exist.

        Behavior:
            - Each column in vtk_table becomes a vtkAbstractArray in mesh.FieldData.
            - If overwrite=True, arrays with the same name are removed first (a warning is issued).
            - If overwrite=False, existing arrays are kept (a warning is issued and the column is skipped).

        Returns:
            None
        """
        field_data = mesh.GetFieldData()
        if not field_data:
            field_data = vtk.vtkFieldData()
            mesh.SetFieldData(field_data)

        for i in range(vtk_table.GetNumberOfColumns()):
            column_name = vtk_table.GetColumnName(i)
            vtk_array = vtk_table.GetColumn(i)

            if field_data.HasArray(column_name):
                if overwrite:
                    print(
                        f"Overwriting existing field data array '{column_name}'.",
                    )
                    field_data.RemoveArray(column_name)
                    field_data.AddArray(vtk_array)
                else:
                    print(
                        f"Field data array '{column_name}' already exists — skipping (overwrite=False).",
                    )
            else:
                field_data.AddArray(vtk_array)

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
