import vtk
import vtk.util.numpy_support as vtknp
import pandas as pd

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
    def dict_to_vtk_table(data_dict):
        """
        Converts a dictionary (keys: str, values: int/float/str) into a vtkTable.
        Each key becomes a column with a single value.

        Parameters:
            data_dict (dict): The input dictionary.

        Returns:
            vtk.vtkTable: The resulting vtkTable.
        """
        vtk_table = vtk.vtkTable()
        for key, value in data_dict.items():
            if isinstance(value, int):
                vtk_array = vtk.vtkIntArray()
            elif isinstance(value, float):
                vtk_array = vtk.vtkDoubleArray()
            else:
                vtk_array = vtk.vtkStringArray()
            vtk_array.SetName(str(key))
            vtk_array.InsertNextValue(value)
            vtk_table.AddColumn(vtk_array)
        return vtk_table

    @staticmethod
    def add_table_to_mesh(mesh, vtk_table):
        """
        Adds a vtkTable as field data to the given mesh.
        Deletes any existing field data arrays with names matching the columns of the vtkTable.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to add field data to.
            vtk_table (vtk.vtkTable): The vtkTable to add.
        """
        field_data = mesh.GetFieldData()
        if not field_data:
            field_data = vtk.vtkFieldData()
            mesh.SetFieldData(field_data)
        for i in range(vtk_table.GetNumberOfColumns()):
            column_name = vtk_table.GetColumnName(i)
            if field_data.HasArray(column_name):
                field_data.RemoveArray(column_name)
            field_data.AddArray(vtk_table.GetColumn(i))

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
