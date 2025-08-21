import numpy as np
import pandas as pd
import re
import vtk.util.numpy_support as vtknp

class AnsysHelper:
    @staticmethod
    def read_bfelist(mesh, file_path):
        """
        Reads and maps ANSYS element load data from a BFELIST file to the VTK mesh's CellData.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The VTK mesh object to update.
            file_path (str): Path to the BFELIST file containing element load data.

        Side Effects:
            Updates the mesh's CellData arrays with mapped element loads.
            Prints a summary message on completion.
        """
        df = AnsysHelper.read_parse_bfelist(file_path)
        mapping_array = AnsysHelper.map_data_to_vtk(mesh.GetCellData(), "ansys_elem_num")
        df = AnsysHelper.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_elem_num")
        from api.vtk_grid_helper import VTKGridHelper
        VTKGridHelper.update_vtk_grid(mesh.GetCellData(), df)
        print(f"Grid updated with element loads from '{file_path}'.")

    @staticmethod
    def convert_bfelist(mesh, bfelist_file, output_excel_file):
        """
        Reads a BFELIST file, updates the mesh's CellData, and writes the loads to an Excel file.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The VTK mesh object to update.
            bfelist_file (str): Path to the BFELIST file.
            output_excel_file (str): Path to the output Excel file.

        Side Effects:
            Updates mesh CellData and writes loads to Excel.
            Prints a summary message on completion.
        """
        import os
        if not os.path.exists(bfelist_file):
            raise FileNotFoundError(f"BFELIST file '{bfelist_file}' not found.")
        AnsysHelper.read_bfelist(mesh, bfelist_file)
        from api.excel_io_helper import ExcelIOHelper
        ExcelIOHelper.write_excel_cell_point_data(mesh, output_excel_file, data_type="cell_data")
        print(f"Loads from '{bfelist_file}' written to '{output_excel_file}'.")

    @staticmethod
    def read_dlist(mesh, file_path):
        """
        Reads and maps ANSYS node data from a DLIST file to the VTK mesh's PointData.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The VTK mesh object to update.
            file_path (str): Path to the DLIST file containing node data.

        Side Effects:
            Updates the mesh's PointData arrays with mapped node data.
            Prints a summary message on completion.
        """
        df = AnsysHelper.read_parse_dlist(file_path)
        mapping_array = AnsysHelper.map_data_to_vtk(mesh.GetPointData(), "ansys_node_num")
        df = AnsysHelper.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_node_num")
        from api.vtk_grid_helper import VTKGridHelper
        VTKGridHelper.update_vtk_grid(mesh.GetPointData(), df)
        print(f"Grid updated with node data from '{file_path}'.")

    @staticmethod
    def convert_dlist(mesh, dlist_file, output_excel_file):
        """
        Reads a DLIST file, updates the mesh's PointData, and writes the node data to an Excel file.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The VTK mesh object to update.
            dlist_file (str): Path to the DLIST file.
            output_excel_file (str): Path to the output Excel file.

        Side Effects:
            Updates mesh PointData and writes node data to Excel.
            Prints a summary message on completion.
        """
        import os
        if not os.path.exists(dlist_file):
            raise FileNotFoundError(f"DLIST file '{dlist_file}' not found.")
        AnsysHelper.read_dlist(mesh, dlist_file)
        from api.excel_io_helper import ExcelIOHelper
        ExcelIOHelper.write_excel_cell_point_data(mesh, output_excel_file, data_type="point_data")
        print(f"Node data from '{dlist_file}' written to '{output_excel_file}'.")

    @staticmethod
    def map_data_to_vtk(data, label):
        """
        Maps the ANSYS numbers (e.g., element or node numbers) to VTK IDs and returns a 2D NumPy array.
        This function assumes that the VTK grid has been loaded and contains the specified label.

        Parameters:
            data (vtk.vtkDataSetAttributes): The VTK data (e.g., CellData or PointData).
            label (str): The name of the array to map (e.g., "ansys_elem_num" or "ansys_node_num").

        Returns:
            np.ndarray: A 2D NumPy array where each row contains:
                        [VTK ID, ANSYS number].

        Raises:
            ValueError: If the specified label is not found in the provided data.
        """
        vtk_array = data.GetArray(label)
        if vtk_array is None:
            raise ValueError(f"Array '{label}' not found in the provided VTK data.")
        ansys_numbers = vtknp.vtk_to_numpy(vtk_array)
        vtk_ids = np.arange(len(ansys_numbers))
        mapping_array = np.column_stack((vtk_ids, ansys_numbers))
        return mapping_array

    @staticmethod
    def replace_ansys_with_vtk_ids(df, mapping_array, label):
        """
        Replaces the 'ansys_elem_num' column in the DataFrame with the corresponding VTK ID values
        using the provided mapping array. The column is renamed to 'vtk_element_id'.

        Args:
            df (pd.DataFrame): The DataFrame containing the 'ansys_elem_num' column.
            mapping_array (np.ndarray): A 2D NumPy array where each row contains [VTK ID, ANSYS element number].

        Returns:
            pd.DataFrame: A new DataFrame with the 'ansys_elem_num' column replaced by VTK IDs and renamed to 'vtk_element_id'.
        """
        ansys_numbers = df[label].to_numpy()
        vtk_ids = mapping_array[:, 0]
        ansys_mapping = mapping_array[:, 1]
        vtk_id_array = vtk_ids[np.searchsorted(ansys_mapping, ansys_numbers)]
        df.insert(0, "vtk_id", vtk_id_array)
        df = df.drop(columns=[label])
        return df

    @staticmethod
    def read_parse_bfelist(file_path):
        """
        Parses element load data from the given file and stores it in a pandas DataFrame.

        Args:
            file_path (str): Path to the file containing element load data.

        Returns:
            pd.DataFrame: A DataFrame with columns:
                - ansys_elem_num: The element ID.
                - load_label: The label associated with the load (e.g., "Js").
                - component: The component index (0, 1, 2).
                - value: The value of the load for the given component.
        """
        element_pattern = re.compile(r"ELEMENT\\s*=\\s*(\\d+)")
        density_pattern = re.compile(r"\\s*([\\d.\-E]+)\\s+([\\d.\-E]+)\\s+([\\d.\-E]+)\\s+([\\d.\-E]+)")
        data = []
        current_element = None
        load_label = None
        with open(file_path, "r") as file:
            for line in file:
                if line.strip():
                    if "CURRENT DENSITIES" in line:
                        load_label = "Js"
                    break
            for line in file:
                element_match = element_pattern.search(line)
                if element_match:
                    current_element = int(element_match.group(1))
                    continue
                if current_element is not None:
                    density_match = density_pattern.match(line)
                    if density_match:
                        for component in range(3):
                            value = float(density_match.group(component + 1))
                            if value != 0.0:
                                data.append({
                                    "ansys_elem_num": current_element,
                                    "label": load_label,
                                    "component": component,
                                    "value": value,
                                })
                        current_element = None
        df = pd.DataFrame(data, columns=["ansys_elem_num", "label", "component", "value"])
        return df

    @staticmethod
    def read_parse_dlist(file_path):
        """
        Parses node data from the given file and stores it in a pandas DataFrame.

        Args:
            file_path (str): Path to the file containing node data.

        Returns:
            pd.DataFrame: A DataFrame with columns:
                - ansys_node_num: The node ID.
                - variable_label: The label of the variable (e.g., "A", "VOLT").
                - component: The component index (0, 1, 2 for vector variables, 0 for scalar variables).
                - value: The value of the variable for the given component.
        """
        node_pattern = re.compile(r"^\\s*(\\d+)\\s+(\\S+)\\s+([\\d.\-E]+)")
        data = []
        with open(file_path, "r") as file:
            for line in file:
                if line.startswith("NODE") or line.startswith("//"):
                    continue
                node_match = node_pattern.match(line)
                if node_match:
                    node_num = int(node_match.group(1))
                    label = node_match.group(2)
                    value = float(node_match.group(3))
                    if label.endswith("X"):
                        component = 0
                        label = label[:-1]
                    elif label.endswith("Y"):
                        component = 1
                        label = label[:-1]
                    elif label.endswith("Z"):
                        component = 2
                        label = label[:-1]
                    else:
                        component = 0
                    data.append({
                        "ansys_node_num": node_num,
                        "label": label,
                        "component": component,
                        "value": value,
                    })
        df = pd.DataFrame(data, columns=["ansys_node_num", "label", "component", "value"])
        return df
