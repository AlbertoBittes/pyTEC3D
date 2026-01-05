import os
import re
from os import PathLike
from typing import TYPE_CHECKING, Union

import numpy as np
import pandas as pd
import vtk.util.numpy_support as vtknp

from mesh_edit import MeshEdit
from excel_io_helper import ExcelIOHelper

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtk


class AnsysHelper:
    @staticmethod
    def read_bfelist(
        mesh: "vtk.vtkUnstructuredGrid", file_path: Union[str, PathLike[str]]
    ) -> None:
        """
        Reads and maps ANSYS element load data from a BFELIST file to the VTK mesh's CellData.

        Args:
            mesh: VTK mesh object to update (CellData will be modified).
            file_path: Path to the BFELIST file containing element load data.

        Returns:
            None

        Raises:
            FileNotFoundError: If the BFELIST file does not exist.
            ValueError: If required ANSYS identifiers are missing in the mesh CellData.

        Side Effects:
            Updates the mesh's CellData arrays and prints a summary message.
        """
        df = AnsysHelper.read_parse_bfelist(file_path)
        mapping_array = AnsysHelper.map_data_to_vtk(
            mesh.GetCellData(), "ansys_elem_num"
        )
        df = AnsysHelper.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_elem_num")
        MeshEdit.update_vtk_grid(mesh.GetCellData(), df)
        print(f"Grid updated with element loads from '{file_path}'.")

    @staticmethod
    def convert_bfelist(
        mesh: "vtk.vtkUnstructuredGrid",
        bfelist_file: Union[str, PathLike[str]],
        output_excel_file: Union[str, PathLike[str]],
    ) -> None:
        """
        Reads a BFELIST file, updates the mesh's CellData, and writes the loads to an Excel file.

        Args:
            mesh: VTK mesh object to update (CellData will be modified).
            bfelist_file: Path to the BFELIST file.
            output_excel_file: Path to the output Excel file to create.

        Returns:
            None

        Raises:
            FileNotFoundError: If the BFELIST file does not exist.
            ValueError: If required ANSYS identifiers are missing in the mesh CellData.

        Side Effects:
            Updates mesh CellData, writes loads to Excel, and prints a summary message.
        """

        if not os.path.exists(bfelist_file):
            raise FileNotFoundError(f"BFELIST file '{bfelist_file}' not found.")
        AnsysHelper.read_bfelist(mesh, bfelist_file)

        ExcelIOHelper.write_excel_cell_point_data(
            mesh, output_excel_file, data_type="elem"
        )
        print(f"Loads from '{bfelist_file}' written to '{output_excel_file}'.")

    @staticmethod
    def read_dlist(
        mesh: "vtk.vtkUnstructuredGrid", file_path: Union[str, PathLike[str]]
    ) -> None:
        """
        Reads and maps ANSYS node data from a DLIST file to the VTK mesh's PointData.

        Args:
            mesh: VTK mesh object to update (PointData will be modified).
            file_path: Path to the DLIST file containing node data.

        Returns:
            None

        Raises:
            FileNotFoundError: If the DLIST file does not exist.
            ValueError: If required ANSYS identifiers are missing in the mesh PointData.

        Side Effects:
            Updates the mesh's PointData arrays and prints a summary message.
        """
        df = AnsysHelper.read_parse_dlist(file_path)
        mapping_array = AnsysHelper.map_data_to_vtk(
            mesh.GetPointData(), "ansys_node_num"
        )
        df = AnsysHelper.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_node_num")

        MeshEdit.update_vtk_grid(mesh.GetPointData(), df)
        print(f"Grid updated with node data from '{file_path}'.")

    @staticmethod
    def convert_dlist(
        mesh: "vtk.vtkUnstructuredGrid",
        dlist_file: Union[str, PathLike[str]],
        output_excel_file: Union[str, PathLike[str]],
    ) -> None:
        """
        Reads a DLIST file, updates the mesh's PointData, and writes the node data to an Excel file.

        Args:
            mesh: VTK mesh object to update (PointData will be modified).
            dlist_file: Path to the DLIST file.
            output_excel_file: Path to the output Excel file to create.

        Returns:
            None

        Raises:
            FileNotFoundError: If the DLIST file does not exist.
            ValueError: If required ANSYS identifiers are missing in the mesh PointData.

        Side Effects:
            Updates mesh PointData, writes node data to Excel, and prints a summary message.
        """
        if not os.path.exists(dlist_file):
            raise FileNotFoundError(f"DLIST file '{dlist_file}' not found.")
        AnsysHelper.read_dlist(mesh, dlist_file)

        ExcelIOHelper.write_excel_cell_point_data(
            mesh, output_excel_file, data_type="node"
        )
        print(f"Node data from '{dlist_file}' written to '{output_excel_file}'.")

    @staticmethod
    def map_data_to_vtk(data: "vtk.vtkDataSetAttributes", label: str) -> np.ndarray:
        """
        Maps the ANSYS numbers (e.g., element or node numbers) to VTK IDs and returns a 2D NumPy array.
        This function assumes that the VTK grid has been loaded and contains the specified label.

        Args:
            data: The VTK data (e.g., CellData or PointData).
            label: Name of the array to map (e.g., "ansys_elem_num" or "ansys_node_num").

        Returns:
            np.ndarray: 2D array where each row contains [VTK ID, ANSYS number].

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
    def replace_ansys_with_vtk_ids(
        df: pd.DataFrame, mapping_array: np.ndarray, label: str
    ) -> pd.DataFrame:
        """
        Replace an ANSYS id column in a DataFrame with the corresponding VTK IDs.

        Args:
            df: DataFrame containing the ANSYS id column.
            mapping_array: 2D array where each row contains [VTK ID, ANSYS number].
            label: Column name in `df` holding the ANSYS numbers to replace.

        Returns:
            pd.DataFrame: Copy with ANSYS numbers replaced by VTK IDs in a `vtk_id` column.
        """
        ansys_numbers = df[label].to_numpy()
        vtk_ids = mapping_array[:, 0]
        ansys_mapping = mapping_array[:, 1]
        vtk_id_array = vtk_ids[np.searchsorted(ansys_mapping, ansys_numbers)]
        df.insert(0, "vtk_id", vtk_id_array)
        df = df.drop(columns=[label])
        return df

    @staticmethod
    def read_parse_bfelist(file_path: Union[str, PathLike[str]]) -> pd.DataFrame:
        """
        Parses element load data from the given file and stores it in a pandas DataFrame.

        Args:
            file_path: Path to the file containing element load data.

        Returns:
            pd.DataFrame: A DataFrame with columns:
                - ansys_elem_num: The element ID.
                - label: The label associated with the load (e.g., "Js").
                - component: The component index (0, 1, 2).
                - value: The value of the load for the given component.

        Raises:
            FileNotFoundError: If `file_path` cannot be opened.
        """
        # Regular expressions to match element and density lines
        element_pattern = re.compile(r"ELEMENT\s*=\s*(\d+)")
        density_pattern = re.compile(
            r"\s*([\d\.\-E]+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)"
        )

        # Initialize variables
        data = []
        current_element = None
        load_label = None

        with open(file_path, "r") as file:
            # Skip empty lines and find the line with "CURRENT DENSITIES"
            for line in file:
                if line.strip():  # Skip empty lines
                    if "CURRENT DENSITIES" in line:
                        load_label = "Js"
                    break  # Exit the loop after finding the first non-empty line

            # Parse the rest of the file
            for line in file:
                # Match element line
                element_match = element_pattern.search(line)
                if element_match:
                    current_element = int(element_match.group(1))
                    continue

                # Match density line
                if current_element is not None:
                    density_match = density_pattern.match(line)
                    if density_match:
                        for component in range(3):  # Components 0, 1, 2
                            value = float(density_match.group(component + 1))
                            if value != 0.0:
                                data.append(
                                    {
                                        "ansys_elem_num": current_element,
                                        "label": load_label,
                                        "component": component,
                                        "value": value,
                                    }
                                )
                        current_element = None  # Reset current element after processing

        # Convert the data to a pandas DataFrame
        df = pd.DataFrame(
            data, columns=["ansys_elem_num", "label", "component", "value"]
        )
        return df

    @staticmethod
    def read_parse_dlist(file_path: Union[str, PathLike[str]]) -> pd.DataFrame:
        """
        Parses node data from the given file and stores it in a pandas DataFrame.

        Args:
            file_path: Path to the file containing node data.

        Returns:
            pd.DataFrame: A DataFrame with columns:
                - ansys_node_num: The node ID.
                - label: The label of the variable (e.g., "A", "VOLT").
                - component: The component index (0, 1, 2 for vector variables, 0 for scalar variables).
                - value: The value of the variable for the given component.

        Raises:
            FileNotFoundError: If `file_path` cannot be opened.
        """
        # Regular expression to match node data lines
        node_pattern = re.compile(r"^\s*(\d+)\s+(\S+)\s+([\d\.\-E]+)")

        # List to store parsed data
        data = []

        with open(file_path, "r") as file:
            for line in file:
                # Skip header lines
                if line.startswith("NODE") or line.startswith("//"):
                    continue

                # Match node data line
                node_match = node_pattern.match(line)
                if node_match:
                    node_num = int(node_match.group(1))
                    label = node_match.group(2)
                    value = float(node_match.group(3))

                    # Determine the component based on the label
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
                        component = 0  # Scalar variable

                    # Append the parsed data to the list
                    data.append(
                        {
                            "ansys_node_num": node_num,
                            "label": label,
                            "component": component,
                            "value": value,
                        }
                    )

        # Convert the data to a pandas DataFrame
        df = pd.DataFrame(
            data, columns=["ansys_node_num", "label", "component", "value"]
        )
        return df
