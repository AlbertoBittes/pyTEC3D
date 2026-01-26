import os
import sys
import copy
import struct
from os import PathLike
from typing import (
    Any,
    BinaryIO,
    Dict,
    Optional,
    Sequence,
    Union,
)

import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp

from tree_gauge import TreeGauge
from model_edit import ModelEdit
from binary import Binary
from excel import Excel
from solver import Solver
from constants import *


# Add the current directory to the Python path to ensure imports work
try:
    import openpyxl
except ImportError:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    local_packages_path = os.path.join(current_dir, "local_packages")
    sys.path.append(local_packages_path)
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    try:
        import openpyxl

        # print("openpyxl imported successfully.")
    except ImportError:
        print(
            "Failed to import openpyxl. Ensure the package is downloaded and the path is correct."
        )


class APEC_M:
    """
    APEC class to manage VTK unstructured grids, metadata, and IO utilities.
    """

    def __init__(self) -> None:
        """
        Initialize the APEC object without requiring a file path.
        """
        self.mesh: vtk.vtkUnstructuredGrid = vtk.vtkUnstructuredGrid()

        self.jobname: Optional[str] = None
        self.folderpath: Optional[str] = None

        self.element_types: pd.DataFrame = pd.DataFrame(
            columns=list(ELEMENT_TYPES_COLUMNS_TYPES.keys())
        )
        self.materials: pd.DataFrame = pd.DataFrame(
            columns=list(MATERIALS_COLUMNS_TYPES.keys())
        )
        self.real_constants: pd.DataFrame = pd.DataFrame(
            columns=list(REAL_CONSTANTS_COLUMNS_TYPES.keys())
        )
        self.time_parameters: Dict[str, Any] = copy.deepcopy(DEFAULT_TIME_PARAMETERS)
        self.apec_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_APEC_SETTINGS)
        self.solver_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self.dof_couplings: pd.DataFrame = pd.DataFrame(
            columns=list(DOF_COUPLINGS_COLUMNS_TYPES.keys())
        )
        self.solver: Solver = Solver(self)

    #
    #
    #
    #
    #
    #  Basic methods: read, initialize and save model
    #
    #
    #
    #
    #
    def read_mesh(
        self, vtk_file_path: Union[str, PathLike[str]], reset_arrays: bool = True
    ) -> None:
        """
        Reads a VTK unstructured grid from a file in XML format.

        Args:
            vtk_file_path: Path to the VTK file in XML format.
            reset_arrays: When True, reset initialized arrays to defaults after loading.
        """
        self.mesh = ModelEdit.read_mesh(vtk_file_path, reset_arrays=reset_arrays)

        mesh_path = os.fspath(vtk_file_path)
        self.jobname = os.path.splitext(os.path.basename(mesh_path))[0]
        self.folderpath = os.path.dirname(mesh_path)

        # Prepare per-mesh folders inside the project directory
        mesh_root = os.path.join(self.folderpath, self.jobname)
        input_folder = os.path.join(mesh_root, "input")
        os.makedirs(input_folder, exist_ok=True)

        self.load_steps_folder = input_folder

        # Ensure relative paths resolve from the mesh directory
        os.chdir(self.folderpath)

    def save_mesh(
        self, output_file: Optional[Union[str, PathLike[str]]] = None
    ) -> None:
        """
        Saves the unstructured grid in self.mesh to a .vtu file.
        Before writing, attaches the element_types DataFrame as field data to the grid.
        Validates that the first column of element_types has unique values, all values have a correspondence
        in the 'ansys_elem_num' array, and there are no undefined element types in the grid.

        Returns:
            None

        Raises:
            ValueError: If output_file is not provided and self.jobname or folderpath is not set.
        """
        if output_file is None:
            if not self.jobname:
                raise ValueError(
                    "Job name is not set. Please set self.jobname before saving the mesh."
                )

            if not self.folderpath:
                raise ValueError(
                    "Folder path is not set. Please set self.folderpath before saving the mesh."
                )

        # Remove arrays listed in ARRAYS_TO_REMOVE using safe_remove_array
        for arr in ARRAYS_TO_REMOVE:
            ModelEdit.remove_array(self.mesh, arr)

        # Define the output file name
        if output_file is None:
            output_file = os.path.join(self.folderpath, f"{self.jobname}.vtu")

        # Write the unstructured grid to the .vtu file
        ModelEdit.save_mesh(self.mesh, output_file, binary=True)

        print(f"Mesh saved to '{output_file}'.")

    #
    #
    #
    #
    #
    #  Excel I/O methods
    #
    #
    #
    #
    #
    def read_excel(
        self,
        excel_file_path: Union[str, PathLike[str]],
        data_to_read: Sequence[str] = tuple(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys()),
    ) -> None:
        """
        Reads data from an Excel file and processes it based on the provided labels.

        Args:
            excel_file_path: Path to the Excel file.
            data_to_read: Iterable of sheet labels to read and process.

        Returns:
            None

        Raises:
            FileNotFoundError: If the Excel file is missing.
        """
        # Check if the Excel file exists
        if not os.path.exists(excel_file_path):
            raise FileNotFoundError(f"Excel file '{excel_file_path}' not found.")

        # Loop through the data_to_read labels and process each one
        for label in data_to_read:
            try:
                self.read_excel_sheet(excel_file_path, label)
            except ValueError as e:
                print(f"Warning: {e}")

    def read_excel_sheet(
        self, excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> None:
        """
        Processes an Excel sheet by determining the appropriate method to call
        (read_excel_df, read_excel_dict, or update_vtk_grid_from_excel) based on a global dictionary.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to process.

        Raises:
            ValueError: If the sheet name is not recognized or mapped method is missing.
        """
        if sheet_name not in SHEET_PROCESSING_MAP:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(SHEET_PROCESSING_MAP.keys())}."
            )
        processing_info = SHEET_PROCESSING_MAP[sheet_name]
        print(f"Reading sheet '{sheet_name}' from '{excel_file}'...")
        if isinstance(processing_info, str):
            func = getattr(Excel, processing_info)
            result = func(excel_file, sheet_name)
            if result is not None:
                setattr(self, sheet_name, result)
        elif isinstance(processing_info, tuple):
            method_name, data_type = processing_info
            func = getattr(Excel, method_name)
            func(self.mesh, excel_file, sheet_name, data_type)

    #
    #
    #
    #
    #
    #  Binary Output methods
    #
    #
    #
    #
    #
    def write_to_binary(
        self,
        output_file: Union[str, PathLike[str]],
        labels_to_write: Optional[Sequence[str]] = None,
    ) -> None:
        """
        Writes the current load step data to a single binary file based on the provided labels.

        Args:
            output_file: Path to the output binary file.
            labels_to_write: Labels specifying what to write. Defaults to all supported labels.

        Returns:
            None

        Raises:
            ValueError: If an unrecognized label is provided.
        """

        # Default labels to write if none are provided
        if labels_to_write is None:
            labels_to_write = list(WRITING_LABELS)

        # Open the file in write mode to overwrite if it already exists
        with open(output_file, "wb") as binary_file:
            binary_file.write(struct.pack("I", len(labels_to_write)))

            for label in labels_to_write:
                try:
                    if label not in WRITING_LABEL_TO_ID:
                        raise ValueError(f"Unrecognized label '{label}'.")

                    label_id = WRITING_LABEL_TO_ID[label]
                    binary_file.write(struct.pack("I", label_id))
                    self._write_label_to_binary(binary_file, label)
                    print(f"{label.replace('_', ' ').capitalize()} written.")
                except ValueError as e:
                    print(f"Warning: {e}")

        print(f"Binary written to '{output_file}'.")

    def write_load_step(
        self,
        load_step_number: int,
        filename: Optional[Union[str, PathLike[str]]] = None,
    ) -> None:
        """
        Write both time and grid data for a specific load step.

        Args:
            load_step_number: Load step index (1-based).
            filename: Optional custom filename for the grid binary; time binary uses default naming.

        Returns:
            None

        Raises:
            ValueError: If folderpath is unset.
        """

        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write load step binaries.")

        # Backward compatibility: writing both time and grid data separately.
        self.write_time_step(load_step_number)
        self.write_loads_bcs(load_step_number, filename=filename)

    def write_time_step(
        self,
        load_step_number: int,
        filename: Optional[Union[str, PathLike[str]]] = None,
    ) -> None:
        """
        Writes a time step binary containing time parameters and APEC settings
        in the same folder as the mesh file.

        Args:
            load_step_number: Load step index (1-based).
            filename: Optional override for the output filename.

        Returns:
            None

        Raises:
            ValueError: If folderpath or jobname (when filename absent) is unset.
        """
        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write time step binary.")
        time_folder = os.path.join(self.folderpath, "_time_")
        if filename:
            timestep_path = os.path.join(time_folder, filename)
            target_dir = os.path.dirname(timestep_path) or time_folder
        else:
            if not self.jobname:
                raise ValueError("Job name is not set. Cannot write time step binary.")
            timestep_path = os.path.join(
                time_folder, f"{self.jobname}_timestep_{load_step_number}.bin"
            )
            target_dir = time_folder

        os.makedirs(target_dir, exist_ok=True)
        self.write_to_binary(
            timestep_path,
            labels_to_write=["time_parameters", "apec_settings", "dof_couplings"],
        )

    def write_loads_bcs(
        self,
        load_step_number: int = 1,
        filename: Optional[Union[str, PathLike[str]]] = None,
    ) -> None:
        """
        Writes a grid binary containing point and cell data in the mesh input folder.

        Args:
            load_step_number: Load step index (1-based). Default is 1.
            filename: Optional override for the output filename.

        Returns:
            None

        Raises:
            ValueError: If folderpath or jobname (when filename absent) is unset.
        """
        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write grid step binary.")

        if filename:
            grid_path = os.path.join(self.folderpath, filename)
            target_dir = os.path.dirname(grid_path) or self.folderpath
        else:
            if not self.jobname:
                raise ValueError("Job name is not set. Cannot write grid step binary.")
            load_steps_folder = os.path.join(self.folderpath, self.jobname, "input")
            os.makedirs(load_steps_folder, exist_ok=True)
            grid_path = os.path.join(
                load_steps_folder, f"{self.jobname}_grid_{load_step_number}.bin"
            )
            target_dir = load_steps_folder

        os.makedirs(target_dir, exist_ok=True)
        self.write_to_binary(
            grid_path,
            labels_to_write=["node", "elem"],
        )

    def write_model_data(
        self, filename: Optional[Union[str, PathLike[str]]] = None
    ) -> None:
        """
        Write model-wide metadata (element types, materials, etc.) to a binary file.

        Args:
            filename: Optional custom output filename; defaults to 'model_data.bin' in folderpath.

        Returns:
            None

        Raises:
            ValueError: If folderpath or jobname (when filename absent) is unset.
        """

        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write model data.")

        model_data_folder = os.path.join(self.folderpath, "_model_data_")
        if filename:
            output_path = os.path.join(model_data_folder, filename)
            target_dir = os.path.dirname(output_path) or model_data_folder
            os.makedirs(target_dir, exist_ok=True)
        else:
            if not self.jobname:
                raise ValueError("Job name is not set. Cannot write model data.")

            os.makedirs(model_data_folder, exist_ok=True)
            output_path = os.path.join(model_data_folder, "model_data.bin")

        self.write_to_binary(
            output_path,
            labels_to_write=[
                "element_types",
                "materials",
                "real_constants",
                "solver_settings",
            ],
        )

    def _write_label_to_binary(self, binary_file: BinaryIO, label: str) -> None:
        """
        Writes a specific label to the binary file using the appropriate function.

        Args:
            binary_file: An open binary file stream.
            label: The label to write (e.g., 'time_parameters', 'materials').

        Raises:
            ValueError: If the label is not recognized or the corresponding data is missing.
        """
        if label not in WRITING_FUNCTIONS_MAP:
            raise ValueError(f"Unrecognized label '{label}'.")

        # Always call the function with (mesh or data, binary_file, label)
        if label in ["elem", "node"]:
            Binary.write_numpy_array_to_binary(self.mesh, binary_file, label)
        else:
            data = getattr(self, label, None)
            if data is None:
                raise ValueError(f"Data for '{label}' is missing.")
            Binary.write_to_binary(data, binary_file, label)

    #
    #
    #
    #
    #
    #  Method to plot the model in ParaView
    #
    #
    #
    #
    #
    def plot(self) -> None:
        """
        Replot the mesh in the active view in Paraview.
        """
        from paraview.simple import TrivialProducer, RenameSource, Show, Render

        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        print("Plotting the mesh...")
        producer = TrivialProducer()
        producer.GetClientSideObject().SetOutput(self.mesh)
        RenameSource(self.jobname, producer)
        Show(producer)
        Render()


# Only execute this block if the script is run directly
if __name__ == "__main__":
    # Example usage
    apec = APEC_M()
    # print("APEC instance created.")
