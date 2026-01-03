import os
import sys
import copy
import struct
from os import PathLike
from typing import Any, BinaryIO, Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp

from tree_gauge_helper import TreeGaugeHelper
from ansys_helper import AnsysHelper
from selector import Selector
from vtk_grid_helper import VTKGridHelper
from binary_helper import BinaryHelper
from excel_io_helper import ExcelIOHelper
from cli_object import CLIObject
from constants import *
from dof_coupling_helper import DOFCouplingHelper


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
        self.materials: pd.DataFrame = pd.DataFrame(columns=list(MATERIALS_COLUMNS_TYPES.keys()))
        self.real_constants: pd.DataFrame = pd.DataFrame(
            columns=list(REAL_CONSTANTS_COLUMNS_TYPES.keys())
        )
        self.coupled_variables: pd.DataFrame = pd.DataFrame(
            columns=list(COUPLED_VARIABLES_COLUMNS_TYPES.keys())
        )
        self.time_parameters: Dict[str, Any] = copy.deepcopy(DEFAULT_TIME_PARAMETERS)
        self.apec_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_APEC_SETTINGS)
        self.solver_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self.transient_analysis: pd.DataFrame = pd.DataFrame(
            columns=list(TRANSIENT_ANALYSIS_COLUMNS_TYPES.keys())
        )
        self.periodic_boundaries: List[Any] = []
        self.dof_couplings: List[Dict[str, Any]] = []
        self.cli: CLIObject = CLIObject(self)

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
    def read_mesh(self, vtk_file_path: Union[str, PathLike[str]], reset_arrays: bool = True) -> None:
        """
        Reads a VTK unstructured grid from a file in XML format.

        Args:
            vtk_file_path: Path to the VTK file in XML format.
            reset_arrays: When True, reset initialized arrays to defaults after loading.
        """
        self.mesh = VTKGridHelper.read_mesh(vtk_file_path, reset_arrays=reset_arrays)

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

    def build_dof_couplings(self, print_summary: bool = False) -> List[Dict[str, Any]]:
        """
        Construct and cache DOF coupling forests from CP_* point-data arrays.

        Args:
            print_summary: When True, print a human-readable summary of the computed forests.

        Returns:
            List[Dict[str, Any]]: Cached coupling forests, one per CP_* id combination.
        """
        return self.refresh_dof_couplings(print_summary=print_summary)

    def refresh_dof_couplings(self, print_summary: bool = False) -> List[Dict[str, Any]]:
        """
        Recompute DOF coupling forests and optionally print a summary.

        Args:
            print_summary: When True, print a human-readable summary of the computed forests.

        Returns:
            List[Dict[str, Any]]: Updated list of coupling forests.
        """
        forests = DOFCouplingHelper.build_dof_couplings(
            self.mesh,
            coupled_variables=self.coupled_variables,
        )
        self.dof_couplings = forests
        if print_summary:
            DOFCouplingHelper.print_dof_couplings(forests)
        return forests

    def save_mesh(self, output_file: Optional[Union[str, PathLike[str]]] = None) -> None:
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
            VTKGridHelper.safe_remove_array(self.mesh, arr)

        # Define the output file name
        if output_file is None:
            output_file = os.path.join(self.folderpath, f"{self.jobname}.vtu")

        # Write the unstructured grid to the .vtu file
        VTKGridHelper.save_mesh(self.mesh, output_file, binary=True)

        print(f"Mesh saved to '{output_file}'.")

    #
    #
    #
    #
    #
    #  Methods for CLI interaction with solver (not used or tested yet).
    #
    #
    #
    #
    #
    def get_cli(self) -> CLIObject:
        """
        Return the CLI helper attached to this instance.
        """
        return self.cli

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

    def read_excel_sheet(self, excel_file: Union[str, PathLike[str]], sheet_name: str) -> None:
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
            func = getattr(self, processing_info)
            result = func(excel_file, sheet_name)
            if result is not None:
                setattr(self, sheet_name, result)
        elif isinstance(processing_info, tuple):
            method_name, data_type = processing_info
            func = getattr(self, method_name)
            func(excel_file, sheet_name, data_type)

    def update_vtk_grid_from_excel(
        self, excel_file: Union[str, PathLike[str]], sheet_name: str, data_type: Optional[str] = None
    ) -> None:
        """
        Load a sheet from Excel and write values into the mesh CellData/PointData.

        Mirrors `ExcelIOHelper.update_vtk_grid_from_excel`, reading `sheet_name` and
        writing into the VTK arrays specified there. If `data_type` is None, the sheet
        name is used to decide whether to update cell or point data.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Sheet containing the grid data.
            data_type: Target VTK data type ('elem' or 'node'); defaults to sheet name.

        Returns:
            None

        Raises:
            ValueError: If `data_type` is invalid or the sheet is missing required columns.
            FileNotFoundError: If the Excel file is missing.
        """
        ExcelIOHelper.update_vtk_grid_from_excel(
            self.mesh, excel_file, sheet_name, data_type
        )

    def read_excel_df(self, excel_file: Union[str, PathLike[str]], sheet_name: str) -> Optional[pd.DataFrame]:
        """
        Read a sheet into a DataFrame with column validation.

        Behaves like `ExcelIOHelper.read_excel_df`, ensuring columns match the expected
        schema defined in `DEFAULT_COLUMNS_ATTRIBUTES_TYPES` and coercing types.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Sheet name to parse.

        Returns:
            Optional[pd.DataFrame]: Parsed DataFrame or None if missing.

        Raises:
            FileNotFoundError: If the Excel file is missing.
            ValueError: If the sheet is not found or columns do not match expectations.
        """
        return ExcelIOHelper.read_excel_df(excel_file, sheet_name)

    def read_excel_dict(self, excel_file: Union[str, PathLike[str]], sheet_name: str) -> Optional[Dict[str, Any]]:
        """
        Read a sheet into a dictionary with key/type validation.

        Behaves like `ExcelIOHelper.read_excel_dict`, ensuring 'SETTING'/'VALUE' columns
        exist, keys match the expected schema, and values are coerced to the declared types.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Sheet name to parse.

        Returns:
            Optional[Dict[str, Any]]: Parsed dictionary or None if missing.

        Raises:
            FileNotFoundError: If the Excel file is missing.
            ValueError: If the sheet is missing, malformed, or contains invalid keys.
        """
        return ExcelIOHelper.read_excel_dict(excel_file, sheet_name)

    def read_periodic_boundaries(
        self, excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Read periodic boundary mappings from Excel.

        Delegates to ExcelIOHelper.read_periodic_boundaries to validate required
        columns and evaluate periodic matrices into NumPy arrays. On success, the
        parsed list is stored on self.periodic_boundaries.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Sheet name containing periodic boundary definitions.

        Returns:
            Optional[List[Dict[str, Any]]]: Parsed periodic boundaries or None if the sheet is missing.

        Raises:
            FileNotFoundError: If the Excel file is missing.
            ValueError: If the sheet is malformed or missing required values.
        """
        periodic_boundaries = ExcelIOHelper.read_periodic_boundaries(
            excel_file, sheet_name
        )
        if periodic_boundaries is not None:
            self.periodic_boundaries = periodic_boundaries
        return periodic_boundaries

    def write_excel_cell_point_data(
        self, output_excel_file: Union[str, PathLike[str]], data_type: str = "elem"
    ) -> None:
        """
        Export loads/boundary conditions to Excel, one sheet per data type.

        Mirrors `ExcelIOHelper.write_excel_cell_point_data`, flattening the structured
        load/bc array into a DataFrame with columns `vtk_id`, `label`, `component`, `value`
        and writing it to `output_excel_file` under the sheet `data_type`.

        Args:
            output_excel_file: Destination Excel path.
            data_type: 'elem' or 'node'.

        Returns:
            None

        Raises:
            ValueError: If `data_type` is invalid.
        """
        ExcelIOHelper.write_excel_cell_point_data(
            self.mesh, output_excel_file, data_type
        )

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
        self, output_file: Union[str, PathLike[str]], labels_to_write: Optional[Sequence[str]] = None
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

    def write_load_step(self, load_step_number: int, filename: Optional[Union[str, PathLike[str]]] = None) -> None:
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

    def write_time_step(self, load_step_number: int, filename: Optional[Union[str, PathLike[str]]] = None) -> None:
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
        if filename:
            timestep_path = os.path.join(self.folderpath, filename)
            target_dir = os.path.dirname(timestep_path) or self.folderpath
        else:
            if not self.jobname:
                raise ValueError("Job name is not set. Cannot write time step binary.")
            timestep_path = os.path.join(
                self.folderpath, f"{self.jobname}_timestep_{load_step_number}.bin"
            )
            target_dir = self.folderpath

        os.makedirs(target_dir, exist_ok=True)
        self.write_to_binary(
            timestep_path,
            labels_to_write=["time_parameters", "apec_settings"],
        )

    def write_loads_bcs(self, load_step_number: int = 1, filename: Optional[Union[str, PathLike[str]]] = None) -> None:
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
            labels_to_write=['node', "elem"],
        )

    def write_model_data(self, filename: Optional[Union[str, PathLike[str]]] = None) -> None:
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

        if filename:
            output_path = os.path.join(self.folderpath, filename)
            target_dir = os.path.dirname(output_path) or self.folderpath
            os.makedirs(target_dir, exist_ok=True)
        else:
            if not self.jobname:
                raise ValueError("Job name is not set. Cannot write model data.")

            os.makedirs(self.folderpath, exist_ok=True)
            output_path = os.path.join(self.folderpath, "model_data.bin")

        self.write_to_binary(
            output_path,
            labels_to_write=[
                "element_types",
                "materials",
                "real_constants",
                "solver_settings",
                "periodic_boundaries",
                "dof_couplings"
            ],
        )

    def generate_transient_analysis(self) -> None:
        """
        Generates and writes load steps for transient analysis based on the rows in self.transient_analysis.

        For each row in the transient_analysis DataFrame:
            - Updates self.time_parameters with the values from the current row.
            - If a valid 'load_step_input_file' is provided, reads the corresponding Excel file and updates
              materials, real_constants, solver_settings, point_data, and cell_data using EXCEL_LOAD_STEP_LABELS.
            - Writes the current load step data to a binary file inside the input folder.
            - Prints progress and warnings for missing input files.

        Side Effects:
            - Updates self.time_parameters for each load step.
            - Reads and processes Excel files for each load step if provided.
            - Writes binary files for each load step inside the input folder.
            - Prints progress information to the console.

        Returns:
            None
        """
        for idx, row in self.transient_analysis.iterrows():
            # Step 1: Update self.time_parameters with the values from the current row
            for column in DEFAULT_TIME_PARAMETERS.keys():
                if column in row:
                    self.time_parameters[column] = row[column]

            # Only read Excel if load_step_input_file is provided and not empty/NaN/"nan"
            input_file = row.get("load_step_input_file", None)
            if (
                input_file
                and pd.notna(input_file)
                and str(input_file).strip()
                and str(input_file).strip().lower() != "nan"
            ):
                self.read_excel(
                    input_file,
                    EXCEL_LOAD_STEP_LABELS,
                )
            else:
                print(
                    f"Skipping Excel read for load step {idx + 1}: no input file provided."
                )

            self.write_time_step(idx + 1)
            self.write_loads_bcs(idx + 1)

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

        # Get the function name from the map
        function_name = WRITING_FUNCTIONS_MAP[label]

        # Always call the function with (mesh or data, binary_file, label)
        if label in ["elem", 'node']:
            getattr(self, function_name)(self.mesh, binary_file, label)
        else:
            data = getattr(self, label, None)
            if data is None:
                raise ValueError(f"Data for '{label}' is missing.")
            getattr(self, function_name)(data, binary_file, label)

    # Wrapper methods for binary output functions to match WRITING_FUNCTIONS_MAP
    def _write_to_binary(self, settings: Any, binary_file: BinaryIO, label: str) -> None:
        """
        Wrapper dispatch to BinaryHelper.write_to_binary.

        Args:
            settings: Data payload to serialize.
            binary_file: Open binary stream.
            label: Label used for serialization.
        """
        BinaryHelper.write_to_binary(settings, binary_file, label)

    def _write_numpy_array_to_binary(self, mesh: vtk.vtkUnstructuredGrid, binary_file: BinaryIO, label: str) -> None:
        """
        Wrapper dispatch to BinaryHelper.write_numpy_array_to_binary.

        Args:
            mesh: Mesh providing point/cell data to serialize.
            binary_file: Open binary stream.
            label: 'elem' or 'node'.
        """
        BinaryHelper.write_numpy_array_to_binary(mesh, binary_file, label)

    #
    #
    #
    #
    # Selection methods
    #
    #
    #
    #
    #

    def select(
        self,
        array_name: str,
        value: Any,
        component: Union[int, str] = 0,
        type: str = 'node',
        mode: str = "s",
    ) -> None:
        """
        Select cells or points whose array values match the provided criteria.

        Mirrors `Selector.select` and supports multiple matching strategies:
            - Scalar equality for single values.
            - Inclusive range checks for tuples of length 2 (low, high).
            - Membership tests for sets.
            - Magnitude checks for vector arrays when `component='magnitude'`.

        The resulting boolean mask is written into `selection_cell` or `selection_point`
        (created if missing). How the new mask combines with the existing selection is
        controlled by `mode`:
            - 's' (set): replace the selection with the new mask.
            - 'a' (add): logical OR with the existing selection.
            - 'r' (reselect): logical AND with the existing selection.
            - 'u' (unselect): remove matches from the existing selection.

        Args:
            array_name: Data array used for filtering.
            value: Target value(s); scalar for equality, tuple(len=2) for range, set for membership.
            component: Component index or 'magnitude' for vector arrays.
            type: 'elem' or 'node' selection mask to update.
            mode: Selection mode — 's' replace, 'a' union, 'r' intersection, 'u' subtraction.

        Returns:
            None

        Raises:
            ValueError: If the target array is missing, mode is invalid, or a range tuple is malformed.
        """
        Selector.select(self.mesh, array_name, value, component, type, mode)

    def select_pos(self, axis: int, value: Any, mode: str = "s") -> None:
        """
        Select points based on their position along a given axis.

        Uses `Selector.select_pos` to update `selection_point` using the same value semantics
        and combination modes as `select`.

        Args:
            axis: Coordinate axis (0,1,2).
            value: Target value(s) (scalar, tuple for range, set for membership).
            mode: Selection mode ('s','a','r','u').

        Returns:
            None

        Raises:
            ValueError: If the mode is invalid or the range tuple length is not 2.
        """
        Selector.select_pos(self.mesh, axis, value, mode)

    def select_all(self, type: Optional[str] = None) -> None:
        """
        Mark all cells and/or points as selected.

        Args:
            type: 'elem', 'node', or None for both.

        Returns:
            None
        """
        Selector.select_all(self.mesh, type)

    def unselect_all(self, type: Optional[str] = None) -> None:
        """
        Clear selections for cells and/or points.

        Args:
            type: 'elem', 'node', or None for both.

        Returns:
            None
        """
        Selector.unselect_all(self.mesh, type)

    def nsel_ext(self) -> None:
        """
        Select external nodes of the currently selected cells.

        Sets `selection_point` for nodes on the exterior surface of selected cells using `Selector.nsel_ext`.
        """
        Selector.nsel_ext(self.mesh)

    def set_value_on_selection(
        self,
        array_name: str,
        value: Any,
        component: int = 0,
        type: str = 'node',
    ) -> None:
        """
        Assign a value to all selected cells or points in the given array.

        Args:
            array_name: Array to modify.
            value: Value to assign to selected entries.
            component: Component index for multi-component arrays.
            type: 'elem' or 'node'.

        Returns:
            None

        Raises:
            ValueError: If the target array is missing.
        """
        Selector.set_value_on_selection(
            self.mesh, array_name, value, component, type
        )

    def nsle(self, mode: str = "s") -> None:
        """
        Select nodes attached to the currently selected cells.

        Args:
            mode: Combination mode for updating `selection_point` ('s','a','r','u').

        Returns:
            None

        Raises:
            ValueError: If required selection arrays are missing.
        """
        Selector.nsle(self.mesh, mode=mode)

    def esln(self, all_nodes: bool = True, mode: str = "s") -> None:
        """
        Select cells based on the current node selection.

        Args:
            all_nodes: When True, select cells only if all their nodes are selected; otherwise if any node is selected.
            mode: Combination mode for updating `selection_cell` ('s','a','r','u').

        Returns:
            None

        Raises:
            ValueError: If required selection arrays are missing.
        """
        Selector.esln(self.mesh, all_nodes=all_nodes, mode=mode)

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
        RenameSource("Mesh", producer)
        Show(producer)
        Render()

    #
    #
    #
    #
    #
    #  ANSYS files processing functions
    #
    #
    #
    #
    #
    def read_bfelist(self, file_path: Union[str, PathLike[str]]) -> None:
        """
        Wrapper for AnsysHelpers.read_bfelist to allow getattr and consistent API.

        Args:
            file_path: Path to the BFELIST file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the BFELIST file does not exist.
            ValueError: If required identifiers are missing in the mesh CellData.
        """
        AnsysHelper.read_bfelist(self.mesh, file_path)

    def convert_bfelist(self, bfelist_file: Union[str, PathLike[str]], output_excel_file: Union[str, PathLike[str]]) -> None:
        """
        Wrapper for AnsysHelpers.convert_bfelist to allow getattr and consistent API.

        Args:
            bfelist_file: Path to the BFELIST file.
            output_excel_file: Destination Excel file path.

        Returns:
            None

        Raises:
            FileNotFoundError: If the BFELIST file does not exist.
            ValueError: If required identifiers are missing in the mesh CellData.
        """
        AnsysHelper.convert_bfelist(self.mesh, bfelist_file, output_excel_file)

    def read_dlist(self, file_path: Union[str, PathLike[str]]) -> None:
        """
        Wrapper for AnsysHelpers.read_dlist to allow getattr and consistent API.

        Args:
            file_path: Path to the DLIST file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the DLIST file does not exist.
            ValueError: If required identifiers are missing in the mesh PointData.
        """
        AnsysHelper.read_dlist(self.mesh, file_path)

    def convert_dlist(self, dlist_file: Union[str, PathLike[str]], output_excel_file: Union[str, PathLike[str]]) -> None:
        """
        Wrapper for AnsysHelpers.convert_dlist to allow getattr and consistent API.

        Args:
            dlist_file: Path to the DLIST file.
            output_excel_file: Destination Excel file path.

        Returns:
            None

        Raises:
            FileNotFoundError: If the DLIST file does not exist.
            ValueError: If required identifiers are missing in the mesh PointData.
        """
        AnsysHelper.convert_dlist(self.mesh, dlist_file, output_excel_file)

    def map_ansys_ids_to_vtk(
        self, data: "vtk.vtkDataSetAttributes", label: str
    ) -> np.ndarray:
        """
        Wrapper around AnsysHelper.map_data_to_vtk for convenience.

        Args:
            data: VTK CellData or PointData containing the ANSYS id array.
            label: Array name to map (e.g., 'ansys_elem_num', 'ansys_node_num').

        Returns:
            np.ndarray: 2D mapping array [vtk_id, ansys_id].
        """
        return AnsysHelper.map_data_to_vtk(data, label)

    #
    #
    #
    #
    # General grid tools
    #
    #
    #
    def compute_and_store_cell_centroids(self) -> None:
        """
        Computes centroids for all cells and stores them in a VTK array 'cell_centroids' in CellData.

        Returns:
            None
        """
        VTKGridHelper.store_cell_centroids(self.mesh)

    def compute_and_store_cell_volumes(self, array_name: str = "cell_volume") -> None:
        """
        Computes the volume (for 3D cells) or area (for 2D cells) for each cell in self.mesh,
        and stores the result as a cell data array named `array_name`.

        Args:
            array_name: Name for the output CellData array.

        Returns:
            None
        """
        VTKGridHelper.store_cell_volumes(self.mesh, array_name=array_name)
        print(
            f"Cell volumes/areas computed and stored in CellData array '{array_name}'."
        )

    def map_local_coordinates(
        self,
        coordinate_system: "CoordinateSystem",
        location: str = "point",
        selected_only: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Map mesh points or cell centroids into the provided coordinate system.

        Args:
            coordinate_system: Coordinate system used for local coordinates.
            location: "point" to map point IDs or "cell" to map cell IDs.
            selected_only: When True, use selection arrays to filter IDs.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (local_coordinates, ids)
        """
        return VTKGridHelper.map_local_coordinates(
            self.mesh,
            coordinate_system,
            location=location,
            selected_only=selected_only,
        )

    def rotate_vector_array_to_coordinate_system(
        self,
        array_name: str,
        coordinate_system: "CoordinateSystem",
        location: str = "point",
        to_local: bool = True,
    ) -> np.ndarray:
        """
        Rotate a 3-component array to/from the provided coordinate system.

        Args:
            array_name: Name of the array to rotate.
            coordinate_system: Coordinate system used for rotation.
            location: "point" or "cell".
            to_local: When True, rotate from global to local; otherwise local to global.

        Returns:
            np.ndarray: The rotated array (NumPy view).
        """
        return VTKGridHelper.rotate_vector_array_to_coordinate_system(
            self.mesh,
            array_name,
            coordinate_system,
            location=location,
            to_local=to_local,
        )

    def build_toroidal_projection_map(
        self,
        input_mesh: "vtk.vtkUnstructuredGrid",
        location: str = "point",
        tolerance: float = 1e-6,
        coordinate_system: Optional["CoordinateSystem"] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to self.mesh.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            location: "point" to map point IDs or "cell" to map cell IDs.
            tolerance: Absolute tolerance for matching radial distance and axial coordinate.
            coordinate_system: Cylindrical coordinate system used for local coordinates.

        Returns:
            np.ndarray: 2D array with input IDs in column 0 and mapped target IDs in the remaining columns.
                If a selection array exists on the input mesh for the chosen location,
                only selected IDs are included.
        """
        return VTKGridHelper.build_toroidal_projection_map(
            input_mesh,
            self.mesh,
            location=location,
            tolerance=tolerance,
            coordinate_system=coordinate_system,
        )

    def toroidal_projection(
        self,
        input_mesh: "vtk.vtkUnstructuredGrid",
        input_array_name: str,
        input_component: int,
        target_array_name: str,
        target_component: int = 0,
        projection_map: Optional[np.ndarray] = None,
        location: str = "point",
        tolerance: float = 1e-6,
        coordinate_system: Optional["CoordinateSystem"] = None,
    ) -> None:
        """
        Project a 2D scalar component onto self.mesh.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            input_array_name: Name of the input array on the input mesh.
            input_component: Component index to read from the input array.
            target_array_name: Name of the target array on self.mesh.
            projection_map: Optional precomputed toroidal projection map.
            location: "point" to map point data or "cell" to map cell data.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            coordinate_system: Cylindrical coordinate system used for local coordinates.
            target_component: Component index to write for scalar copy mode.

        Returns:
            None
        """
        VTKGridHelper.toroidal_projection(
            input_mesh,
            input_array_name,
            input_component,
            self.mesh,
            target_array_name,
            projection_map=projection_map,
            location=location,
            tolerance=tolerance,
            coordinate_system=coordinate_system,
            target_component=target_component,
        )

    def toroidal_projection_inverse(
        self,
        input_mesh: "vtk.vtkUnstructuredGrid",
        input_array_name: str,
        input_component: int,
        target_array_name: str,
        target_component: int = 0,
        projection_map: Optional[np.ndarray] = None,
        location: str = "point",
        tolerance: float = 1e-6,
        coordinate_system: Optional["CoordinateSystem"] = None,
    ) -> None:
        """
        Average scalar fields from a 3D input mesh onto self.mesh.

        Args:
            input_mesh: 3D unstructured grid providing the field.
            input_array_name: Name of the input array on the input mesh.
            input_component: Component index to read from the input array.
            target_array_name: Name of the target array on self.mesh.
            target_component: Component index to write for scalar targets.
            projection_map: Optional precomputed toroidal projection map.
            location: "point" to map point data or "cell" to map cell data.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            coordinate_system: Cylindrical coordinate system used for local coordinates.

        Returns:
            None
        """
        VTKGridHelper.toroidal_projection_inverse(
            input_mesh,
            input_array_name,
            input_component,
            self.mesh,
            target_array_name,
            target_component=target_component,
            projection_map=projection_map,
            location=location,
            tolerance=tolerance,
            coordinate_system=coordinate_system,
        )

    def remove_array(self, array_name: str) -> None:
        """
        Removes an array from CellData or PointData using VTKGridHelpers.safe_remove_array.

        Args:
            array_name: The name of the array to remove.

        Returns:
            None
        """
        VTKGridHelper.safe_remove_array(self.mesh, array_name)


    def add_boundary_array(self) -> str:
        """
        Add a new sequential `boundary_id_N` array to PointData and return its name.

        The array is initialized with zeros and appended after the highest existing
        `boundary_id_*` to avoid collisions.
        """
        return VTKGridHelper.add_boundary_array(self.mesh)

    def add_midside_nodes_mask_array(self, array_name: str = "MidsideNodes") -> "vtk.vtkIntArray":
        """
        Add a midside-node mask array to PointData (1 at midside nodes, 0 otherwise).

        Args:
            array_name: Name for the mask array (default "MidsideNodes").

        Returns:
            vtkIntArray: Added mask array.
        """
        return TreeGaugeHelper.add_midside_nodes_mask_array(self.mesh, array_name=array_name)

    def add_array_to_grid(
        self,
        name: str,
        location: str = "point",
        components: int = 1,
        fill: Any = 0,
        dtype: Optional[int] = vtk.VTK_DOUBLE,
        values: Optional[np.ndarray] = None,
    ) -> "vtk.vtkDataArray":
        """
        Wrapper to add an array to the mesh via VTKGridHelper for a consistent API.

        Args:
            name: Array name.
            location: 'cell' or 'point'.
            components: Number of components.
            fill: Default fill value if no values provided.
            dtype: VTK dtype for the array.
            values: Explicit values to assign.

        Returns:
            vtkDataArray: The created or existing array.
        """
        return VTKGridHelper.add_array_to_grid(
            self.mesh,
            name=name,
            location=location,
            components=components,
            fill=fill,
            dtype=dtype,
            values=values,
        )

    def initialize_arrays(self, arrays_to_init: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
        """
        Initialize required arrays on the mesh (CellData/PointData).

        Args:
            arrays_to_init: Custom array spec; defaults to constants.VTK_ARRAYS_TO_INIT when None.
        """
        VTKGridHelper.initialize_vtk_arrays(self.mesh, arrays_to_init)

    def reset_arrays(self, arrays_to_reset: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
        """
        Reset arrays on the mesh to their default fill values.

        Args:
            arrays_to_reset: Custom array spec; defaults to constants.VTK_ARRAYS_TO_INIT when None.
        """
        VTKGridHelper.reset_vtk_arrays_to_default(self.mesh, arrays_to_reset)

    def validate_dataframe_with_vtk_data(self, df: pd.DataFrame, array_name: str) -> None:
        """
        Validate a DataFrame's first column against a VTK array for uniqueness and coverage.

        Args:
            df: DataFrame whose first column will be validated.
            array_name: VTK array name to validate against (CellData).

        Returns:
            None

        Raises:
            ValueError: If validation fails or the VTK array is missing.
        """
        VTKGridHelper.validate_dataframe_with_vtk_data(
            self.mesh.GetCellData(), df, array_name
        )

    def add_component_from_selection(self, component_name: str, type: str = "cell") -> "vtk.vtkDataArray":
        """
        Create a component mask array from the current selection for cells/points.

        Args:
            component_name: Name for the component array to create.
            type: 'cell' to use CellData selection, 'point' for PointData selection.

        Returns:
            vtkDataArray: The created component mask array.

        Raises:
            ValueError: If selection arrays are missing or type is invalid.
        """
        return VTKGridHelper.add_component(self.mesh, component_name, type)

    def update_vtk_grid_from_dataframe(self, vtk_data: "vtk.vtkDataSetAttributes", df: pd.DataFrame) -> None:
        """
        Wrapper for VTKGridHelper.update_vtk_grid to write DataFrame values into arrays.

        Args:
            vtk_data: Target VTK DataSetAttributes (CellData or PointData).
            df: DataFrame with columns 'vtk_id', 'label', 'component', 'value'.

        Returns:
            None

        Raises:
            ValueError: If any referenced VTK array is missing.
        """
        VTKGridHelper.update_vtk_grid(vtk_data, df)
    #
    #
    #
    #
    #
    #
    # caller for tree-cotree gauge helper
    #
    #
    #
    #
    #
    def build_tree_cotree_gauge(
        self,
        array_name: str = "Az",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Dict[str, int]:
        """
        Build the tree–cotree decomposition for midside edges in the mesh using TreeGaugeHelper.
        It analyzes the specified array in PointData (default "Az") to determine constrained edges
        (where the array value is not NaN). These are used to generate a set of constrained
        corner nodes, which are passed to the tree-cotree builder.

        After building the tree–cotree decomposition, all nodes with value 1 in the tree mask array
        (tree_name, default "TREE") will have their value in the selected array (array_name, default "Az") set to 0.

        Args:
            array_name: Name of the PointData array to use for constrained edges (default "Az").
            tree_name: Name for the output tree mask array (default "TREE").
            cotree_name: Name for the output cotree mask array (default "COTREE").

        Returns:
            Dict[str, int]: Information about the tree/cotree decomposition (from TreeGaugeHelper).

        Raises:
            ValueError: If required array (array_name) is missing in PointData.

        Side Effects:
            Adds/updates tree and cotree mask arrays in mesh PointData.
            Sets values in array_name to 0 for all nodes with value 1 in tree_name.
            Prints summary information about the decomposition.
        """
        arr = self.mesh.GetPointData().GetArray(array_name)
        if arr is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        arr_np = vtknp.vtk_to_numpy(arr)
        npts = self.mesh.GetNumberOfPoints()

        # Instead of adding constrained Az to forced_cotree_mask,
        # generate constrained corner node set from Az and pass to tree-cotree builder
        constrained_corners = TreeGaugeHelper.constrained_corners_from_midside_array(
            self.mesh, array_name
        )

        info = TreeGaugeHelper.build_tree_cotree_on_grid(
            self.mesh,
            constrained_point_ids=constrained_corners,
            tree_name=tree_name,
            cotree_name=cotree_name,
        )
        print(f"Tree–cotree gauge built. Info: {info}")

        # Set value in array_name to 0 for all nodes with value 1 in tree_name
        tree_arr = self.mesh.GetPointData().GetArray(tree_name)
        if tree_arr is not None:
            tree_np = vtknp.vtk_to_numpy(tree_arr)
            arr_np[tree_np == 1] = 0
            arr.Modified()

        return info


# Only execute this block if the script is run directly
if __name__ == "__main__":
    # Example usage
    apec = APEC_M()
    # print("APEC instance created.")
