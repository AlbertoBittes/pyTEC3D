import os
import sys
from abc import ABC, abstractmethod
import numpy as np
import struct
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp

from tree_gauge_helper import TreeGaugeHelper
from ansys_helper import AnsysHelper
from selector import Selector
from vtk_tables_helper import VTKTablesHelper
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
    APEC class to manage VTK unstructured grids and element types.
    """

    def __init__(self):
        """
        Initialize the APEC object without requiring a file path.
        """
        self.mesh = vtk.vtkUnstructuredGrid()

        self.jobname = None
        self.folderpath = None

        self.element_types = pd.DataFrame(
            columns=list(ELEMENT_TYPES_COLUMNS_TYPES.keys())
        )
        self.materials = pd.DataFrame(columns=list(MATERIALS_COLUMNS_TYPES.keys()))
        self.real_constants = pd.DataFrame(
            columns=list(REAL_CONSTANTS_COLUMNS_TYPES.keys())
        )
        self.time_parameters = DEFAULT_TIME_PARAMETERS
        self.apec_parameters = {}
        self.solver_settings = DEFAULT_SOLVER_SETTINGS
        self.transient_analysis = pd.DataFrame(
            columns=list(TRANSIENT_ANALYSIS_COLUMNS_TYPES.keys())
        )
        self.process = None
        self.num_processors = 1
        self.periodic_boundaries = []
        self.dof_couplings = []

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
    def read_mesh(self, vtk_file_path):
        """
        Reads a VTK unstructured grid from a file in XML format.

        Parameters:
            vtk_file_path (str): Path to the VTK file in XML format.
        """
        # Check if the mesh file exists
        if not os.path.exists(vtk_file_path):
            raise FileNotFoundError(f"Mesh file '{vtk_file_path}' not found.")

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file_path)
        reader.Update()
        self.mesh = reader.GetOutput()

        self.jobname = os.path.splitext(os.path.basename(vtk_file_path))[0]
        self.folderpath = os.path.dirname(vtk_file_path)

        # Prepare per-mesh folders inside the project directory
        mesh_root = os.path.join(self.folderpath, self.jobname)
        os.makedirs(mesh_root, exist_ok=True)

        self.load_steps_folder = os.path.join(mesh_root, "load_steps")
        os.makedirs(self.load_steps_folder, exist_ok=True)

        # Ensure relative paths resolve from the mesh directory
        os.chdir(self.folderpath)

        VTKGridHelper.initialize_vtk_arrays(self.mesh)

        Selector.ensure_selection_arrays(self.mesh)

        self.compute_and_store_cell_centroids()

        TreeGaugeHelper.add_midside_nodes_mask_array(self.mesh)

        print(f"VTK unstructured grid loaded from '{vtk_file_path}'.")
        # print a summary of the mesh
        print(f"Number of points: {self.mesh.GetNumberOfPoints()}")
        print(f"Number of cells: {self.mesh.GetNumberOfCells()}")
        print(
            f"Number of point data arrays: {self.mesh.GetPointData().GetNumberOfArrays()}"
        )
        print(
            f"Number of cell data arrays: {self.mesh.GetCellData().GetNumberOfArrays()}"
        )
        print(
            f"Number of field data arrays: {self.mesh.GetFieldData().GetNumberOfArrays()}"
        )

    def build_dof_couplings(self, print_summary: bool = False):
        """
        Construct and cache DOF coupling forests from CP_* point-data arrays.

        Parameters
        ----------
        print_summary : bool, optional
            When True, print a human-readable summary of the computed forests.

        Returns
        -------
        list[dict]
            The cached list of coupling forests, one per coupling set.
        """
        return self.refresh_dof_couplings(print_summary=print_summary)

    def refresh_dof_couplings(self, print_summary: bool = False):
        forests = DOFCouplingHelper.build_dof_couplings(self.mesh)
        self.dof_couplings = forests
        if print_summary:
            DOFCouplingHelper.print_dof_couplings(forests)
        return forests

    def save_mesh(self):
        """
        Saves the unstructured grid in self.mesh to a .vtu file named after self.jobname.
        Before writing, attaches the element_types DataFrame as field data to the grid.
        Validates that the first column of element_types has unique values, all values have a correspondence
        in the 'ansys_elem_num' array, and there are no undefined element types in the grid.

        Raises:
            ValueError: If self.jobname is not set, if the folder path is not set,
                        or if the validation fails.
        """
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
        output_file = os.path.join(self.folderpath, f"{self.jobname}.vtu")

        # Write the unstructured grid to the .vtu file
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(self.mesh)
        writer.SetDataModeToBinary()  # Set the writer to binary mode
        writer.Write()

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
    def get_cli(self):
        """
        Returns a CLIObject instance bound to this APEC_M instance.
        """
        return CLIObject(self)

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
        excel_file_path,
        data_to_read=list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys()),
    ):
        """
        Reads data from an Excel file and processes it based on the provided labels.

        Parameters:
            excel_file_path (str): The path to the Excel file.
            data_to_read (list): A list of data types to read and process.
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

    def read_excel_sheet(self, excel_file, sheet_name):
        """
        Processes an Excel sheet by determining the appropriate method to call
        (read_excel_df, read_excel_dict, or update_vtk_grid_from_excel) based on a global dictionary.

        Parameters:
            excel_file (str): The path to the Excel file.
            sheet_name (str): The name of the sheet to process.

        Raises:
            ValueError: If the sheet name is not recognized.
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

    def update_vtk_grid_from_excel(self, excel_file, sheet_name, data_type=None):
        """
        Delegates to ExcelIOHelper static method to update VTK grid from Excel.
        """
        ExcelIOHelper.update_vtk_grid_from_excel(
            self.mesh, excel_file, sheet_name, data_type
        )

    def read_excel_df(self, excel_file, sheet_name):
        """
        Wrapper for ExcelIOHelper.read_excel_df to allow getattr to work in read_excel_sheet.
        """
        return ExcelIOHelper.read_excel_df(excel_file, sheet_name)

    def read_excel_dict(self, excel_file, sheet_name):
        """
        Wrapper for ExcelIOHelper.read_excel_dict to allow getattr to work in read_excel_sheet.
        """
        return ExcelIOHelper.read_excel_dict(excel_file, sheet_name)

    def write_excel_cell_point_data(self, output_excel_file, data_type="cell_data"):
        """
        Delegates writing cell/point data to ExcelIOHelper static method.
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
    def write_to_binary(self, output_file, labels_to_write=None):
        """
        Writes the current load step data to a single binary file based on the provided labels.

        Parameters:
            output_file (str): The path to the output binary file.
            labels_to_write (list): A list of labels specifying what to write. Defaults to all supported labels.
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

    def write_load_step(self, load_step_number: int):
        if not self.jobname:
            raise ValueError("Job name is not set. Cannot write load step binaries.")

        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write load step binaries.")

        load_steps_folder = os.path.join(self.folderpath, self.jobname, "load_steps")
        os.makedirs(load_steps_folder, exist_ok=True)

        timestep_path = os.path.join(
            load_steps_folder, f"{self.jobname}_timestep_{load_step_number}.bin"
        )
        grid_data_path = os.path.join(
            load_steps_folder, f"{self.jobname}_grid_data_{load_step_number}.bin"
        )

        self.write_to_binary(timestep_path, labels_to_write=["time_parameters"])
        self.write_to_binary(
            grid_data_path,
            labels_to_write=["point_data", "cell_data"],
        )

    def write_model_data(self):
        if not self.jobname:
            raise ValueError("Job name is not set. Cannot write model data.")

        if not self.folderpath:
            raise ValueError("Folder path is not set. Cannot write model data.")

        os.makedirs(self.folderpath, exist_ok=True)

        output_path = os.path.join(self.folderpath, "model_data.bin")
        self.write_to_binary(
            output_path,
            labels_to_write=[
                "element_types",
                "materials",
                "real_constants",
                "periodic_boundaries",
                "dof_couplings",
            ],
        )

    def generate_transient_analysis(self):
        """
        Generates and writes load steps for transient analysis based on the rows in self.transient_analysis.

        For each row in the transient_analysis DataFrame:
            - Updates self.time_parameters with the values from the current row.
            - If a valid 'load_step_input_file' is provided, reads the corresponding Excel file and updates
              materials, real_constants, solver_settings, point_data, and cell_data using EXCEL_LOAD_STEP_LABELS.
            - Writes the current load step data to a binary file in the load_steps_folder.
            - Prints progress and warnings for missing input files.

        Side Effects:
            - Updates self.time_parameters for each load step.
            - Reads and processes Excel files for each load step if provided.
            - Writes binary files for each load step in the load_steps_folder.
            - Prints progress information to the console.
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

            self.write_load_step(idx + 1)

    def _write_label_to_binary(self, binary_file, label):
        """
        Writes a specific label to the binary file using the appropriate function.

        Parameters:
            binary_file (file object): An open binary file stream.
            label (str): The label to write (e.g., 'time_parameters', 'materials').

        Raises:
            ValueError: If the label is not recognized or the corresponding data is missing.
        """
        if label not in WRITING_FUNCTIONS_MAP:
            raise ValueError(f"Unrecognized label '{label}'.")

        # Get the function name from the map
        function_name = WRITING_FUNCTIONS_MAP[label]

        # Always call the function with (mesh or data, binary_file, label)
        if label in ["cell_data", "point_data"]:
            getattr(self, function_name)(self.mesh, binary_file, label)
        else:
            data = getattr(self, label, None)
            if data is None:
                raise ValueError(f"Data for '{label}' is missing.")
            getattr(self, function_name)(data, binary_file, label)

    # Wrapper methods for binary output functions to match WRITING_FUNCTIONS_MAP
    def _write_settings_to_binary(self, settings, binary_file, label):
        BinaryHelper.write_settings_to_binary(settings, binary_file, label)

    def _write_dataframe_to_binary(self, df, binary_file, label):
        BinaryHelper.write_dataframe_to_binary(df, binary_file, label)

    def _write_numpy_array_to_binary(self, mesh, binary_file, label):
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
        self, array_name, value, component=0, cell_or_point="point_data", mode="s"
    ):
        """
        Delegates selection to Selector static method.
        """
        Selector.select(self.mesh, array_name, value, component, cell_or_point, mode)

    def select_pos(self, axis, value, mode="s"):
        """
        Delegates position-based selection to Selector static method.
        """
        Selector.select_pos(self.mesh, axis, value, mode)

    def select_all(self, cell_or_point=None):
        """
        Delegates select_all to Selector static method.
        """
        Selector.select_all(self.mesh, cell_or_point)

    def unselect_all(self, cell_or_point=None):
        """
        Delegates unselect_all to Selector static method.
        """
        Selector.unselect_all(self.mesh, cell_or_point)

    def nsel_ext(self):
        """
        Delegates external node selection to Selector static method.
        """
        Selector.nsel_ext(self.mesh)

    def set_value_on_selection(
        self, array_name, value, component=0, cell_or_point="point_data"
    ):
        """
        Delegates set_value_on_selection to Selector static method.
        """
        Selector.set_value_on_selection(
            self.mesh, array_name, value, component, cell_or_point
        )

    def nsle(self, mode="s"):
        """Select nodes attached to currently selected elements."""
        Selector.nsle(self.mesh, mode=mode)

    def esln(self, all_nodes=True, mode="s"):
        """Select elements based on current node selection."""
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
    def plot(self):
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
    def read_bfelist(self, file_path):
        """
        Wrapper for AnsysHelpers.read_bfelist to allow getattr and consistent API.
        """
        AnsysHelper.read_bfelist(self.mesh, file_path)

    def convert_bfelist(self, bfelist_file, output_excel_file):
        """
        Wrapper for AnsysHelpers.convert_bfelist to allow getattr and consistent API.
        """
        AnsysHelper.convert_bfelist(self.mesh, bfelist_file, output_excel_file)

    def read_dlist(self, file_path):
        """
        Wrapper for AnsysHelpers.read_dlist to allow getattr and consistent API.
        """
        AnsysHelper.read_dlist(self.mesh, file_path)

    def convert_dlist(self, dlist_file, output_excel_file):
        """
        Wrapper for AnsysHelpers.convert_dlist to allow getattr and consistent API.
        """
        AnsysHelper.convert_dlist(self.mesh, dlist_file, output_excel_file)

    #
    #
    #
    #
    # General grid tools
    #
    #
    #
    def compute_and_store_cell_centroids(self):
        """
        Computes centroids for all cells and stores them in a VTK array 'cell_centroids' in CellData.
        """
        VTKGridHelper.store_cell_centroids(self.mesh)

    def compute_and_store_cell_volumes(self, array_name="cell_volume"):
        """
        Computes the volume (for 3D cells) or area (for 2D cells) for each cell in self.mesh,
        and stores the result as a cell data array named `array_name`.
        """
        VTKGridHelper.store_cell_volumes(self.mesh, array_name=array_name)
        print(
            f"Cell volumes/areas computed and stored in CellData array '{array_name}'."
        )

    def remove_array(self, array_name):
        """
        Removes an array from CellData or PointData using VTKGridHelpers.safe_remove_array.

        Parameters:
            array_name (str): The name of the array to remove.
        """
        VTKGridHelper.safe_remove_array(self.mesh, array_name)


    def add_boundary_array(self):
        """Add a new boundary_id_N array to the mesh and return its name."""
        return VTKGridHelper.add_boundary_array(self.mesh)
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
        self, array_name="Az", tree_name="TREE", cotree_name="COTREE"
    ):
        """
        Build the tree–cotree decomposition for midside edges in the mesh using TreeGaugeHelper.

        This method excludes all conducting regions (elements with KEYOPT1 not equal to 1, 2, or 3)
        by forcing all midside nodes of those elements into the cotree. It also analyzes the specified
        array in PointData (default "Az") to determine constrained edges (where the array value is not NaN).
        These are used to generate a set of constrained corner nodes, which are passed to the tree-cotree builder.

        After building the tree–cotree decomposition, all nodes with value 1 in the tree mask array
        (tree_name, default "TREE") will have their value in the selected array (array_name, default "Az") set to 0.

        Parameters:
            array_name (str): Name of the array in PointData to use for constrained edges (default "Az").
            tree_name (str): Name for the output tree mask array (default "TREE").
            cotree_name (str): Name for the output cotree mask array (default "COTREE").

        Raises:
            ValueError: If the specified array is not found in PointData.

        Side Effects:
            Adds/updates tree and cotree mask arrays in mesh PointData.
            Sets values in array_name to 0 for all nodes with value 1 in tree_name.
            Prints summary information about the decomposition.

        Returns:
            info (dict): Information about the tree/cotree decomposition (from TreeGaugeHelper).
        """
        arr = self.mesh.GetPointData().GetArray(array_name)
        if arr is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        arr_np = vtknp.vtk_to_numpy(arr)
        npts = self.mesh.GetNumberOfPoints()

        # Step 1: Find element type IDs with KEYOPT1 not in [1, 2, 3]
        invalid_keyopt_mask = ~self.element_types["KEYOPT1"].isin([0, 1, 2, 3])
        invalid_elemtype_ids = set(
            self.element_types.loc[invalid_keyopt_mask, "ansys_etype"]
        )

        # Step 2: Select cells with those element type IDs
        cell_data = self.mesh.GetCellData()
        ansys_elem_array = cell_data.GetArray("ansys_etype")
        if ansys_elem_array is None:
            raise ValueError("CellData array 'ansys_etype' not found.")
        ansys_elem_np = vtknp.vtk_to_numpy(ansys_elem_array)
        # Save current selection to restore later
        selection_cell_array = cell_data.GetArray("selection_cell")
        selection_cell_backup = vtknp.vtk_to_numpy(selection_cell_array).copy()
        # Select cells
        Selector.select(
            self.mesh,
            "ansys_etype",
            invalid_elemtype_ids,
            cell_or_point="cell_data",
            mode="s",
        )

        # Step 3: Extract selected cells and find midside nodes
        selected_cell_ids = np.where(
            vtknp.vtk_to_numpy(cell_data.GetArray("selection_cell")) == 1
        )[0]
        if len(selected_cell_ids) == 0:
            print("No cells with KEYOPT1 not in [0, 1,2,3] found.")
            forced_cotree_mask = np.zeros(npts, dtype=np.int32)
        else:
            # Extract selection
            id_array = vtknp.numpy_to_vtkIdTypeArray(selected_cell_ids, deep=True)
            selection_node = vtk.vtkSelectionNode()
            selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            selection_node.SetSelectionList(id_array)
            selection = vtk.vtkSelection()
            selection.AddNode(selection_node)
            extract_selection = vtk.vtkExtractSelection()
            extract_selection.SetInputData(0, self.mesh)
            extract_selection.SetInputData(1, selection)
            extract_selection.Update()
            selected_mesh = extract_selection.GetOutput()

            # Get MidsideNodes array from extracted mesh
            midside_arr = selected_mesh.GetPointData().GetArray("MidsideNodes")
            if midside_arr is None:
                raise ValueError("Extracted mesh missing 'MidsideNodes' array.")
            midside_np = vtknp.vtk_to_numpy(midside_arr)
            # Map extracted mesh point IDs back to original mesh point IDs
            orig_ids_arr = selected_mesh.GetPointData().GetArray("vtkOriginalPointIds")
            if orig_ids_arr is None:
                raise ValueError("Extracted mesh missing 'vtkOriginalPointIds' array.")
            orig_ids_np = vtknp.vtk_to_numpy(orig_ids_arr)
            # Build forced_cotree_mask: 1 at original midside node ids
            forced_cotree_mask = np.zeros(npts, dtype=np.int32)
            midside_orig_ids = orig_ids_np[midside_np == 1]
            forced_cotree_mask[midside_orig_ids] = 1

        # Restore previous selection
        vtknp.vtk_to_numpy(cell_data.GetArray("selection_cell"))[
            :
        ] = selection_cell_backup
        cell_data.GetArray("selection_cell").Modified()

        # Instead of adding constrained Az to forced_cotree_mask,
        # generate constrained corner node set from Az and pass to tree-cotree builder
        constrained_corners = TreeGaugeHelper.constrained_corners_from_midside_array(
            self.mesh, array_name
        )

        # constrained_corners_mask = np.zeros(npts, dtype=np.int32)
        # constrained_corners_mask[list(constrained_corners)] = 1
        # constrained_corners_arr = vtknp.numpy_to_vtk(
        #     constrained_corners_mask, deep=1, array_type=vtk.VTK_INT
        # )
        # constrained_corners_arr.SetName("ConstrainedCornerNodes")
        # self.mesh.GetPointData().AddArray(constrained_corners_arr)

        info = TreeGaugeHelper.build_tree_cotree_on_grid(
            self.mesh,
            constrained_point_ids=constrained_corners,
            forced_cotree_mask=forced_cotree_mask,
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
