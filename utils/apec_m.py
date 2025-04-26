import pandas as pd
from abc import ABC, abstractmethod
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import struct
import os
import re


class APEC:
    """
    APEC class to manage VTK unstructured grids and element types.
    """

    def __init__(self):
        """
        Initialize the APEC object without requiring a file path.
        """
        # self.csv_file_path = None
        self.element_types_collector = None
        self.mesh = vtk.vtkUnstructuredGrid()

        self.initialize_dataframes()
        self.jobname = None
        self.folderpath = None
        self.solve_settings = {}
        self.solve_settings["time"] = 1
        self.solve_settings["nsubsteps"] = 1
        self.solve_settings["timint"] = False
        self.solve_settings["kbc"] = True
        self.solve_settings["outres"] = True
        self.solve_settings["solver"] = "Direct"
        self.solve_settings["preconditioner"] = "LU"
        self.solve_settings["theta"] = 1.0

    def initialize_dataframes(self):
        """
        Initialize the DataFrames with the proper columns.
        """
        self.element_types = pd.DataFrame(
            columns=[
                "ansys_etype",
                "element_type",
                "keyopt1",
                "keyopt2",
                "keyopt3",
                "keyopt4",
                "keyopt5",
                "keyopt6",
                "keyopt7",
                "keyopt8",
                "keyopt9",
            ]
        )
        self.materials = pd.DataFrame(
            columns=[
                "ansys_material_type",
                "murx",
                "mury",
                "murz",
                "rsvx",
                "rsvy",
                "rsvz",
                "perx",
                "pery",
                "perz",
            ]
        )
        self.time_data = pd.DataFrame(
            columns=[
                "time_step_num",
                "time",
                "number_of_substeps",
                "transient_effects",
                "ramped_loads",
                "write_substep_results",
            ]
        )
        self.apec_parameters = pd.DataFrame(columns=["Setting", "Value"])

        # print("DataFrames initialized with proper columns.")

    def read_dataframe_from_csv(self, csv_file_path):
        """
        Generic function to read a DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.

        Returns:
            pd.DataFrame: The DataFrame read from the CSV file.
        """
        if not os.path.exists(csv_file_path):
            raise FileNotFoundError(f"CSV file '{csv_file_path}' not found.")
        return pd.read_csv(csv_file_path)

    def read_element_types(self, csv_file_path):
        """
        Read the element types DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.element_types = self.read_dataframe_from_csv(csv_file_path)
        print(f"Element types loaded from '{csv_file_path}'.")

    def read_materials(self, csv_file_path):
        """
        Read the materials DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.materials = self.read_dataframe_from_csv(csv_file_path)
        print(f"Materials loaded from '{csv_file_path}'.")

    def read_time_data(self, csv_file_path):
        """
        Read the time data DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.time_data = self.read_dataframe_from_csv(csv_file_path)
        print(f"Time data loaded from '{csv_file_path}'.")

    def read_apec_settings(self, csv_file_path):
        """
        Read the APEC settings DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.apec_settings = self.read_dataframe_from_csv(csv_file_path)
        print(f"APEC settings loaded from '{csv_file_path}'.")

    # def read_element_types(self, csv_file_path):
    #     """
    #     Load element types from the CSV file into the ElementTypesCollector.

    #     Parameters:
    #         csv_file_path (str): Path to the CSV file for initializing the ElementTypesCollector.
    #     """
    #     self.csv_file_path = csv_file_path
    #     self.element_types_collector = ElementTypesCollector(self.csv_file_path)

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

        # Set the current working directory
        os.chdir(self.folderpath)

        # Create a folder named 'loads_bcs' in the current working directory
        self.loads_bcs_folder = os.path.join(os.getcwd(), "loads_bcs")
        if not os.path.exists(self.loads_bcs_folder):
            os.makedirs(self.loads_bcs_folder)

        self.distribute_dofs()

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

    # def read_pipeline(self):
    #     # Get the active object in the pipeline
    #     self.mesh = GetActiveSource()

    #     # Check if an active object exists
    #     if self.pipeline_object is None:
    #         raise RuntimeError("No active object in the pipeline. Please select an object.")

    def write_loads_to_binary(self, array_name, time_step=None):
        """
        Writes a binary file for loads from a VTK cell data array.
        Handles multi-component arrays by iterating over components and writing the component number.
        Skips NaN values using NumPy masking for optimization.
        """
        # Get the cell data array
        cell_data = self.mesh.GetCellData()
        load_array = cell_data.GetArray(array_name)
        if load_array is None:
            raise ValueError(f"Array '{array_name}' not found in the VTK file.")

        # Convert the VTK array to a NumPy array
        load_array_np = vtknp.vtk_to_numpy(load_array)

        # Get the number of components and cells
        num_components = load_array.GetNumberOfComponents()
        num_cells = self.mesh.GetNumberOfCells()

        if time_step is None:
            output_file = os.path.join(
                self.loads_bcs_folder, f"{self.jobname}_loads.bin"
            )
        else:
            output_file = os.path.join(
                self.loads_bcs_folder, f"{self.jobname}_loads_{time_step}.bin"
            )

        # Open the output file in binary mode
        with open(output_file, "wb") as file:
            # # Write the length of the load name
            # load_name_length = len(array_name)
            # file.write(struct.pack("Q", load_name_length))  # 'Q' is for unsigned long long (size_t)

            # # Write the load name
            # file.write(array_name.encode("utf-8"))

            # # Write the number of components
            # file.write(struct.pack("I", num_components))  # 'I' is for unsigned int

            # Iterate over components
            for component in range(num_components):
                print(f"Processing component {component} of array '{array_name}'")

                # Extract the values for the current component
                component_values = load_array_np[:, component]
                mask = ~np.isnan(component_values)

                # Get the filtered cell IDs and values
                filtered_cell_ids = np.arange(num_cells)[mask]
                filtered_values = component_values[mask]

                # Write the data for this component
                for cell_id, value in zip(filtered_cell_ids, filtered_values):
                    if np.isnan(value):
                        continue
                    # Write the element number (cell ID)
                    file.write(
                        struct.pack("I", cell_id)
                    )  # 'I' is for uint32_t (dof_id_type)

                    # Write the component number
                    file.write(struct.pack("I", component))  # 'I' is for unsigned int

                    # Write the load value
                    file.write(struct.pack("d", value))  # 'd' is for double (Real)

    def write_boundary_conditions_to_binary(self, array_name, time_step=None):
        """
        Writes a binary file for boundary conditions from a VTK point data array.
        Handles multi-component arrays by appending "_x", "_y", "_z", etc., to the variable name.
        Skips NaN values using NumPy masking for optimization.
        """
        # Get the point data array
        point_data = self.mesh.GetPointData()
        bc_array = point_data.GetArray(array_name)
        if bc_array is None:
            raise ValueError(f"Array '{array_name}' not found in the VTK file.")

        # Convert the VTK array to a NumPy array
        bc_array_np = vtknp.vtk_to_numpy(bc_array)

        # Get the number of components and points
        num_components = bc_array.GetNumberOfComponents()
        num_points = self.mesh.GetNumberOfPoints()

        if time_step is None:
            output_file = os.path.join(self.loads_bcs_folder, f"{self.jobname}_bcs.bin")
        else:
            output_file = os.path.join(
                self.loads_bcs_folder, f"{self.jobname}_bcs_{time_step}.bin"
            )

        # Open the output file in binary mode
        with open(output_file, "wb") as file:
            # Iterate over components
            for component in range(num_components):
                # Generate the label for this component
                if num_components == 1:
                    # If there is only one component, use the array name directly
                    component_label = array_name
                else:
                    # Append '_x', '_y', '_z', etc., for multi-component arrays
                    component_label = (
                        f"{array_name}_{chr(120 + component)}"  # 'x', 'y', 'z', etc.
                    )
                # Create a mask to filter out NaN values for this component
                component_values = bc_array_np[:, component]
                mask = ~np.isnan(component_values)

                # Get the filtered point IDs and values
                filtered_point_ids = np.arange(num_points)[mask]
                filtered_values = component_values[mask]

                # Write the data for this component
                for point_id, value in zip(filtered_point_ids, filtered_values):
                    if np.isnan(value):
                        continue
                    # Write the point ID (node number)
                    file.write(
                        struct.pack("I", point_id)
                    )  # 'I' is for uint32_t (dof_id_type)

                    # Write the label length
                    label_length = len(component_label)
                    file.write(struct.pack("B", label_length))  # 'B' is for uint8_t

                    # Write the label
                    file.write(component_label.encode("utf-8"))

                    # Write the value
                    file.write(struct.pack("d", value))  # 'd' is for double (Real)

    def distribute_dofs(self):
        """
        Initialize the object's vtkUnstructuredGrid by adding variables as point data
        and element loads as cell data.
        """

        self.element_types_collector = ElementTypesCollector(self)

        element_types = self.element_types_collector.collect_type_classes()

        # Add point data arrays
        for element_id, element in element_types.items():
            variables = element.getVariables()
            for var_name, num_components in variables.items():
                array = vtk.vtkDoubleArray()
                array.SetName(var_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(self.mesh.GetNumberOfPoints())

                # Initialize all components to NaN
                for component in range(num_components):
                    array.FillComponent(component, np.nan)

                self.mesh.GetPointData().AddArray(array)

        # Add cell data arrays
        for element_id, element in element_types.items():
            element_loads = element.getElementLoads()
            for load_name, num_components in element_loads.items():
                array = vtk.vtkDoubleArray()
                array.SetName(load_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(self.mesh.GetNumberOfCells())

                # Initialize all components to NaN
                for component in range(num_components):
                    array.FillComponent(component, np.nan)

                self.mesh.GetCellData().AddArray(array)

    def plot(self):
        """
        Replot the mesh in the active view in Paraview.
        """
        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        print("Plotting the mesh...")
        producer = TrivialProducer()
        producer.GetClientSideObject().SetOutput(self.mesh)
        RenameSource("Mesh", producer)
        # producer.GetExecutive().GetOutputInformation(0).Set(
        #     vtk.vtkDataObject.DATA_PIECE_NAME(), "Mesh"
        # )
        # for table, name in [
        #     (self.element_types, "Element Types"),
        #     (self.materials, "Materials"),
        #     (self.time_data, "Time Data"),
        #     (self.apec_settings, "APEC Settings"),
        # ]:
        #     # print(f"Adding table to pipeline: {name}")
        #     table_producer = TrivialProducer()
        #     table_producer.GetClientSideObject().SetOutput(table)
        #     RenameSource(name, table_producer)
        #     Show(table_producer)

        Show(producer)
        Render()

    def parse_element_bfelist_data(self, file_path):
        """
        Parses element load data from the given file and associates a label based on the file's content.
        If the first row contains "CURRENT DENSITIES", the label is set to "Js".

        Args:
            file_path (str): Path to the file containing element load data.

        Returns:
            tuple: A tuple containing:
                - load_label (str): The label associated with the load (e.g., "Js").
                - loads_dict (dict): A dictionary where keys are element IDs and values are dictionaries
                of component values (0, 1, 2), defaulting to np.nan if not present.
        """
        import re
        import numpy as np

        # Regular expressions to match element and density lines
        element_pattern = re.compile(r"ELEMENT=\s+(\d+)")
        density_pattern = re.compile(
            r"([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)\s+([-+]?\d*\.\d+|\d+)"
        )

        # Initialize variables
        current_densities = {}
        current_element = None
        load_label = None

        with open(file_path, "r") as file:
            # Check the first line for "CURRENT DENSITIES"
            first_line = file.readline()
            if "CURRENT DENSITIES" in first_line:
                load_label = "Js"

            # Parse the rest of the file
            for line in file:
                # Match element line
                element_match = element_pattern.match(line)
                if element_match:
                    current_element = int(element_match.group(1))
                    current_densities[current_element] = {
                        0: np.nan,  # Default to NaN for Component 0
                        1: np.nan,  # Default to NaN for Component 1
                        2: np.nan,  # Default to NaN for Component 2
                    }
                    continue

                # Match density line
                if current_element is not None:
                    density_match = density_pattern.match(line)
                    if density_match:
                        current_densities[current_element][0] = float(
                            density_match.group(1)
                        )  # Component 0
                        current_densities[current_element][1] = float(
                            density_match.group(2)
                        )  # Component 1
                        current_densities[current_element][2] = float(
                            density_match.group(3)
                        )  # Component 2

        return load_label, current_densities

    def update_vtk_grid_bfelist(self, load_label, element_data):
        """
        Updates the VTK grid's cell data array corresponding to the given load label
        with the provided element data. Defaults to np.nan if a value is not present.

        Args:
            load_label (str): The label of the load (e.g., "Js").
            element_data (dict): A dictionary where keys are element IDs and values are dictionaries
                                of component values (0, 1, 2).
        """

        # Get the VTK array corresponding to the load label
        load_array = self.mesh.GetCellData().GetArray(load_label)
        if load_array is None:
            raise ValueError(f"Array '{load_label}' not found in the VTK file.")

        # Get the "ansys_elem_num" array to map element IDs
        ansys_elem_num_array = self.mesh.GetCellData().GetArray("ansys_elem_num")
        if ansys_elem_num_array is None:
            raise ValueError("Array 'ansys_elem_num' not found in the VTK file.")

        # Update the load array based on the element data
        for i in range(self.mesh.GetNumberOfCells()):
            elem_num = int(ansys_elem_num_array.GetValue(i))
            if elem_num in element_data:
                js_values = element_data[elem_num]
                for component, value in js_values.items():
                    load_array.SetComponent(
                        i, component, value if not np.isnan(value) else np.nan
                    )
            # else:
            #     # If the element is not in the data, set all components to NaN
            #     for component in range(load_array.GetNumberOfComponents()):
            #         load_array.SetComponent(i, component, np.nan)

    def read_bfelist(self, file_path):
        """
        Combines parsing element load data and updating the VTK grid.
        Parses the element load data from the given file and populates the VTK grid.

        Args:
            file_path (str): Path to the file containing element load data.

        Returns:
            str: The label of the load (e.g., "Js").
        """
        # Parse the element load data
        load_label, element_data = self.parse_element_bfelist_data(file_path)

        # Update the VTK grid with the parsed data
        self.update_vtk_grid_bfelist(load_label, element_data)

        print(f"Grid updated with element loads from '{file_path}'.")
        return load_label

    def read_write_bfelist(self, file_path, time_step=None):
        """
        Combines reading element load data, updating the VTK grid, and writing the loads to a binary file.

        Args:
            file_path (str): Path to the file containing element load data.
            time_step (int, optional): Time step to include in the binary file name. Defaults to None.
        """
        # Step 1: Read and update the VTK grid with element load data, and get the load label
        load_label = self.read_bfelist(file_path)

        # Step 2: Write the updated loads to a binary file
        self.write_loads_to_binary(load_label, time_step)

        print(f"Processed and wrote loads from '{file_path}' to binary file.")

    def read_bfelist_files(self):
        """
        Processes a series of bfelist files for the current jobname and writes one binary load file for each time step.

        Uses self.jobname as the basename to identify the files.

        Raises:
            FileNotFoundError: If no bfelist files are found for the current jobname.
        """

        # Ensure jobname is set
        if not self.jobname:
            raise ValueError(
                "Jobname is not set. Please set self.jobname before calling this function."
            )

        # Get the directory containing the files
        directory = os.getcwd()

        # Initialize a list to store matching files
        bfelist_files = []

        # Iterate over files in the directory
        for file_name in os.listdir(directory):
            # Check if the file name matches the expected pattern
            if (
                file_name.startswith(self.jobname)
                and file_name.endswith(".out")
                and "_bfelist_" in file_name
            ):
                bfelist_files.append(file_name)

        # Sort the files to ensure they are processed in order
        bfelist_files.sort()

        if not bfelist_files:
            raise FileNotFoundError(
                f"No bfelist files found for jobname '{self.jobname}'."
            )

        print(f"Found {len(bfelist_files)} bfelist files for jobname '{self.jobname}'.")

        # Process each file and write the corresponding binary load file
        for time_step, file_name in enumerate(bfelist_files, start=1):
            file_path = os.path.join(directory, file_name)
            print(f"Processing file '{file_path}' for time step {time_step}...")
            self.read_write_bfelist(file_path, time_step)

        print(f"Processed all bfelist files for jobname '{self.jobname}'.")

    def parse_node_dlist_data(self, file_path):
        """
        Parses node data from the given file.

        Args:
            file_path (str): Path to the file containing node data.

        Returns:
            dict: A dictionary where keys are node IDs and values are dictionaries
                of variable labels with their corresponding components and values.
                Components are 0, 1, 2 for vector variables and 0 for scalar variables.
        """

        # Regular expression to match node data lines
        node_pattern = re.compile(r"^\s*(\d+)\s+(\S+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)")

        # Dictionary to store parsed node data
        node_data = {}

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
                    real = float(node_match.group(3))

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

                    if node_num not in node_data:
                        node_data[node_num] = {}

                    if label not in node_data[node_num]:
                        node_data[node_num][label] = {}

                    node_data[node_num][label][component] = real

        return node_data

    def update_vtk_grid_dlist(self, node_data):
        """
        Updates the VTK grid's point data arrays based on the parsed node data.

        Args:
            node_data (dict): A dictionary where keys are node IDs and values are dictionaries
                            of variable labels with their corresponding components and values.
                            Scalars are represented with component 0.
        """
        # Get the "ansys_node_num" array to map node IDs
        ansys_node_num_array = self.mesh.GetPointData().GetArray("ansys_node_num")
        if ansys_node_num_array is None:
            raise ValueError("Array 'ansys_node_num' not found in the VTK file.")

        # Iterate over the points in the VTK grid
        for i in range(self.mesh.GetNumberOfPoints()):
            node_num = int(ansys_node_num_array.GetValue(i))
            if node_num in node_data:
                node_values = node_data[node_num]

                for label, components in node_values.items():
                    # Get the corresponding VTK array
                    vtk_array = self.mesh.GetPointData().GetArray(label)
                    if vtk_array is None:
                        raise ValueError(f"Array '{label}' not found in the VTK file.")

                    # Update the VTK array with the parsed values
                    for component, value in components.items():
                        vtk_array.SetComponent(i, component, value)

    def read_dlist(self, file_path):
        """
        Combines parsing node data and updating the VTK grid.
        Parses the node data from the given file and populates the VTK grid.

        Args:
            file_path (str): Path to the file containing node data.

        Returns:
            list: A list of variable labels found in the file.
        """
        # Parse the node data
        node_data = self.parse_node_dlist_data(file_path)

        # Update the VTK grid with the parsed data
        self.update_vtk_grid_dlist(node_data)

        # Extract variable labels from the parsed data
        variable_labels = set()
        for node_values in node_data.values():
            variable_labels.update(node_values.keys())

        print(f"Grid updated with node data from '{file_path}'.")
        return list(variable_labels)

    def read_write_dlist(self, file_path, time_step=None):
        """
        Combines reading node data, updating the VTK grid, and writing the boundary conditions to a binary file.

        Args:
            file_path (str): Path to the file containing node data.
            time_step (int, optional): Time step to include in the binary file name. Defaults to None.
        """
        # Step 1: Read and update the VTK grid with node data, and get the variable labels
        variable_labels = self.read_dlist(file_path)

        # Step 2: Write the updated boundary conditions to a binary file for each variable label
        for label in variable_labels:
            self.write_boundary_conditions_to_binary(label, time_step)

        print(
            f"Processed and wrote boundary conditions from '{file_path}' to binary file."
        )

    def read_dlist_files(self):
        """
        Processes a series of dlist files for the current jobname and writes one binary boundary condition file for each time step.

        Uses self.jobname as the basename to identify the files.

        Raises:
            FileNotFoundError: If no dlist files are found for the current jobname.
        """
        # Ensure jobname is set
        if not self.jobname:
            raise ValueError(
                "Jobname is not set. Please set self.jobname before calling this function."
            )

        # Get the directory containing the files
        directory = os.getcwd()

        # Initialize a list to store matching files
        dlist_files = []

        # Iterate over files in the directory
        for file_name in os.listdir(directory):
            # Check if the file name matches the expected pattern
            if (
                file_name.startswith(self.jobname)
                and file_name.endswith(".out")
                and "_dlist_" in file_name
            ):
                dlist_files.append(file_name)

        # Sort the files to ensure they are processed in order
        dlist_files.sort()

        if not dlist_files:
            raise FileNotFoundError(
                f"No dlist files found for jobname '{self.jobname}'."
            )

        print(f"Found {len(dlist_files)} dlist files for jobname '{self.jobname}'.")

        # Process each file and write the corresponding binary boundary condition file
        for time_step, file_name in enumerate(dlist_files, start=1):
            file_path = os.path.join(directory, file_name)
            print(f"Processing file '{file_path}' for time step {time_step}...")
            self.read_write_dlist(file_path, time_step)

        print(f"Processed all dlist files for jobname '{self.jobname}'.")


class ElementType(ABC):
    """
    Abstract base class for all element types.
    """

    @abstractmethod
    def getVariables(self):
        pass

    @abstractmethod
    def getElementLoads(self):
        pass


class SOLID97Type(ElementType):
    def getVariables(self):
        return {"A": 3, "V": 1}

    def getElementLoads(self):
        return {"Js": 3}


class SOLID236Type(ElementType):
    def getVariables(self):
        return {"Az": 1, "V": 1}

    def getElementLoads(self):
        return {"Js": 3}


class PLANE13Type(ElementType):
    def getVariables(self):
        return {"A": 3}

    def getElementLoads(self):
        return {"Js": 3}


class PLANE233Type(ElementType):
    def getVariables(self):
        return {"A": 3}

    def getElementLoads(self):
        return {"Js": 3}


class ElementTypesCollector:
    """
    Class to collect and manage element types from a CSV file.
    """

    def __init__(self, apec):
        """
        Initializes the instance with the given APEC data.
        Args:
            apec: The APEC data to be associated with the instance.
        """

        self.apec = apec
        self.element_types_map = {}

    def read_element_types(self):

        for _, row in apec.element_types.iterrows():
            element_type = row["element_type"]
            ansys_etype = row["ansys_etype"]
            self.element_types_map[ansys_etype] = self.build_element_types(element_type)

    def build_element_types(self, element_type):
        """
        Maps a given element type to its corresponding ANSYS element type object and
        stores it in the `element_types_map` dictionary.
        Parameters:
            element_type (str): The name of the element type to map.
                                Supported values are "SOLID97", "SOLID236", "PLANE13", and "PLANE233".
        Raises:
            ValueError: If the provided `element_type` is not recognized.
        """
        element_type_string = element_type + "Type"

        try:
            return globals()[element_type_string]()
        except KeyError:
            raise ValueError(f"Unknown element type: {element_type}")

    def collect_type_classes(self):
        """
        Collects all class names in the current script that have 'Type' in their name
        and creates instances of those classes.

        Returns:
            dict: A dictionary where keys are class names (strings) and values are instances of those classes.
        """
        class_instances_with_type = {
            name: cls()
            for name, cls in globals().items()
            if isinstance(cls, type) and "Type" in name
        }
        return class_instances_with_type

    def get_all_loads(self):
        all_loads = {}
        for element_type in self.element_types_map.values():
            all_loads.update(element_type.getElementLoads())
        return all_loads


apec = APEC()
