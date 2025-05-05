import pandas as pd
from abc import ABC, abstractmethod
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import struct
import os
import re
import subprocess
# import logging

VARIABLES_GLOBAL_NUMBERS = {"A": 0, "VOLT": 1, "Az": 2}
VARIABLES_LOADS_NUMBERS = {"Js": 0}

DEFAULT_SETTINGS = {
    "ksp_type": "kspreonly",
    "pc_type": "lu",
    "initial_linear_tolerance": 1e-6,
    "max_linear_iterations": 100,
    "minimum_linear_tolerance": 1e-10,
    "linear_tolerance_multiplier": 1.0,
    "theta": 0.5,
    "write_full_matrix": 0,
    "timesolver_tolerance": 1e-6,
    "analytic_jacobians": 0,
    "verify_analytic_jacobians": 0,
    "numerical_jacobian_h": 1e-6,
    "quiet": 0,
    "max_nonlinear_it": 25,
    "min_steplength": 1e-6,
    "relative_step_tolerance": 1e-6,
    "require_residual_reduction": 0,
    "reduce_deltat_on_diffsolver_failure": 1,
    "continue_after_max_iterations": 0,
    "continue_after_backtrack_failure": 0,
    "timesolver_core": 1,
    "print_solution_norms": 0,
    "print_solutions": 0,
    "print_residual_norms": 0,
    "print_residuals": 0,
    "print_jacobian_norms": 0,
    "print_jacobians": 0,
    "print_element_solutions": 0,
    "print_element_residuals": 0,
    "print_element_jacobians": 0,
    "constrain_in_solver": 0,
    "use_petsc_snes": 0,
    "min_step_length": 1e-6,
}


class APEC_M:
    """
    APEC class to manage VTK unstructured grids and element types.
    """

    def __init__(self):
        """
        Initialize the APEC object without requiring a file path.
        """
        self.mesh = vtk.vtkUnstructuredGrid()

        self.initialize_dataframes()
        self.jobname = None
        self.folderpath = None
        self.solver_settings = DEFAULT_SETTINGS

        self.process = None
        self.num_processors = 1

    def start_process(self, command, num_processors=1):
        """
        Start a process and set it as an attribute. If num_processors > 1, use mpirun.

        Args:
            command (list): The command to execute as a list of strings.
            num_processors (int): The number of processors to use. Defaults to 1.
        """
        self.num_processors = num_processors
        command = command.split(" ")
        # Prepend mpirun if num_processors > 1
        if self.num_processors > 1:
            command = ["mpirun", "-np", str(self.num_processors)] + command

        self.process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        print(f"Process started: {' '.join(command)}")
        while True:
            line = self.process.stdout.readline()
            if not line:
                break
            print(line.strip())
            if line.strip().endswith(">>>"):  # Example stop condition for Python prompt
                break
        # line = self.process.stdout.readline()
        # print(">>", line.strip())

    #     # Start a thread to read the process output
    #     threading.Thread(target=self._read_output, daemon=True).start()

    # def _read_output(self):
    #     """
    #     Continuously read and print the process output.
    #     """
    #     if self.process and self.process.stdout:
    #         for line in iter(self.process.stdout.readline, ""):
    #             print(line, end="")  # Print the output line by line

    def exec(self, command):
        """
        Send a command to the running process.

        Args:
            command (str): The command to send to the process.
        """
        if self.process and self.process.stdin:
            self.process.stdin.write(command + "\n")
            self.process.stdin.flush()
            print(f"Command sent: {command}")
        else:
            print("No process is running.")
        # Read lines until ready for next input (may need a better stop condition)
        while True:
            line = self.process.stdout.readline()
            if not line:
                break
            print(line.strip())
            if line.strip().endswith(">>>"):  # Example stop condition for Python prompt
                break

    def exit(self):
        """
        Stop the running process.
        """
        if self.process:
            self.exec("-exit")  # Send the exit command to the process
            self.process.wait()  # Wait for the process to terminate
            print("Process stopped.")
            self.process = None
        else:
            print("No process is running.")

    def start(self, num_processors=1):
        """
        Run an interactive session with the process.

        Args:
            num_processors (int): The number of processors to use. Defaults to 1.
        """
        try:
            # Start the process
            self.start_process(["apec_m", "-i"], num_processors=num_processors)

            print("Interactive session started with `apec_m -i`.")
            print(
                "Enter commands (e.g., -model arg, -solu, -solve arg, -finish, -exit)."
            )
            print("Type 'exit' to terminate the session.")

            # Interactive loop for sending commands
            while True:
                command = input("> ")  # Prompt for user input
                if command.strip().lower() == "exit":
                    self.exit()
                    break
                self.exec(command)

        except FileNotFoundError:
            print("Error: `apec_m` executable not found. Ensure it is in your PATH.")
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            if self.process and self.process.poll() is None:
                self.process.terminate()  # Ensure the process is terminated

    def initialize_dataframes(self):
        """
        Initialize the DataFrames with the proper columns.
        """
        self.materials = pd.DataFrame(
            columns=[
                "ansys_material_type",
                "MURX",
                "MURY",
                "MURZ",
                "RSVX",
                "RSVY",
                "RSVZ",
                "PERX",
                "PERY",
                "PERZ",
                "MGXX",
                "MGYY",
                "MGZZ",
                "MURX_B-H",
                "MURY_B-H",
                "MURZ_B-H"
            ]
        )
        self.real_constants = pd.DataFrame(
            columns=[
                "ansys_real_constant",
                "CARE",
                "TURN",
                "LENG",
                "DIRZ",
                "FILL",
                "SC",
                "NC",
                "RAD",
                "TX",
                "TY",
                "TZ",
                "R",
                "SYM",
                "VELOX",
                "VELOY",
                "VELOZ",
                "OMEGAX",
                "OMEGAY",
                "OMEGAZ",
                "XLOC",
                "YLOC",
                "ZLOC",
            ]
        )
        self.time_parameters = {
            "time_step_number": 1,
            "time": 0.0,
            "number_of_substeps": 1,
            "transient_effects": 1,
            "ramped_loads": 1,
            "write_substep_results": 0,
        }
        self.apec_parameters = {}

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

    def dataframe_to_vtk_table(self, df):
        """
        Converts a Pandas DataFrame into a vtkTable.

        Parameters:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            vtk.vtkTable: The converted vtkTable.
        """
        vtk_table = vtk.vtkTable()

        for column_name in df.columns:
            # Determine the data type of the column
            if pd.api.types.is_numeric_dtype(df[column_name]):
                vtk_array = vtk.vtkDoubleArray()
            else:
                vtk_array = vtk.vtkStringArray()

            vtk_array.SetName(column_name)

            # Add data to the vtk array
            for value in df[column_name]:
                if isinstance(vtk_array, vtk.vtkDoubleArray):
                    vtk_array.InsertNextValue(float(value))
                else:
                    vtk_array.InsertNextValue(str(value))

            # Add the vtk array to the vtkTable
            vtk_table.AddColumn(vtk_array)

        return vtk_table

    def add_table_to_mesh(self, vtk_table):
        """
        Adds a vtkTable as field data to self.mesh.
        Deletes any existing field data arrays with names matching the columns of the vtkTable.

        Parameters:
            vtk_table (vtk.vtkTable): The vtkTable to add.
        """
        # Get the field data from the unstructured grid
        field_data = self.mesh.GetFieldData()
        if not field_data:
            field_data = vtk.vtkFieldData()
            self.mesh.SetFieldData(field_data)

        # Iterate over the columns in the vtkTable
        for i in range(vtk_table.GetNumberOfColumns()):
            column_name = vtk_table.GetColumnName(i)

            # Check if a field array with the same name exists and remove it
            if field_data.HasArray(column_name):
                field_data.RemoveArray(column_name)

            # Add the column as a field array
            field_data.AddArray(vtk_table.GetColumn(i))

    def read_element_types(self, csv_file_path):
        """
        Reads the element types from a CSV file, converts them to a vtkTable,
        and populates the field data of self.mesh.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        # Step 1: Read the CSV file into a DataFrame
        df = self.read_dataframe_from_csv(csv_file_path)

        # Step 2: Check that the DataFrame has the expected columns
        expected_columns = [
            "ansys_etype",
            "element_type",
            "KEYOPT1",
            "KEYOPT2",
            "KEYOPT3",
            # "KEYOPT4",
            # "KEYOPT5",
            # "KEYOPT6",
            # "KEYOPT7",
            # "KEYOPT8",
        ]
        if not all(column in df.columns for column in expected_columns):
            raise ValueError(
                f"CSV file must contain the following columns: {expected_columns}"
            )

        # Step 3: Convert the DataFrame to a vtkTable
        vtk_table = self.dataframe_to_vtk_table(df)

        # Step 4: Add the vtkTable to the mesh as field data
        self.add_table_to_mesh(vtk_table)

        print(
            f"Element types loaded and added to mesh field data from '{csv_file_path}'."
        )

    def get_dataframe_from_field_data(self):
        """
        Extracts a pandas DataFrame from the field data of self.mesh.

        Returns:
            pd.DataFrame: A DataFrame containing the field data.

        Raises:
            ValueError: If the field data is missing in the mesh.
        """
        # Get the field data from the mesh
        field_data = self.mesh.GetFieldData()
        if not field_data:
            raise ValueError("Field data is missing in the mesh.")

        # Initialize a dictionary to store the data
        data = {}

        # Iterate over all arrays in the field data
        for i in range(field_data.GetNumberOfArrays()):
            vtk_array = field_data.GetArray(i)
            column_name = vtk_array.GetName()

            # Convert the VTK array to a NumPy array and store it in the dictionary
            data[column_name] = vtknp.vtk_to_numpy(vtk_array)

        # Convert the dictionary to a pandas DataFrame
        df = pd.DataFrame(data)

        return df

    def read_materials(self, csv_file_path):
        """
        Read the materials DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.materials = self.read_dataframe_from_csv(csv_file_path)
        print(f"Materials loaded from '{csv_file_path}'.")

    def read_real_constants(self, csv_file_path):
        """
        Read the real constants DataFrame from a CSV file.

        Parameters:
            csv_file_path (str): Path to the CSV file.
        """
        self.real_constants = self.read_dataframe_from_csv(csv_file_path)
        print(f"Real constants loaded from '{csv_file_path}'.")

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
        self.load_steps_folder = os.path.join(os.getcwd(), "load_steps")
        if not os.path.exists(self.load_steps_folder):
            os.makedirs(self.load_steps_folder)

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

    # def read_solver_settings_from_json(self, json_file_path):
    #     """
    #     Reads solver settings from a JSON file and updates the current settings.

    #     Args:
    #         json_file_path (str): Path to the JSON file containing the solver settings.

    #     Raises:
    #         FileNotFoundError: If the JSON file does not exist.
    #     """

    #     if not os.path.exists(json_file_path):
    #         raise FileNotFoundError(f"JSON file '{json_file_path}' not found.")

    #     with open(json_file_path, "r") as json_file:
    #         self.solver_settings = json.load(json_file)
    #     print(f"Solver settings loaded from '{json_file_path}'.")

    # def write_solver_settings_to_json(self, json_file_path):
    #     """
    #     Writes the solver settings to a JSON file.

    #     Args:
    #         json_file_path (str): Path to the JSON file where the solver settings will be saved.
    #     """

    #     with open(json_file_path, "w") as json_file:
    #         json.dump(self.solver_settings, json_file, indent=4)
    #     print(f"Solver settings written to '{json_file_path}'.")

    def write_settings_to_binary(self, settings, file_path):
        """
        Writes a settings dictionary to a binary file. If the file exists, appends to it.

        Parameters:
            settings (dict): The settings dictionary to write.
            file_path (str): The path to the binary file.
        """
        # Open the file in append mode if it exists, otherwise create a new file
        mode = "ab" if os.path.exists(file_path) else "wb"
        with open(file_path, mode) as binary_file:
            # Write the number of settings
            binary_file.write(struct.pack("i", len(settings)))

            for key, value in settings.items():
                # Write the size of the string label and the label itself
                binary_file.write(struct.pack("i", len(key)))
                binary_file.write(key.encode("utf-8"))

                # Determine the type of the value and write it
                if isinstance(value, int):
                    binary_file.write(struct.pack("i", 0))  # Type: int
                    binary_file.write(struct.pack("i", value))
                elif isinstance(value, float):
                    binary_file.write(struct.pack("i", 1))  # Type: double
                    binary_file.write(struct.pack("d", value))
                elif isinstance(value, str):
                    binary_file.write(struct.pack("i", 2))  # Type: string
                    binary_file.write(struct.pack("i", len(value)))
                    binary_file.write(value.encode("utf-8"))
                else:
                    raise ValueError(f"Unsupported value type: {type(value)}")

    def write_dataframe_to_binary(df, file_path):
        """
        Writes a pandas DataFrame to a binary file.

        The binary file format includes:
            - Number of rows (unsigned int)
            - Number of columns (unsigned int)
            - Data values (row by row):
                - First column as unsigned int (ID)
                - Remaining columns as double (real numbers)

        Parameters:
            df (pd.DataFrame): The DataFrame to write.
            file_path (str): The path to the binary file.
        """
        # Ensure the first column is treated as an unsigned integer ID
        if not pd.api.types.is_integer_dtype(df.iloc[:, 0]):
            raise ValueError("The first column must contain integer IDs.")

        # Open the file in write mode
        with open(file_path, "wb") as binary_file:
            # Write the number of rows and columns
            num_rows = len(df)
            num_columns = len(df.columns)
            binary_file.write(struct.pack("I", num_rows))  # Number of rows (unsigned int)
            binary_file.write(struct.pack("I", num_columns))  # Number of columns (unsigned int)

            # Write the data row by row
            for _, row in df.iterrows():
                # Write the first column as an unsigned integer (ID)
                binary_file.write(struct.pack("I", int(row.iloc[0])))

                # Write the remaining columns as double (real numbers)
                for value in row.iloc[1:]:
                    binary_file.write(struct.pack("d", float(value)))

    def write_loads_to_binary(self, output_file):
        """
        Writes a binary file for loads from a VTK cell data array.
        Handles multi-component arrays by creating structured arrays for load numbers, components,
        cell IDs, and values, and writes them efficiently.

        Writes:
            - Number of load rows (unsigned int)
            - For each load row:
                - cell_id (unsigned int)
                - load_number (unsigned int)
                - component (unsigned int)
                - value (double)
        """
        # Get the cell data array
        cell_data = self.mesh.GetCellData()

        # Get the load labels from the VARIABLES_LOADS_NUMBERS mapping
        loads_labels = VARIABLES_LOADS_NUMBERS.keys()

        # Open the output file in append mode if it exists, otherwise create a new file
        mode = "ab" if os.path.exists(output_file) else "wb"
        with open(output_file, mode) as file:
            total_rows = 0  # To keep track of the total number of rows to write
            all_stacked_data = []  # Temporary storage for structured arrays

            # Define the structured array data type
            dtype = np.dtype([
                ('cell_id', np.uint32),
                ('load_number', np.uint32),
                ('component', np.uint32),
                ('value', np.float64)
            ])

            # Iterate over each load label
            for load_label in loads_labels:
                load_array = cell_data.GetArray(load_label)
                if load_array is None:
                    continue  # Skip if the array does not exist in the VTK file

                # Convert the VTK array to a NumPy array
                load_array_np = vtknp.vtk_to_numpy(load_array)

                # Get the number of components and cells
                num_components = load_array.GetNumberOfComponents()
                num_cells = self.mesh.GetNumberOfCells()

                # Iterate over components to filter data
                for component in range(num_components):
                    # Extract the values for the current component
                    component_values = load_array_np[:, component]

                    # Create a mask to filter out NaN values
                    mask = ~np.isnan(component_values)

                    # Get the filtered cell IDs and values
                    filtered_cell_ids = np.arange(num_cells, dtype=np.uint32)[mask]
                    filtered_values = component_values[mask]

                    # Create arrays for load number and component
                    load_number_array = np.full(len(filtered_cell_ids), VARIABLES_LOADS_NUMBERS[load_label], dtype=np.uint32)
                    component_array = np.full(len(filtered_cell_ids), component, dtype=np.uint32)

                    # Create a structured array for the current component
                    structured_data = np.zeros(len(filtered_cell_ids), dtype=dtype)
                    structured_data['cell_id'] = filtered_cell_ids
                    structured_data['load_number'] = load_number_array
                    structured_data['component'] = component_array
                    structured_data['value'] = filtered_values

                    # Append the structured array to the list
                    all_stacked_data.append(structured_data)

                    # Update the total number of rows to write
                    total_rows += len(filtered_cell_ids)

            # Concatenate all structured arrays into a single array
            if all_stacked_data:
                final_data = np.concatenate(all_stacked_data)

                # Write the total number of rows to the file
                file.write(struct.pack("I", total_rows))  # 'I' is for unsigned int

                # Write the final structured data to the file
                file.write(final_data.tobytes())

    def write_boundary_conditions_to_binary(self, output_file):
        """
        Writes a binary file for boundary conditions from a VTK point data array.
        Handles multi-component arrays by creating structured arrays for variable numbers, components,
        point IDs, and values, and writes them efficiently.

        Writes:
            - Number of boundary condition rows (unsigned int)
            - For each boundary condition row:
                - point_id (unsigned int)
                - variable_number (unsigned int, from VARIABLES_GLOBAL_NUMBERS mapping)
                - component (unsigned int)
                - value (double)
        """
        # Get the point data array
        point_data = self.mesh.GetPointData()

        # Get the variable labels from the VARIABLES_GLOBAL_NUMBERS mapping
        variable_labels = VARIABLES_GLOBAL_NUMBERS.keys()

        # Open the output file in append mode if it exists, otherwise create a new file
        mode = "ab" if os.path.exists(output_file) else "wb"
        with open(output_file, mode) as file:
            total_rows = 0  # To keep track of the total number of rows to write
            all_stacked_data = []  # Temporary storage for structured arrays

            # Define the structured array data type
            dtype = np.dtype([
                ('point_id', np.uint32),
                ('variable_number', np.uint32),
                ('component', np.uint32),
                ('value', np.float64)
            ])

            # Iterate over each variable label
            for variable_label in variable_labels:
                variable_array = point_data.GetArray(variable_label)
                if variable_array is None:
                    continue  # Skip if the array does not exist in the VTK file

                # Convert the VTK array to a NumPy array
                variable_array_np = vtknp.vtk_to_numpy(variable_array)

                # Get the number of components and points
                num_components = variable_array.GetNumberOfComponents()
                num_points = self.mesh.GetNumberOfPoints()

                # Iterate over components to filter data
                for component in range(num_components):
                    # Extract the values for the current component
                    component_values = variable_array_np[:, component]

                    # Create a mask to filter out NaN values
                    mask = ~np.isnan(component_values)

                    # Get the filtered point IDs and values
                    filtered_point_ids = np.arange(num_points, dtype=np.uint32)[mask]
                    filtered_values = component_values[mask]

                    # Create arrays for variable number and component
                    variable_number_array = np.full(len(filtered_point_ids), VARIABLES_GLOBAL_NUMBERS[variable_label], dtype=np.uint32)
                    component_array = np.full(len(filtered_point_ids), component, dtype=np.uint32)

                    # Create a structured array for the current component
                    structured_data = np.zeros(len(filtered_point_ids), dtype=dtype)
                    structured_data['point_id'] = filtered_point_ids
                    structured_data['variable_number'] = variable_number_array
                    structured_data['component'] = component_array
                    structured_data['value'] = filtered_values

                    # Append the structured array to the list
                    all_stacked_data.append(structured_data)

                    # Update the total number of rows to write
                    total_rows += len(filtered_point_ids)

            # Concatenate all structured arrays into a single array
            if all_stacked_data:
                final_data = np.concatenate(all_stacked_data)

                # Write the total number of rows to the file
                file.write(struct.pack("I", total_rows))  # 'I' is for unsigned int

                # Write the final structured data to the file
                file.write(final_data.tobytes())

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

    def read_parse_bfelist(self, file_path):
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
                            data.append(
                                {
                                    "ansys_elem_num": current_element,
                                    "load_label": load_label,
                                    "component": component,
                                    "value": value,
                                }
                            )
                        current_element = None  # Reset current element after processing

        # Convert the data to a pandas DataFrame
        df = pd.DataFrame(
            data, columns=["ansys_elem_num", "load_label", "component", "value"]
        )
        return df

    def map_ansys_elem_num_to_vtk(self):
        """
        Maps the ANSYS element numbers to VTK cell IDs and returns a 2D NumPy array.
        This function assumes that the VTK grid has been loaded and contains the "ansys_elem_num" array.

        Returns:
            np.ndarray: A 2D NumPy array where each row contains:
                        [VTK cell ID, ANSYS element number].
        """
        # Get the "ansys_elem_num" array from the VTK grid
        ansys_elem_num_array = self.mesh.GetCellData().GetArray("ansys_elem_num")
        if ansys_elem_num_array is None:
            raise ValueError("Array 'ansys_elem_num' not found in the VTK file.")

        # Convert the VTK array to a NumPy array
        ansys_elem_nums = vtknp.vtk_to_numpy(ansys_elem_num_array)

        # Create a 2D NumPy array with VTK IDs and ANSYS element numbers
        vtk_ids = np.arange(len(ansys_elem_nums))  # VTK IDs are sequential from 0
        mapping_array = np.column_stack((vtk_ids, ansys_elem_nums))

        return mapping_array

    def replace_ansys_with_vtk_ids(self, df, mapping_array):
        """
        Replaces the 'ansys_elem_num' column in the DataFrame with the corresponding VTK ID values
        using the provided mapping array. The column is renamed to 'vtk_element_id'.

        Args:
            df (pd.DataFrame): The DataFrame containing the 'ansys_elem_num' column.
            mapping_array (np.ndarray): A 2D NumPy array where each row contains [VTK ID, ANSYS element number].

        Returns:
            pd.DataFrame: A new DataFrame with the 'ansys_elem_num' column replaced by VTK IDs and renamed to 'vtk_element_id'.
        """
        # Extract the 'ansys_elem_num' column as a NumPy array
        ansys_elem_nums = df["ansys_elem_num"].to_numpy()

        # Extract the VTK IDs using NumPy indexing
        vtk_ids = mapping_array[
            :, 0
        ]  # VTK IDs are in the first column of the mapping array
        vtk_id_array = vtk_ids[ansys_elem_nums]

        # Replace the column in the DataFrame and rename it to 'vtk_element_id'
        df.insert(0, "vtk_element_id", vtk_id_array)  # Insert 'vtk_element_id' as the first column
        df = df.drop(columns=["ansys_elem_num"])

        return df

    def update_vtk_grid_loads(self, df):
        """
        Updates the VTK grid's cell data arrays using the DataFrame resulting from `replace_ansys_with_vtk_ids`.
        The function assigns values to the appropriate components at the corresponding VTK IDs.

        Args:
            df (pd.DataFrame): The DataFrame containing columns 'vtk_element_id', 'load_label', 'component', and 'value'.
        """
        # Get the distinct load labels
        load_labels = df["load_label"].unique()

        for load_label in load_labels:
            # Filter the DataFrame for the current load label
            load_df = df[df["load_label"] == load_label]

            # Get the corresponding VTK array
            vtk_array = self.mesh.GetCellData().GetArray(load_label)
            if vtk_array is None:
                raise ValueError(f"Array '{load_label}' not found in the VTK file.")

            # Convert the VTK array to a NumPy array
            vtk_array_np = vtknp.vtk_to_numpy(vtk_array)

            # Iterate over distinct component values
            for component in load_df["component"].unique():
                # Filter rows for the current component
                component_mask = load_df["component"] == component
                vtk_ids = load_df.loc[component_mask, "vtk_element_id"].to_numpy()
                values = load_df.loc[component_mask, "value"].to_numpy()

                # Assign values to the appropriate indices and component
                vtk_array_np[vtk_ids, component] = values

            # Convert the NumPy array back to a VTK array
            vtk_array.Modified()

    def read_csv_loads(self, csv_file_path):
        """
        Reads a CSV file containing load data and updates the VTK grid's cell data arrays.
        The CSV file must have the following columns in order:
            - vtk_element_id: The VTK cell ID.
            - load_label: The label associated with the load (e.g., "Js").
            - component: The component index (0, 1, 2).
            - value: The value of the load for the given component.

        Args:
            csv_file_path (str): Path to the CSV file containing the load data.

        Raises:
            FileNotFoundError: If the CSV file does not exist.
        """
        # Check if the CSV file exists
        if not os.path.exists(csv_file_path):
            raise FileNotFoundError(f"CSV file '{csv_file_path}' not found.")

        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file_path)

        # Ensure the required columns are present
        required_columns = ["vtk_element_id", "load_label", "component", "value"]
        if not all(column in df.columns for column in required_columns):
            raise ValueError(
                f"CSV file must contain the following columns: {required_columns}"
            )

        # Call update_vtk_grid_loads to update the VTK grid
        self.update_vtk_grid_loads(df)

        print(f"Grid updated with loads from '{csv_file_path}'.")

    def read_bfelist(self, file_path):
        """
        Combines parsing element load data, mapping ANSYS element numbers to VTK IDs,
        replacing ANSYS element numbers with VTK IDs in the DataFrame, and updating the VTK grid.

        Args:
            file_path (str): Path to the file containing element load data.

        Returns:
            str: The label of the load (e.g., "Js").
        """
        # Step 1: Parse the element load data into a DataFrame
        df = self.read_parse_bfelist(file_path)

        # Step 2: Map ANSYS element numbers to VTK IDs
        mapping_array = self.map_ansys_elem_num_to_vtk()

        # Step 3: Replace ANSYS element numbers with VTK IDs in the DataFrame
        df = self.replace_ansys_with_vtk_ids(df, mapping_array)

        # Step 4: Update the VTK grid with the modified DataFrame
        self.update_vtk_grid_loads(df)

        print(f"Grid updated with element loads from '{file_path}'.")

    def read_parse_dlist(self, file_path):
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
                            "variable_label": label,
                            "component": component,
                            "value": value,
                        }
                    )

        # Convert the data to a pandas DataFrame
        df = pd.DataFrame(
            data, columns=["ansys_node_num", "variable_label", "component", "value"]
        )
        return df

    def map_ansys_node_num_to_vtk(self):
        """
        Maps the ANSYS node numbers to VTK point IDs and returns a 2D NumPy array.
        This function assumes that the VTK grid has been loaded and contains the "ansys_node_num" array.

        Returns:
            np.ndarray: A 2D NumPy array where each row contains:
                        [VTK point ID, ANSYS node number].
        """
        # Get the "ansys_node_num" array from the VTK grid
        ansys_node_num_array = self.mesh.GetPointData().GetArray("ansys_node_num")
        if ansys_node_num_array is None:
            raise ValueError("Array 'ansys_node_num' not found in the VTK file.")

        # Convert the VTK array to a NumPy array
        ansys_node_nums = vtknp.vtk_to_numpy(ansys_node_num_array)

        # Create a 2D NumPy array with VTK point IDs and ANSYS node numbers
        vtk_ids = np.arange(len(ansys_node_nums))  # VTK point IDs are sequential from 0
        mapping_array = np.column_stack((vtk_ids, ansys_node_nums))

        return mapping_array

    def replace_ansys_with_vtk_point_ids(self, df, mapping_array):
        """
        Replaces the 'ansys_node_num' column in the DataFrame with the corresponding VTK point ID values
        using the provided mapping array. The column is renamed to 'vtk_node_id'.

        Args:
            df (pd.DataFrame): The DataFrame containing the 'ansys_node_num' column.
            mapping_array (np.ndarray): A 2D NumPy array where each row contains [VTK point ID, ANSYS node number].

        Returns:
            pd.DataFrame: A new DataFrame with the 'ansys_node_num' column replaced by VTK point IDs and renamed to 'vtk_node_id'.
        """
        # Extract the 'ansys_node_num' column as a NumPy array
        ansys_node_nums = df["ansys_node_num"].to_numpy()

        # Extract the VTK point IDs using NumPy indexing
        vtk_ids = mapping_array[
            :, 0
        ]  # VTK point IDs are in the first column of the mapping array
        vtk_id_array = vtk_ids[np.searchsorted(mapping_array[:, 1], ansys_node_nums)]

        # Replace the column in the DataFrame and rename it to 'vtk_node_id'
        df.insert(0, "vtk_node_id", vtk_id_array)  # Insert 'vtk_node_id' as the first column
        df = df.drop(columns=["ansys_node_num"])  # Drop the 'ansys_node_num' column

        return df

    def read_csv_bcs(self, csv_file_path):
        """
        Reads a CSV file containing boundary condition data and updates the VTK grid's point data arrays.
        The CSV file must have the following columns in order:
            - vtk_node_id: The VTK point ID.
            - variable_label: The label associated with the boundary condition (e.g., "A", "VOLT").
            - component: The component index (0, 1, 2 for vector variables, 0 for scalar variables).
            - value: The value of the boundary condition for the given component.

        Args:
            csv_file_path (str): Path to the CSV file containing the boundary condition data.

        Raises:
            FileNotFoundError: If the CSV file does not exist.
            ValueError: If the required columns are not present in the CSV file.
        """
        # Check if the CSV file exists
        if not os.path.exists(csv_file_path):
            raise FileNotFoundError(f"CSV file '{csv_file_path}' not found.")

        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file_path)

        # Ensure the required columns are present
        required_columns = ["vtk_node_id", "variable_label", "component", "value"]
        if not all(column in df.columns for column in required_columns):
            raise ValueError(
                f"CSV file must contain the following columns: {required_columns}"
            )

        # Call update_vtk_grid_bcs to update the VTK grid
        self.update_vtk_grid_bcs(df)

        print(f"Grid updated with boundary conditions from '{csv_file_path}'.")

    def update_vtk_grid_bcs(self, df):
        """
        Updates the VTK grid's point data arrays using the DataFrame resulting from `replace_ansys_with_vtk_point_ids`.
        The function assigns values to the appropriate components at the corresponding VTK point IDs.

        Args:
            df (pd.DataFrame): The DataFrame containing columns 'vtk_node_id', 'variable_label', 'component', and 'value'.
        """
        # Get the distinct variable labels
        variable_labels = df["variable_label"].unique()

        for variable_label in variable_labels:
            # Filter the DataFrame for the current variable label
            variable_df = df[df["variable_label"] == variable_label]

            # Get the corresponding VTK array
            vtk_array = self.mesh.GetPointData().GetArray(variable_label)
            if vtk_array is None:
                raise ValueError(f"Array '{variable_label}' not found in the VTK file.")

            # Convert the VTK array to a NumPy array
            vtk_array_np = vtknp.vtk_to_numpy(vtk_array)

            # Iterate over distinct component values
            for component in variable_df["component"].unique():
                # Filter rows for the current component
                component_mask = variable_df["component"] == component
                vtk_ids = variable_df.loc[component_mask, "vtk_node_id"].to_numpy()
                values = variable_df.loc[component_mask, "value"].to_numpy()

                # Assign values to the appropriate indices and component
                vtk_array_np[vtk_ids, component] = values

            # Mark the VTK array as modified
            vtk_array.Modified()

    def read_dlist(self, file_path):
        """
        Combines parsing node data, mapping ANSYS node numbers to VTK point IDs,
        replacing ANSYS node numbers with VTK point IDs in the DataFrame, and updating the VTK grid.

        Args:
            file_path (str): Path to the file containing node data.

        Returns:
            list: A list of variable labels found in the file.
        """
        # Step 1: Parse the node data into a DataFrame
        df = self.read_parse_dlist(file_path)

        # Step 2: Map ANSYS node numbers to VTK point IDs
        mapping_array = self.map_ansys_node_num_to_vtk()

        # Step 3: Replace ANSYS node numbers with VTK point IDs in the DataFrame
        df = self.replace_ansys_with_vtk_point_ids(df, mapping_array)

        # Step 4: Update the VTK grid with the modified DataFrame
        self.update_vtk_grid_bcs(df)

        # Return the list of variable labels
        return df["variable_label"].unique().tolist()

    def write_load_step(self, output_file):
        """
        Writes the current load step data to a single binary file.

        This includes:
            - Time parameters.
            - Solver settings.
            - Materials DataFrame.
            - Loads.
            - Boundary conditions.

        Parameters:
            output_file (str): The path to the output binary file.
        """
        # Open the file in write mode to overwrite if it already exists
        with open(output_file, "wb") as binary_file:
            # Write time parameters
            self.write_settings_to_binary(self.time_parameters, binary_file)
            print("Time parameters written.")

            # Write solver settings
            self.write_settings_to_binary(self.solver_settings, binary_file)
            print("Solver settings written.")

            # Write materials DataFrame
            self.write_dataframe_to_binary(self.materials, binary_file)
            print("Materials written.")

            # Write real constants DataFrame
            self.write_dataframe_to_binary(self.real_constants, binary_file)
            print("Real constants written.")

            # Write loads to binary
            self.write_loads_to_binary(binary_file)
            print("Loads written.")

            # Write boundary conditions to binary
            self.write_boundary_conditions_to_binary(binary_file)
            print("Boundary conditions written.")

    def read_load_step(
        self,
        solver_settings_json="solver_settings.json",
        materials_csv="materials.csv",
        real_constants_csv="real_constants.csv",
        loads_csv="loads_1.csv",
        bcs_csv="bcs_1.csv",
    ):
        """
        Initializes the ApecM object by reading data from various files.

        Parameters:
            solver_settings_json (str): Path to the JSON file containing solver settings (default: "solver_settings.json").
            materials_csv (str): Path to the CSV file containing materials data (default: "materials.csv").
            real_constants_csv (str): Path to the CSV file containing real constants data (default: "real_constants.csv").
            time_data_csv (str): Path to the CSV file containing time parameters (default: "time_data.csv").
            loads_csv (str): Path to the CSV file containing loads data (default: "loads.csv").
            bcs_csv (str): Path to the CSV file containing boundary conditions data (default: "bcs.csv").
        """

        # Step 1: Read materials from CSV
        self.read_materials(materials_csv)
        print(f"Materials loaded from '{materials_csv}'.")

        # Step 2: Read real constants from CSV
        self.read_real_constants(real_constants_csv)
        print(f"Real constants loaded from '{real_constants_csv}'.")

        # Step 3: Read loads from CSV
        self.read_csv_loads(loads_csv)
        print(f"Loads loaded from '{loads_csv}' and applied to the mesh.")

        # Step 4: Read boundary conditions from CSV
        self.read_csv_bcs(bcs_csv)
        print(f"Boundary conditions loaded from '{bcs_csv}' and applied to the mesh.")

        # Step 5: Read solver settings from JSON
        self.read_solver_settings_from_json(solver_settings_json)
        print(f"Solver settings loaded from '{solver_settings_json}'.")

    def write_multiple_load_steps(
        self,
        time_data_csv="time_data.csv",
        solver_settings_json="solver_settings.json",
        materials_csv="materials.csv",
        real_constants_csv="real_constants.csv",
        loads_csv="loads_1.csv",
        bcs_csv="bcs_1.csv",
        output_folder=None,
    ):
        """
        Reads time_data.csv and writes multiple load steps based on its rows.
        Dynamically reads loads_1.csv, loads_2.csv, etc., and bcs_1.csv, bcs_2.csv, etc.

        Parameters:
            time_data_csv (str): Path to the CSV file containing time parameters (default: "time_data.csv").
            solver_settings_json (str): Path to the JSON file containing solver settings (default: "solver_settings.json").
            materials_csv (str): Path to the CSV file containing materials data (default: "materials.csv").
            real_constants_csv (str): Path to the CSV file containing real constants data (default: "real_constants.csv").
            loads_csv (str): Path to the CSV file containing initial loads data (default: "loads_1.csv").
            bcs_csv (str): Path to the CSV file containing initial boundary conditions data (default: "bcs_1.csv").
            output_folder (str): Path to the folder where load step binary files will be saved (default: "load_steps").
        """
        # Ensure the output folder exists
        if output_folder is None:
            output_folder = self.load_steps_folder
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        time_data = pd.read_csv(time_data_csv)
        # Ensure the required columns are present
        required_columns = self.time_parameters.keys()
        if not all(column in time_data.columns for column in required_columns):
            raise ValueError(
                f"Time data file must contain the following columns: {required_columns}"
            )

        # Step 1: Read solver settings
        self.read_solver_settings_from_json(solver_settings_json)
        print(f"Solver settings loaded from '{solver_settings_json}'.")

        # Step 2: Read materials from CSV
        self.read_materials(materials_csv)
        print(f"Materials loaded from '{materials_csv}'.")

        # Step 3: Read real constants from CSV
        self.read_real_constants(real_constants_csv)
        print(f"Real constants loaded from '{real_constants_csv}'.")
        # Step 2: Read time_data.csv
        if not os.path.exists(time_data_csv):
            raise FileNotFoundError(f"Time data file '{time_data_csv}' not found.")

        # Step 4: Loop over each row in the time_data DataFrame
        for index, row in time_data.iterrows():
            # Update time parameters with the current row
            self.time_parameters.update(row.to_dict())

            # Dynamically determine the file names for loads and boundary conditions
            if index == 0:
                current_loads_csv = loads_csv
                current_bcs_csv = bcs_csv
            else:
                current_loads_csv = loads_csv.replace("1", str(index + 1))
                current_bcs_csv = bcs_csv.replace("1", str(index + 1))

            # Check if the files exist
            if not os.path.exists(current_loads_csv):
                raise FileNotFoundError(
                    f"Loads file '{current_loads_csv}' not found for step {index + 1}."
                )
            if not os.path.exists(current_bcs_csv):
                raise FileNotFoundError(
                    f"Boundary condition file '{current_bcs_csv}' not found for step {index + 1}."
                )

            # Read and apply loads and boundary conditions
            self.read_csv_loads(current_loads_csv)
            self.read_csv_bcs(current_bcs_csv)

            # Generate the output file name for the current load step
            output_file = os.path.join(output_folder, f"load_step_{index + 1}.bin")

            # Write the current load step to a binary file
            self.write_load_step(output_file)

            print(f"Load step {index + 1} written to '{output_file}'.")

# Only execute this block if the script is run directly
if __name__ == "__main__":
    # Example usage
    apec = APEC_M()
    # print("APEC instance created.")
