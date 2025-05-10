import pandas as pd
from abc import ABC, abstractmethod
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import struct
import os
import sys
import re
import subprocess

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

try:
    # Relative import for package execution
    from .constants import *
except ImportError:
    # Absolute import for standalone execution
    from constants import *


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

        self.element_types = pd.DataFrame(columns=ELEMENT_TYPES_COLUMNS)
        self.materials = pd.DataFrame(columns=MATERIALS_COLUMNS)
        self.real_constants = pd.DataFrame(columns=REAL_CONSTANTS_COLUMNS)
        self.time_parameters = DEFAULT_TIME_PARAMETERS
        self.apec_parameters = {}
        self.solver_settings = DEFAULT_SOLVER_SETTINGS
        self.transient_analysis = pd.DataFrame(columns=TRANSIENT_ANALYSIS_COLUMNS)

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
            if pd.api.types.is_integer_dtype(df[column_name]):
                vtk_array = vtk.vtkIntArray()
            elif pd.api.types.is_float_dtype(df[column_name]):
                vtk_array = vtk.vtkDoubleArray()
            else:
                vtk_array = vtk.vtkStringArray()

            vtk_array.SetName(column_name)

            # Add data to the vtk array
            for value in df[column_name]:
                if isinstance(vtk_array, vtk.vtkIntArray):
                    vtk_array.InsertNextValue(int(value))
                elif isinstance(vtk_array, vtk.vtkDoubleArray):
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

        # Create a folder named 'loads_steps' in the current working directory
        self.load_steps_folder = os.path.join(os.getcwd(), "load_steps")
        if not os.path.exists(self.load_steps_folder):
            os.makedirs(self.load_steps_folder)

        # Check if the arrays already exist in the mesh
        cell_data = self.mesh.GetCellData()
        point_data = self.mesh.GetPointData()

        # Add a cell data array "Js" with 3 components, filled with np.nan, if it doesn't exist
        if cell_data.GetArray("Js") is None:
            num_cells = self.mesh.GetNumberOfCells()
            Js = np.full((num_cells, 3), np.nan)
            cell_data.AddArray(
                vtknp.numpy_to_vtk(Js, deep=True, array_type=vtk.VTK_DOUBLE)
            )
            cell_data.GetArray("Js").SetName("Js")

        # Add a point data array "A" with 3 components, filled with np.nan, if it doesn't exist
        if point_data.GetArray("A") is None:
            num_points = self.mesh.GetNumberOfPoints()
            A = np.full((num_points, 3), np.nan)
            point_data.AddArray(
                vtknp.numpy_to_vtk(A, deep=True, array_type=vtk.VTK_DOUBLE)
            )
            point_data.GetArray("A").SetName("A")

        # Add a point data array "VOLT" with 1 component, filled with np.nan, if it doesn't exist
        if point_data.GetArray("VOLT") is None:
            num_points = self.mesh.GetNumberOfPoints()
            VOLT = np.full(num_points, np.nan)
            point_data.AddArray(
                vtknp.numpy_to_vtk(VOLT, deep=True, array_type=vtk.VTK_DOUBLE)
            )
            point_data.GetArray("VOLT").SetName("VOLT")

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

    def validate_dataframe_with_vtk_data(self, vtk_data, df, array_name):
        """
        Validates that the first column of a DataFrame has unique values,
        all values have a correspondence in the specified VTK data array,
        and there are no undefined values in the VTK data array.

        Parameters:
            vtk_data (vtk.vtkDataSetAttributes): The VTK cell or point data to validate against.
            df (pd.DataFrame): The DataFrame to validate.
            array_name (str): The name of the VTK data array to check.

        Raises:
            ValueError: If the DataFrame is empty, the first column is not unique,
                        the VTK data array is missing, or the checks fail.
        """
        # Check if the DataFrame is empty
        if df.empty:
            raise ValueError("The DataFrame is empty.")

        # Validate that the first column has unique values
        first_column = df.iloc[:, 0]
        if not first_column.is_unique:
            raise ValueError(
                "The first column of the DataFrame must have unique values."
            )

        # Retrieve the VTK data array
        vtk_array = vtk_data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"The '{array_name}' array is missing in the VTK data.")

        # Convert the VTK array to a NumPy array
        vtk_array_np = vtknp.vtk_to_numpy(vtk_array)

        # Check if all values in the first column have a correspondence in the VTK array
        if not set(first_column).issubset(set(vtk_array_np)):
            raise ValueError(
                f"Some values in the first column of the DataFrame do not have a correspondence in the '{array_name}' array."
            )

        # Check for undefined values in the VTK array
        undefined_values = set(vtk_array_np) - set(first_column)
        if undefined_values:
            raise ValueError(
                f"The VTK data contains undefined values: {undefined_values}."
            )

        print(f"Validation passed for '{array_name}' array and DataFrame.")

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

        # Validate the element_types DataFrame against the 'ansys_elem_num' array in cell data
        cell_data = self.mesh.GetCellData()
        self.validate_dataframe_with_vtk_data(
            cell_data, self.element_types, "ansys_etype"
        )

        # Convert the element_types DataFrame to a vtkTable
        vtk_table = self.dataframe_to_vtk_table(self.element_types)

        # Attach the vtkTable as field data to the mesh
        self.add_table_to_mesh(vtk_table)

        # Define the output file name
        output_file = os.path.join(self.folderpath, f"{self.jobname}.vtu")

        # Write the unstructured grid to the .vtu file
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(self.mesh)
        writer.SetDataModeToBinary()  # Set the writer to binary mode
        writer.Write()

        print(f"Mesh saved to '{output_file}'.")

    def read_excel_df(self, excel_file, sheet_name):
        """
        Reads a specific sheet from an Excel file into a pandas DataFrame and assigns it to the corresponding attribute.
        Ensures that the columns in the sheet match the default columns defined in DEFAULT_COLUMNS_ATTRIBUTES.

        Parameters:
            excel_file (str): The path to the Excel file.
            sheet_name (str): The name of the sheet to read.

        Raises:
            ValueError: If the sheet name is not in DEFAULT_COLUMNS_ATTRIBUTES or if the columns do not match.
            FileNotFoundError: If the Excel file does not exist.
        """
        # Check if the sheet name is valid
        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES.keys())}."
            )

        # Read the Excel file
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")

        # Check if the sheet exists in the Excel file
        try:
            # Read the specified sheet into a DataFrame
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            # Raise a warning if the sheet does not exist
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
            return

        # Get the expected columns for the sheet
        expected_columns = DEFAULT_COLUMNS_ATTRIBUTES[sheet_name]

        # Check if the columns in the DataFrame match the expected columns
        if not all(column in expected_columns for column in df.columns):
            raise ValueError(
                f"Columns in sheet '{sheet_name}' do not match the expected columns. "
                f"Expected: {expected_columns}, Found: {list(df.columns)}"
            )

        # Assign the DataFrame to the corresponding attribute
        setattr(self, sheet_name, df)

        print(
            f"Data from sheet '{sheet_name}' successfully loaded into attribute '{sheet_name}'."
        )

    def read_excel_dict(self, excel_file, sheet_name):
        """
        Reads a specific sheet from an Excel file into a dictionary attribute.
        Ensures that the keys in the 'SETTING' column match the default keys defined in DEFAULT_COLUMNS_ATTRIBUTES_DICT.

        Parameters:
            excel_file (str): The path to the Excel file.
            sheet_name (str): The name of the sheet to read.

        Raises:
            ValueError: If the sheet name is not in DEFAULT_COLUMNS_ATTRIBUTES_DICT or if the keys do not match.
            FileNotFoundError: If the Excel file does not exist.
        """
        # Check if the sheet name is valid
        if sheet_name not in DEFAULT_COLUMNS_ATTRIBUTES.keys():
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES.keys())}."
            )

        # Read the Excel file
        if not os.path.exists(excel_file):
            raise FileNotFoundError(f"Excel file '{excel_file}' not found.")

        # Check if the sheet exists in the Excel file
        try:
            # Read the specified sheet into a DataFrame
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            # Raise a warning if the sheet does not exist
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
            return

        # Ensure the required columns are present
        required_columns = ["SETTING", "VALUE"]
        if not all(column in df.columns for column in required_columns):
            raise ValueError(
                f"Sheet '{sheet_name}' must contain the following columns: {required_columns}."
            )

        # Get the expected keys for the dictionary
        expected_keys = DEFAULT_COLUMNS_ATTRIBUTES[sheet_name].keys()

        # Validate the keys in the 'SETTING' column
        invalid_keys = [key for key in df["SETTING"] if key not in expected_keys]
        if invalid_keys:
            raise ValueError(
                f"Invalid keys found in 'SETTING' column of sheet '{sheet_name}': {invalid_keys}. "
                f"Expected keys: {list(expected_keys)}."
            )

        # Populate the dictionary
        target_dict = getattr(self, sheet_name)
        for _, row in df.iterrows():
            setting = row["SETTING"]
            value = row["VALUE"]

            # Convert the value to the appropriate type based on the default dictionary
            if isinstance(target_dict[setting], int):
                value = int(value)
            elif isinstance(target_dict[setting], float):
                value = float(value)
            elif isinstance(target_dict[setting], str):
                value = str(value)
            else:
                raise ValueError(f"Unsupported value type for key '{setting}'.")

            # Assign the value to the dictionary
            target_dict[setting] = value

        print(
            f"Data from sheet '{sheet_name}' successfully loaded into dictionary '{sheet_name}'."
        )

    def write_excel_loads_bcs(self, output_excel_file, data_type="loads"):
        """
        Writes the loads or boundary conditions to an Excel file.

        Parameters:
            output_excel_file (str): The path to the output Excel file.
            data_type (str): The type of data to write ("loads" or "boundary_conditions").
                            Defaults to "loads".

        Raises:
            ValueError: If the data_type is invalid.
        """
        # Determine the VTK data and mapping based on the data_type
        if data_type == "loads":
            vtk_data = self.mesh.GetCellData()
            mapping = VARIABLES_LOADS_NUMBERS
        elif data_type == "boundary_conditions":
            vtk_data = self.mesh.GetPointData()
            mapping = VARIABLES_GLOBAL_NUMBERS
        else:
            raise ValueError(
                "Invalid data_type. Must be 'loads' or 'boundary_conditions'."
            )

        # Retrieve the final_data NumPy array
        final_data = self.create_numpy_array_for_binary(vtk_data, mapping)

        # Convert the NumPy array to a pandas DataFrame
        df = pd.DataFrame(final_data)

        # Replace the "load_bc_number" column with the corresponding "label"
        reverse_mapping = {v: k for k, v in mapping.items()}  # Reverse the mapping
        df["label"] = df["load_bc_number"].map(reverse_mapping)  # Map numbers to labels

        # Reorder columns to make "label" the second column and remove "load_bc_number"
        df = df[["vtk_id", "label", "component", "value"]]

        # Check if the Excel file already exists
        if os.path.exists(output_excel_file):
            # Open the existing Excel file and add a new sheet
            with pd.ExcelWriter(
                output_excel_file, engine="openpyxl", mode="a"
            ) as writer:
                df.to_excel(writer, sheet_name=data_type, index=False)
        else:
            # Create a new Excel file and write the DataFrame to it
            with pd.ExcelWriter(output_excel_file, engine="openpyxl") as writer:
                df.to_excel(writer, sheet_name=data_type, index=False)
        print(f"{data_type.capitalize()} written to '{output_excel_file}'.")

    def write_settings_to_binary(self, settings, binary_file):
        """
        Writes a settings dictionary to an open binary file stream.

        Parameters:
            settings (dict): The settings dictionary to write.
            binary_file (file object): An open binary file stream.

        Raises:
            ValueError: If a value in the settings dictionary has an unsupported type.
        """
        # Write the number of settings
        binary_file.write(struct.pack("i", len(settings)))  # Number of settings (int)

        for key, value in settings.items():
            # Write the size of the string label and the label itself
            binary_file.write(struct.pack("i", len(key)))  # Length of the key (int)
            binary_file.write(key.encode("utf-8"))  # Key as a UTF-8 encoded string

            # Determine the type of the value and write it
            if isinstance(value, int):
                binary_file.write(struct.pack("i", 0))  # Type: int (0)
                binary_file.write(struct.pack("i", value))  # Value as int
            elif isinstance(value, float):
                binary_file.write(struct.pack("i", 1))  # Type: double (1)
                binary_file.write(struct.pack("d", value))  # Value as double
            elif isinstance(value, str):
                binary_file.write(struct.pack("i", 2))  # Type: string (2)
                binary_file.write(struct.pack("i", len(value)))  # Length of the string
                binary_file.write(
                    value.encode("utf-8")
                )  # Value as UTF-8 encoded string
            else:
                raise ValueError(f"Unsupported value type: {type(value)}")

    def write_dataframe_to_binary(self, df, binary_file, key):
        """
        Writes a pandas DataFrame to an open binary file stream using predefined types for each column.

        The binary file format includes:
            - Number of rows (unsigned int)
            - Number of columns (unsigned int)
            - Column metadata:
                - Column name (string)
                - Column type (int: 0=int, 1=double, 2=string)
            - Data values (row by row):
                - Each column is written based on its predefined type in `DEFAULT_COLUMNS_ATTRIBUTES_TYPES`.

        Parameters:
            df (pd.DataFrame): The DataFrame to write.
            binary_file (file object): An open binary file stream.
            key (str): The key to retrieve the column types from DEFAULT_COLUMNS_ATTRIBUTES_TYPES.

        Raises:
            ValueError: If the key is invalid or if a column in the DataFrame does not have a predefined type.
        """
        # Validate the DataFrame if it is self.materials or self.real_constants
        if df is self.materials:
            self.validate_dataframe_with_vtk_data(
                self.mesh.GetCellData(), df, "ansys_material_type"
            )
        elif df is self.real_constants:
            self.validate_dataframe_with_vtk_data(
                self.mesh.GetCellData(), df, "ansys_real_constant"
            )
        # Retrieve the column types dictionary
        if key not in DEFAULT_COLUMNS_ATTRIBUTES_TYPES:
            raise ValueError(
                f"Invalid key '{key}'. Must be one of {list(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys())}."
            )
        column_types = DEFAULT_COLUMNS_ATTRIBUTES_TYPES[key]

        # Ensure all columns in the DataFrame have a corresponding type in column_types
        for column in df.columns:
            if column not in column_types:
                raise ValueError(
                    f"Column '{column}' does not have a predefined type in column_types."
                )

        # Write the number of rows and columns
        num_rows = len(df)
        num_columns = len(df.columns)
        binary_file.write(struct.pack("I", num_rows))  # Number of rows (unsigned int)
        binary_file.write(
            struct.pack("I", num_columns)
        )  # Number of columns (unsigned int)

        # Write column metadata (name and type)
        for column in df.columns:
            column_name = column
            column_type = column_types[column]

            # Write the column name
            binary_file.write(
                struct.pack("I", len(column_name))
            )  # Length of the column name
            binary_file.write(
                column_name.encode("utf-8")
            )  # Column name as UTF-8 encoded string

            # Write the column type (0=int, 1=double, 2=string)
            if column_type == int:
                binary_file.write(struct.pack("i", 0))  # Type: int
            elif column_type == float:
                binary_file.write(struct.pack("i", 1))  # Type: double
            elif column_type == str:
                binary_file.write(struct.pack("i", 2))  # Type: string
            else:
                raise ValueError(
                    f"Unsupported column type '{column_type}' for column '{column}'."
                )

        # Write the data row by row
        for _, row in df.iterrows():
            for column in df.columns:
                value = row[column]
                column_type = column_types[column]

                # Write the value based on its predefined type
                if column_type == int:
                    binary_file.write(struct.pack("i", int(value)))  # Signed integer
                elif column_type == float:
                    binary_file.write(
                        struct.pack("d", float(value))
                    )  # Double (float64)
                elif column_type == str:
                    encoded_value = value.encode("utf-8")
                    binary_file.write(
                        struct.pack("I", len(encoded_value))
                    )  # String length
                    binary_file.write(encoded_value)  # String bytes
                else:
                    raise ValueError(
                        f"Unsupported column type '{column_type}' for column '{column}'."
                    )

    def create_numpy_array_for_binary(self, cell_point_data, loads_bc_labels):
        """
        Creates a NumPy structured array for loads from a VTK cell data array.

        The structured array contains:
            - cell_id (unsigned int)
            - load_number (unsigned int, from VARIABLES_LOADS_NUMBERS mapping)
            - component (unsigned int)
            - value (double)

        Returns:
            np.ndarray: A structured NumPy array containing the loads or bcs.
        """

        total_rows = 0  # To keep track of the total number of rows
        all_stacked_data = []  # Temporary storage for structured arrays

        # Define the structured array data type
        dtype = np.dtype(
            [
                ("vtk_id", np.uint32),
                ("load_bc_number", np.uint32),
                ("component", np.uint32),
                ("value", np.float64),
            ]
        )

        # Iterate over each load label
        for load_bc_label in loads_bc_labels.keys():
            load_bc_array = cell_point_data.GetArray(load_bc_label)
            if load_bc_array is None:
                continue  # Skip if the array does not exist in the VTK file

            # Convert the VTK array to a NumPy array
            load_bc_array_np = vtknp.vtk_to_numpy(load_bc_array)

            # Get the number of components and cells
            num_components = load_bc_array.GetNumberOfComponents()
            num_cells = cell_point_data.GetArray(0).GetNumberOfTuples()

            # Iterate over components to filter data
            for component in range(num_components):
                # Extract the values for the current component
                if load_bc_array_np.ndim == 1:
                    # For 1D arrays, use the array directly
                    component_values = load_bc_array_np
                else:
                    # For 2D arrays, access the specific component
                    component_values = load_bc_array_np[:, component]

                # Create a mask to filter out NaN values
                mask = ~np.isnan(component_values)

                # Get the filtered cell IDs and values
                filtered_cell_points_ids = np.arange(num_cells, dtype=np.uint32)[mask]
                filtered_values = component_values[mask]

                # Create arrays for load number and component
                load_bc_number_array = np.full(
                    len(filtered_cell_points_ids),
                    loads_bc_labels[load_bc_label],
                    dtype=np.uint32,
                )
                component_array = np.full(
                    len(filtered_cell_points_ids), component, dtype=np.uint32
                )

                # Create a structured array for the current component
                structured_data = np.zeros(len(filtered_cell_points_ids), dtype=dtype)
                structured_data["vtk_id"] = filtered_cell_points_ids
                structured_data["load_bc_number"] = load_bc_number_array
                structured_data["component"] = component_array
                structured_data["value"] = filtered_values

                # Append the structured array to the list
                all_stacked_data.append(structured_data)

                # Update the total number of rows
                total_rows += len(filtered_cell_points_ids)

        # Concatenate all structured arrays into a single array
        if all_stacked_data:
            final_data = np.concatenate(all_stacked_data)
            return final_data

        # Return an empty array if no data was found
        return np.zeros(0, dtype=dtype)

    def get_loads_numpy_for_binary(self):
        """
        Retrieves a structured NumPy array for loads using the cell data of the VTK mesh.

        Returns:
            np.ndarray: A structured NumPy array containing the loads.
        """
        # Use the cell data of the VTK mesh and the load labels from VARIABLES_LOADS_NUMBERS
        cell_data = self.mesh.GetCellData()

        # Call create_numpy_array_for_binary to generate the structured array
        return self.create_numpy_array_for_binary(cell_data, VARIABLES_LOADS_NUMBERS)

    def get_bcs_numpy_for_binary(self):
        """
        Retrieves a structured NumPy array for boundary conditions using the point data of the VTK mesh.

        Returns:
            np.ndarray: A structured NumPy array containing the boundary conditions.
        """
        # Use the point data of the VTK mesh and the boundary condition labels from VARIABLES_GLOBAL_NUMBERS
        point_data = self.mesh.GetPointData()

        # Call create_numpy_array_for_binary to generate the structured array
        return self.create_numpy_array_for_binary(point_data, VARIABLES_GLOBAL_NUMBERS)

    def write_numpy_array_to_binary(self, binary_file, numpy_array):
        """
        Writes a NumPy structured array to an open binary file stream.

        Parameters:
            binary_file (file object): An open binary file stream.
            numpy_array (np.ndarray): The NumPy structured array to write.

        Raises:
            ValueError: If the input array is not a NumPy structured array.
        """
        if not isinstance(numpy_array, np.ndarray):
            raise ValueError("Input must be a NumPy array.")

        # Write the total number of rows to the file
        binary_file.write(struct.pack("I", len(numpy_array)))  # 'I' is for unsigned int

        # Write the array data to the file
        binary_file.write(numpy_array.tobytes())

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

    def map_data_to_vtk(self, data, label):
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
        # Get the array with the specified label
        vtk_array = data.GetArray(label)
        if vtk_array is None:
            raise ValueError(f"Array '{label}' not found in the provided VTK data.")

        # Convert the VTK array to a NumPy array
        ansys_numbers = vtknp.vtk_to_numpy(vtk_array)

        # Create a 2D NumPy array with VTK IDs and ANSYS numbers
        vtk_ids = np.arange(len(ansys_numbers))  # VTK IDs are sequential from 0
        mapping_array = np.column_stack((vtk_ids, ansys_numbers))

        return mapping_array

    def replace_ansys_with_vtk_ids(self, df, mapping_array, label):
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
        ansys_numbers = df[label].to_numpy()

        # Extract the VTK IDs using NumPy indexing
        vtk_ids = mapping_array[
            :, 0
        ]  # VTK IDs are in the first column of the mapping array
        ansys_mapping = mapping_array[:, 1]  # ANSYS numbers are in the second column

        vtk_id_array = vtk_ids[np.searchsorted(ansys_mapping, ansys_numbers)]

        # Replace the column in the DataFrame and rename it to 'vtk_element_id'
        df.insert(
            0, "vtk_id", vtk_id_array
        )  # Insert 'vtk_element_id' as the first column
        df = df.drop(columns=[label])

        return df

    def update_vtk_grid(self, data, df):
        """
        Updates the VTK grid's cell data arrays using the DataFrame resulting from `replace_ansys_with_vtk_ids`.
        The function assigns values to the appropriate components at the corresponding VTK IDs.

        Args:
            df (pd.DataFrame): The DataFrame containing columns 'vtk_id', 'label', 'component', and 'value'.
        """
        # Get the distinct load labels
        loads_bc_labels = df["label"].unique()

        for load_bc_label in loads_bc_labels:
            # Filter the DataFrame for the current load label
            load_bc_df = df[df["label"] == load_bc_label]

            # Get the corresponding VTK array
            vtk_array = data.GetArray(load_bc_label)
            if vtk_array is None:
                raise ValueError(f"Array '{load_bc_label}' not found in the VTK file.")

            # Convert the VTK array to a NumPy array
            vtk_array_np = vtknp.vtk_to_numpy(vtk_array)

            # Check if the NumPy array is 1D
            if vtk_array_np.ndim == 1:
                # Iterate over distinct component values
                for component in load_bc_df["component"].unique():
                    # Filter rows for the current component
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()

                    # Assign values directly to the 1D array
                    vtk_array_np[vtk_ids] = values
            else:
                # Iterate over distinct component values
                for component in load_bc_df["component"].unique():
                    # Filter rows for the current component
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()

                    # Assign values to the appropriate indices and component
                    vtk_array_np[vtk_ids, component] = values

            # Convert the NumPy array back to a VTK array
            vtk_array.Modified()

    def update_vtk_grid_from_excel(self, excel_file, sheet_name, data_type):
        """
        Reads a sheet from an Excel file and updates the VTK grid using the data.

        Parameters:
            excel_file (str): The path to the Excel file.
            sheet_name (str): The name of the sheet to read.
            data_type (str): The type of VTK data to update ('cell' or 'point').

        Raises:
            ValueError: If the data_type is invalid or the sheet does not exist.
        """
        # Check if the sheet exists in the Excel file
        try:
            df = pd.read_excel(excel_file, sheet_name=sheet_name, engine="openpyxl")
        except ValueError:
            print(f"Warning: Sheet '{sheet_name}' not found in '{excel_file}'.")
            return

        # Determine the VTK data type to update
        if data_type == "cell":
            vtk_data = self.mesh.GetCellData()
        elif data_type == "point":
            vtk_data = self.mesh.GetPointData()
        else:
            raise ValueError(
                f"Invalid data_type '{data_type}'. Must be 'cell' or 'point'."
            )

        # Update the VTK grid with the DataFrame
        self.update_vtk_grid(vtk_data, df)
        print(f"VTK grid updated with data from sheet '{sheet_name}'.")

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
        mapping_array = self.map_data_to_vtk(self.mesh.GetCellData(), "ansys_elem_num")

        # Step 3: Replace ANSYS element numbers with VTK IDs in the DataFrame
        df = self.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_elem_num")

        # Step 4: Update the VTK grid with the modified DataFrame
        self.update_vtk_grid(self.mesh.GetCellData(), df)

        print(f"Grid updated with element loads from '{file_path}'.")

    def convert_bfelist(self, bfelist_file, output_excel_file):
        """
        Reads a BFELIST file, updates the VTK grid with the loads, and writes the loads to an Excel file.

        Parameters:
            bfelist_file (str): Path to the BFELIST file.
            output_excel_file (str): Path to the output Excel file.
        """
        # Check if the BFELIST file exists
        if not os.path.exists(bfelist_file):
            raise FileNotFoundError(f"BFELIST file '{bfelist_file}' not found.")

        # Read the BFELIST file and update the VTK grid
        self.read_bfelist(bfelist_file)

        # Write the loads to the Excel file
        self.write_excel_loads_bcs(output_excel_file, data_type="loads")

        print(f"Loads from '{bfelist_file}' written to '{output_excel_file}'.")

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
        mapping_array = self.map_data_to_vtk(self.mesh.GetPointData(), "ansys_node_num")

        # Step 3: Replace ANSYS node numbers with VTK point IDs in the DataFrame
        df = self.replace_ansys_with_vtk_ids(df, mapping_array, "ansys_node_num")

        # Step 4: Update the VTK grid with the modified DataFrame
        self.update_vtk_grid(self.mesh.GetPointData(), df)

        print(f"Grid updated with node data from '{file_path}'.")

    def convert_dlist(self, dlist_file, output_excel_file):
        """
        Reads a DLIST file, updates the VTK grid with the node data, and writes the data to an Excel file.

        Parameters:
            dlist_file (str): Path to the DLIST file.
            output_excel_file (str): Path to the output Excel file.
        """
        # Check if the DLIST file exists
        if not os.path.exists(dlist_file):
            raise FileNotFoundError(f"DLIST file '{dlist_file}' not found.")

        # Read the DLIST file and update the VTK grid
        self.read_dlist(dlist_file)

        # Write the node data to the Excel file
        self.write_excel_loads_bcs(output_excel_file, data_type="boundary_conditions")

        print(f"Node data from '{dlist_file}' written to '{output_excel_file}'.")

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

        # Check if the sheet name is valid
        if sheet_name not in SHEET_PROCESSING_MAP:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(SHEET_PROCESSING_MAP.keys())}."
            )

        # Get the processing method or method with additional arguments
        processing_info = SHEET_PROCESSING_MAP[sheet_name]

        # Call the appropriate method
        if isinstance(processing_info, str):
            # Call read_excel_df or read_excel_dict
            getattr(self, processing_info)(excel_file, sheet_name)
        elif isinstance(processing_info, tuple):
            # Call update_vtk_grid_from_excel with additional arguments
            method_name, data_type = processing_info
            getattr(self, method_name)(excel_file, sheet_name, data_type)

        print(f"Read sheet '{sheet_name}' from '{excel_file}'.")

    def read_excel(
        self,
        excel_file_path,
        data_to_read=[
            "element_types",
            "materials",
            "real_constants",
            "solver_settings",
            "time_parameters",
            "loads",
            "boundary_conditions",
        ],
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

    def write_load_step(self, output_file, labels_to_write=None):
        """
        Writes the current load step data to a single binary file based on the provided labels.

        Parameters:
            output_file (str): The path to the output binary file.
            labels_to_write (list): A list of labels specifying what to write. Defaults to all supported labels.
        """

        # Default labels to write if none are provided
        if labels_to_write is None:
            labels_to_write = list(WRITING_FUNCTIONS_MAP.keys())

        # Open the file in write mode to overwrite if it already exists
        with open(output_file, "wb") as binary_file:
            for label in labels_to_write:
                try:
                    # Call the appropriate function based on the label
                    self._write_label_to_binary(binary_file, label)
                    print(f"{label.replace('_', ' ').capitalize()} written.")
                except ValueError as e:
                    print(f"Warning: {e}")

        print(f"Load step written to '{output_file}'.")

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

        # Call the appropriate function dynamically
        if label in ["time_parameters", "solver_settings"]:
            # Write settings dictionaries
            data = getattr(self, label, None)
            if data is None:
                raise ValueError(f"Data for '{label}' is missing.")
            getattr(self, function_name)(data, binary_file)
        elif label in ["materials", "real_constants"]:
            # Write DataFrames
            data = getattr(self, label, None)
            if data is None or data.empty:
                raise ValueError(f"DataFrame for '{label}' is empty or missing.")
            getattr(self, function_name)(data, binary_file, label)
        elif label in ["loads", "boundary_conditions"]:
            # Write NumPy arrays
            if label == "loads":
                data = self.get_loads_numpy_for_binary()
            elif label == "boundary_conditions":
                data = self.get_bcs_numpy_for_binary()
            if data is None or len(data) == 0:
                raise ValueError(f"NumPy array for '{label}' is empty or missing.")
            getattr(self, function_name)(binary_file, data)

    def generate_transient_analysis(
        self,
        excel_file_path,
    ):

        if not os.path.exists(excel_file_path):
            raise FileNotFoundError(f"Excel file '{excel_file_path}' not found.")

        self.read_excel_sheet(excel_file_path, "transient_analysis")
        print(f"Transient analysis loaded from '{excel_file_path}'.")

        for idx, row in self.transient_analysis.iterrows():
            # Step 1: Update self.time_parameters with the values from the current row
            for column in TRANSIENT_ANALYSIS_COLUMNS:
                if column in row:
                    self.time_parameters[column] = row[column]

            # Check if the file exists before reading
            self.read_excel(
                row["load_step_input_file"],
                [
                    "materials",
                    "real_constants",
                    "solver_settings",
                    "loads",
                    "boundary_conditions",
                ],
            )

            self.write_load_step(
                os.path.join(self.load_steps_folder, f"load_step_{idx}.bin")
            )
            print(f"Load step {idx} written to binary file.")


# Only execute this block if the script is run directly
if __name__ == "__main__":
    # Example usage
    apec = APEC_M()
    # print("APEC instance created.")
