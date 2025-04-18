import pandas as pd
from abc import ABC, abstractmethod
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import struct
import os


class APEC:
    """
    APEC class to manage VTK unstructured grids and element types.
    """

    def __init__(self):
        """
        Initialize the APEC object without requiring a file path.
        """
        self.csv_file_path = None
        self.element_types_collector = None
        self.mesh = vtk.vtkUnstructuredGrid()
        self.jobname = None
        self.folderpath = None

    def read_element_types(self, csv_file_path):
        """
        Load element types from the CSV file into the ElementTypesCollector.

        Parameters:
            csv_file_path (str): Path to the CSV file for initializing the ElementTypesCollector.
        """
        self.csv_file_path = csv_file_path
        self.element_types_collector = ElementTypesCollector(self.csv_file_path)

    def read_mesh(self, vtk_file_path):
        """
        Reads a VTK unstructured grid from a file in XML format.

        Parameters:
            vtk_file_path (str): Path to the VTK file in XML format.
        """
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file_path)
        reader.Update()
        self.mesh = reader.GetOutput()

        self.jobname = os.path.splitext(os.path.basename(vtk_file_path))[0]
        self.folderpath = os.path.dirname(vtk_file_path)

        print(f"VTK unstructured grid loaded from '{vtk_file_path}'.")
        # print a summary of the mesh
        print(f"Number of points: {self.mesh.GetNumberOfPoints()}")
        print(f"Number of cells: {self.mesh.GetNumberOfCells()}")
        print(f"Number of point data arrays: {self.mesh.GetPointData().GetNumberOfArrays()}")
        print(f"Number of cell data arrays: {self.mesh.GetCellData().GetNumberOfArrays()}")
        print(f"Number of field data arrays: {self.mesh.GetFieldData().GetNumberOfArrays()}")

    # def read_pipeline(self):
    #     # Get the active object in the pipeline
    #     self.mesh = GetActiveSource()

    #     # Check if an active object exists
    #     if self.pipeline_object is None:
    #         raise RuntimeError("No active object in the pipeline. Please select an object.")

    def write_loads_to_binary(self, array_name):
        """
        Optimized function to write VTK cell data to a binary file using NumPy masking.
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

        # Create a mask to filter out NaN values
        if num_components == 1:
            # Scalar array: directly check for NaN
            mask = ~np.isnan(load_array_np)
        else:
            # Multi-component array: check for NaN across all components
            mask = ~np.isnan(load_array_np).any(axis=1)

        # Get the filtered cell IDs and values
        filtered_cell_ids = np.arange(num_cells)[mask]
        filtered_values = load_array_np[mask]

        # Open the output file in binary mode
        output_file = os.path.join(self.folderpath, f"{self.jobname}_loads_{array_name}.bin")
        with open(output_file, "wb") as file:
            # Write the length of the load name
            load_name_length = len(array_name)
            file.write(struct.pack("Q", load_name_length))  # 'Q' is for unsigned long long (size_t)

            # Write the load name
            file.write(array_name.encode("utf-8"))

            # Write the number of components
            file.write(struct.pack("I", num_components))  # 'I' is for unsigned int

            # Write the filtered data
            for cell_id, values in zip(filtered_cell_ids, filtered_values):
                for component, value in enumerate(values if num_components > 1 else [values]):
                    # Write the element number (cell ID)
                    file.write(struct.pack("I", cell_id))  # 'I' is for uint32_t (dof_id_type)

                    # Write the component number
                    file.write(struct.pack("I", component))  # 'I' is for unsigned int

                    # Write the load value
                    file.write(struct.pack("d", value))  # 'd' is for double (Real)


    def write_boundary_conditions_to_binary(self, array_name):
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
        output_file = os.path.join(self.folderpath, f"{self.jobname}_bcs.bin")

        # Open the output file in binary mode
        with open(output_file, "wb") as file:
            # Iterate over components
            for component in range(num_components):
                # Generate the label for this component (e.g., "A_x", "A_y", "A_z")
                component_label = f"{array_name}_{chr(120 + component)}"  # 'x', 'y', 'z', etc.

                # Create a mask to filter out NaN values for this component
                component_values = bc_array_np[:, component]
                mask = ~np.isnan(component_values)

                # Get the filtered point IDs and values
                filtered_point_ids = np.arange(num_points)[mask]
                filtered_values = component_values[mask]

                # Write the data for this component
                for point_id, value in zip(filtered_point_ids, filtered_values):
                    # Write the point ID (node number)
                    file.write(struct.pack("I", point_id))  # 'I' is for uint32_t (dof_id_type)

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
        if self.element_types_collector is None:
            raise ValueError("Element types have not been loaded. Call 'load_element_types' first.")

        element_types = self.element_types_collector.get_element_types()

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

    def replot(self):
        """
        Replot the mesh in the active view.
        """
        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        print("Replotting the mesh...")
        producer = TrivialProducer()
        producer.GetClientSideObject().SetOutput(self.mesh)
        Show(producer)
        Render()

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
        return {"A_z": 1}

    def getElementLoads(self):
        return {"Js": 1}


class PLANE233Type(ElementType):
    def getVariables(self):
        return {"A_z": 1}

    def getElementLoads(self):
        return {"Js": 1}


class ElementTypesCollector:
    """
    Class to collect and manage element types from a CSV file.
    """

    def __init__(self, csv_file_path):
        self.element_types_map = {}
        self.read_element_types(csv_file_path)

    def read_element_types(self, csv_file_path):
        df = pd.read_csv(csv_file_path)

        for _, row in df.iterrows():
            element_type = row["element_type"]
            if element_type == "SOLID97":
                self.element_types_map[row["ansys_element_type"]] = SOLID97Type()
            elif element_type == "SOLID236":
                self.element_types_map[row["ansys_element_type"]] = SOLID236Type()
            elif element_type == "PLANE13":
                self.element_types_map[row["ansys_element_type"]] = PLANE13Type()
            elif element_type == "PLANE233":
                self.element_types_map[row["ansys_element_type"]] = PLANE233Type()
            else:
                raise ValueError(f"Unknown element type: {element_type}")

    def get_element_types(self):
        return self.element_types_map

    def get_all_loads(self):
        all_loads = {}
        for element_type in self.element_types_map.values():
            all_loads.update(element_type.getElementLoads())
        return all_loads

apec = APEC()
