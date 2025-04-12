import pandas as pd
from abc import ABC, abstractmethod
import vtk
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

    def load_element_types(self, csv_file_path):
        """
        Load element types from the CSV file into the ElementTypesCollector.

        Parameters:
            csv_file_path (str): Path to the CSV file for initializing the ElementTypesCollector.
        """
        self.csv_file_path = csv_file_path
        self.element_types_collector = ElementTypesCollector(self.csv_file_path)

    def write_multiple_cell_data_to_file(self, unstructured_grid, input_basename, array_names, num=None):
        """
        Writes multiple cell data arrays from a vtkUnstructuredGrid object to a single binary file.
        """
        output_file_name = f"{input_basename}_loads"
        if num is not None:
            output_file_name += f"_{num}"
        output_file_name += ".bin"

        with open(output_file_name, 'wb') as file:
            for array_name in array_names:
                cell_data = unstructured_grid.GetCellData()
                data_array = cell_data.GetArray(array_name)
                if data_array is None:
                    raise ValueError(f"Array '{array_name}' not found in cell data.")

                name_length = len(array_name)
                file.write(struct.pack('B', name_length))
                file.write(array_name.encode('utf-8'))

                num_components = data_array.GetNumberOfComponents()
                file.write(struct.pack('I', num_components))

                num_elements = data_array.GetNumberOfTuples()
                file.write(struct.pack('I', num_elements))

                for cell_id in range(num_elements):
                    file.write(struct.pack('I', cell_id))
                    values = data_array.GetTuple(cell_id)
                    file.write(struct.pack(f'{num_components}d', *values))

        print(f"All cell data arrays {array_names} written to '{output_file_name}'.")

    def process_loads(self, vtk_file_path, array_names, num=None):
        """
        Reads a VTK unstructured grid from a file and writes its cell data arrays to a single binary file.
        """
        input_basename = os.path.splitext(os.path.basename(vtk_file_path))[0]
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file_path)
        reader.Update()
        unstructured_grid = reader.GetOutput()
        self.write_multiple_cell_data_to_file(unstructured_grid, input_basename, array_names, num)

    def write_multiple_boundary_conditions(self, vtk_grid, output_basename, array_names, num=None):
        """
        Writes multiple point data arrays from a VTK unstructured grid to a single binary file.
        """
        output_file = f"{output_basename}_bc"
        if num is not None:
            output_file += f"_{num}"
        output_file += ".dat"

        with open(output_file, "wb") as file:
            for array_name in array_names:
                point_data = vtk_grid.GetPointData()
                data_array = point_data.GetArray(array_name)
                if data_array is None:
                    raise ValueError(f"Point data array '{array_name}' not found in the VTK grid.")

                num_components = data_array.GetNumberOfComponents()

                for node_id in range(data_array.GetNumberOfTuples()):
                    values = data_array.GetTuple(node_id)

                    for component_idx, value in enumerate(values):
                        if np.isnan(value):
                            continue

                        if num_components > 1:
                            component_label = f"{array_name}_{['x', 'y', 'z'][component_idx]}"
                        else:
                            component_label = array_name

                        file.write(struct.pack("Q", node_id))
                        label = component_label.encode("utf-8")
                        label_length = len(label)
                        file.write(struct.pack("B", label_length))
                        file.write(label)
                        file.write(struct.pack("d", value))

        print(f"Boundary conditions for arrays {array_names} written to {output_file}.")

    def process_boundary_conditions(self, vtk_file_path, array_names, num=None):
        """
        Reads a VTK unstructured grid from a file and writes its point data arrays to a single binary file.
        """
        output_basename = os.path.splitext(os.path.basename(vtk_file_path))[0]
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file_path)
        reader.Update()
        vtk_grid = reader.GetOutput()
        self.write_multiple_boundary_conditions(vtk_grid, output_basename, array_names, num)

    def initialize_mesh(self, grid):
        """
        Initialize a vtkUnstructuredGrid object by adding variables as point data
        and element loads as cell data.
        """
        if self.element_types_collector is None:
            raise ValueError("Element types have not been loaded. Call 'load_element_types' first.")

        element_types = self.element_types_collector.get_element_types()

        for element_id, element in element_types.items():
            variables = element.getVariables()
            for var_name, num_components in variables.items():
                array = vtk.vtkDoubleArray()
                array.SetName(var_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(grid.GetNumberOfPoints())
                grid.GetPointData().AddArray(array)

        for element_id, element in element_types.items():
            element_loads = element.getElementLoads()
            for load_name, num_components in element_loads.items():
                array = vtk.vtkDoubleArray()
                array.SetName(load_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(grid.GetNumberOfCells())
                grid.GetCellData().AddArray(array)


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


if __name__ == "__main__":
    apec = APEC()

    # Load element types
    csv_file_path = "/home/bitteal/PlasmaArc/simulations/3D_loop_cond/element_types.csv"
    apec.load_element_types(csv_file_path)

    # Create a VTK unstructured grid (example grid)
    grid = vtk.vtkUnstructuredGrid()

    # Initialize the mesh with data
    apec.initialize_mesh(grid)

    # Print point data and cell data for verification
    print("Point Data Arrays:")
    for i in range(grid.GetPointData().GetNumberOfArrays()):
        print(grid.GetPointData().GetArrayName(i))

    print("Cell Data Arrays:")
    for i in range(grid.GetCellData().GetNumberOfArrays()):
        print(grid.GetCellData().GetArrayName(i))