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
        self.mesh = None

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

        # Give a name to the mesh
        file_name = os.path.basename(vtk_file_path)
        file_name_without_extension = os.path.splitext(file_name)[0]
        self.mesh.SetName(file_name_without_extension)
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

    def write_multiple_cell_data_to_file(self, input_basename, array_names, num=None):
        """
        Writes multiple cell data arrays from the object's vtkUnstructuredGrid to a single binary file.
        """
        output_file_name = f"{input_basename}_loads"
        if num is not None:
            output_file_name += f"_{num}"
        output_file_name += ".bin"

        with open(output_file_name, "wb") as file:
            for array_name in array_names:
                cell_data = self.mesh.GetCellData()
                data_array = cell_data.GetArray(array_name)
                if data_array is None:
                    raise ValueError(f"Array '{array_name}' not found in cell data.")

                name_length = len(array_name)
                file.write(struct.pack("B", name_length))
                file.write(array_name.encode("utf-8"))

                num_components = data_array.GetNumberOfComponents()
                file.write(struct.pack("I", num_components))

                num_elements = data_array.GetNumberOfTuples()
                file.write(struct.pack("I", num_elements))

                for cell_id in range(num_elements):
                    file.write(struct.pack("I", cell_id))
                    values = data_array.GetTuple(cell_id)
                    file.write(struct.pack(f"{num_components}d", *values))

        print(f"All cell data arrays {array_names} written to '{output_file_name}'.")

    def write_multiple_boundary_conditions(
        self, output_basename, array_names, num=None
    ):
        """
        Writes multiple point data arrays from the object's vtkUnstructuredGrid to a single binary file.
        """
        output_file = f"{output_basename}_bc"
        if num is not None:
            output_file += f"_{num}"
        output_file += ".dat"

        with open(output_file, "wb") as file:
            for array_name in array_names:
                point_data = self.mesh.GetPointData()
                data_array = point_data.GetArray(array_name)
                if data_array is None:
                    raise ValueError(
                        f"Point data array '{array_name}' not found in the VTK grid."
                    )

                num_components = data_array.GetNumberOfComponents()

                for node_id in range(data_array.GetNumberOfTuples()):
                    values = data_array.GetTuple(node_id)

                    for component_idx, value in enumerate(values):
                        if np.isnan(value):
                            continue

                        if num_components > 1:
                            component_label = (
                                f"{array_name}_{['x', 'y', 'z'][component_idx]}"
                            )
                        else:
                            component_label = array_name

                        file.write(struct.pack("Q", node_id))
                        label = component_label.encode("utf-8")
                        label_length = len(label)
                        file.write(struct.pack("B", label_length))
                        file.write(label)
                        file.write(struct.pack("d", value))

        print(f"Boundary conditions for arrays {array_names} written to {output_file}.")

    def distribute_dofs(self):
        """
        Initialize the object's vtkUnstructuredGrid by adding variables as point data
        and element loads as cell data.
        """
        if self.element_types_collector is None:
            raise ValueError(
                "Element types have not been loaded. Call 'load_element_types' first."
            )

        element_types = self.element_types_collector.get_element_types()

        for element_id, element in element_types.items():
            variables = element.getVariables()
            for var_name, num_components in variables.items():
                array = vtk.vtkDoubleArray()
                array.SetName(var_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(self.mesh.GetNumberOfPoints())
                self.mesh.GetPointData().AddArray(array)

        for element_id, element in element_types.items():
            element_loads = element.getElementLoads()
            for load_name, num_components in element_loads.items():
                array = vtk.vtkDoubleArray()
                array.SetName(load_name)
                array.SetNumberOfComponents(num_components)
                array.SetNumberOfTuples(self.mesh.GetNumberOfCells())
                self.mesh.GetCellData().AddArray(array)

    def plot(self):
        """
        Replot the mesh in the active view.
        """
        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        print("Replotting the mesh...")
        producer = TrivialProducer()
        producer.GetClientSideObject().SetOutput(self.mesh)
        producer.GetExecutive().GetOutputInformation(0).Set(
            vtk.vtkDataObject.DATA_PIECE_NAME(), "Mesh"
        )
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


apec = APEC()
