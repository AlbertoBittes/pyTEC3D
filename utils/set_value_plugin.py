import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input", port_index=0)
@smproperty.input(name="Selection", port_index=1)
class SetValueOnSelection(vtk.vtkPythonAlgorithmBase):
    def __init__(self):
        self.__data_type = "Point Data"
        self.__array_name = ""
        self.__is_vector = False
        self.__component = 0
        self.__value = 1.0
        vtk.vtkPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1)

    def RequestData(self, request, inInfo, outInfo):
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        selection_data = vtk.vtkUnstructuredGrid.GetData(inInfo[1])
        output_data = vtk.vtkUnstructuredGrid.GetData(outInfo)
        output_data.ShallowCopy(input_data)

        if self.__data_type == "Point Data":
            original_ids = vtk_to_numpy(selection_data.GetPointData().GetArray("vtkOriginalPointIds"))
            data_array = vtk_to_numpy(input_data.GetPointData().GetArray(self.__array_name))
        else:
            original_ids = vtk_to_numpy(selection_data.GetCellData().GetArray("vtkOriginalCellIds"))
            data_array = vtk_to_numpy(input_data.GetCellData().GetArray(self.__array_name))

        if original_ids is None or data_array is None:
            raise ValueError("The selection or data array is not available.")

        if self.__is_vector:
            data_array[original_ids, self.__component] = self.__value
        else:
            data_array[original_ids] = self.__value

        if self.__data_type == "Point Data":
            vtk_array = numpy_to_vtk(data_array)
            vtk_array.SetName(self.__array_name)
            input_data.GetPointData().AddArray(vtk_array)
        else:
            vtk_array = numpy_to_vtk(data_array)
            vtk_array.SetName(self.__array_name)
            input_data.GetCellData().AddArray(vtk_array)

        return 1

    @smproperty.stringvector(name="Data Type", default_values="Point Data")
    @smdomain.xml("<StringListDomain><String value='Point Data'/><String value='Cell Data'/></StringListDomain>")
    @smproperty.documentation("Select the type of data to modify: 'Point Data' or 'Cell Data'.")
    def SetDataType(self, data_type):
        if data_type not in ["Point Data", "Cell Data"]:
            raise ValueError("Invalid data type. Available values are 'Point Data' and 'Cell Data'.")
        if self.__data_type != data_type:
            self.__data_type = data_type
            self.Modified()

    @smproperty.intvector(name="Is Vector", default_values=0)
    @smdomain.xml("<BooleanDomain/>")
    @smproperty.documentation("Specify if the quantity is a vector (True) or a scalar (False).")
    def SetIsVector(self, is_vector):
        self.__is_vector = bool(is_vector)
        self.Modified()

    @smproperty.stringvector(name="Array Name", default_values="")
    @smdomain.xml("<ArrayListDomain attribute_type='Scalars'/>")
    @smproperty.documentation("Specify the name of the array to modify.")
    def SetArrayName(self, array_name):
        input_data = self.GetInputDataObject(0, 0)
        if self.__data_type == "Point Data":
            data_array = input_data.GetPointData().GetArray(array_name)
        else:
            data_array = input_data.GetCellData().GetArray(array_name)

        if data_array is None:
            raise ValueError(f"The array '{array_name}' does not exist in the {self.__data_type.lower()} data.")

        self.__array_name = array_name
        self.Modified()

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.xml("""
        <RangeDomain min='0' max='2'/>
        <EnableIf name='IsVector' value='1'/>
    """)
    @smproperty.documentation("Specify the component index of the array to modify. Only applicable if the quantity is a vector.")
    def SetComponent(self, component):
        if not self.__is_vector:
            raise ValueError("Component can only be set if the quantity is a vector.")
        
        input_data = self.GetInputDataObject(0, 0)
        if self.__data_type == "Point Data":
            data_array = input_data.GetPointData().GetArray(self.__array_name)
        else:
            data_array = input_data.GetCellData().GetArray(self.__array_name)

        if data_array is None:
            raise ValueError(f"The array '{self.__array_name}' does not exist in the {self.__data_type.lower()} data.")

        if component < 0 or component >= data_array.GetNumberOfComponents():
            raise ValueError(f"The component index {component} is out of range for the array '{self.__array_name}'.")

        self.__component = component
        self.Modified()