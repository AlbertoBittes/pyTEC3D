import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

# same imports as earlier.
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase


# ------------------------------------------------------------------------------
# A reader example.
# ------------------------------------------------------------------------------
def createModifiedCallback(anobject):
    import weakref

    weakref_obj = weakref.ref(anobject)
    anobject = None

    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()

    return _markmodified


@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input Grid", port_index=0)
@smproperty.input(name="Extracted Selection", port_index=1)
class SetValueOnSelection(VTKPythonAlgorithmBase):
    def __init__(self):
        self.__component = 0
        self.__value = 0.0
        self.__array_names = []
        self.__arrayselection = vtkDataArraySelection()
        self.__arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1)

    def RequestInformation(self, request, inInfo, outInfo):
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        if input_data is None:
            return 0

        point_data_arrays = [
            input_data.GetPointData().GetArrayName(i)
            for i in range(input_data.GetPointData().GetNumberOfArrays())
        ]
        cell_data_arrays = [
            input_data.GetCellData().GetArrayName(i)
            for i in range(input_data.GetCellData().GetNumberOfArrays())
        ]
        self.__array_names = point_data_arrays + cell_data_arrays

        for array_name in self.__array_names:
            self.__arrayselection.AddArray(array_name)
            # self.__arrayselection.EnableArray(array_name)  # Enable the array if needed

        self.__arrayselection.Modified()

        return 1

    def RequestData(self, request, inInfo, outInfo):
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        selection_data = vtk.vtkUnstructuredGrid.GetData(inInfo[1])
        self.__arrayselection.Modified()
        for i in range(self.__arrayselection.GetNumberOfArrays()):
            array_name = self.__arrayselection.GetArrayName(i)
            if not self.__arrayselection.ArrayIsEnabled(array_name):
                continue

            # Try to get the array from point data first
            vtk_data_array = input_data.GetPointData().GetArray(array_name)
            if vtk_data_array is not None:
                data_array = vtk_to_numpy(vtk_data_array)
                original_ids = vtk_to_numpy(
                    selection_data.GetPointData().GetArray("vtkOriginalPointIds")
                )
                is_point_data = True
            else:
                # If not found in point data, try to get it from cell data
                vtk_data_array = input_data.GetCellData().GetArray(array_name)
                if vtk_data_array is not None:
                    data_array = vtk_to_numpy(vtk_data_array)
                    original_ids = vtk_to_numpy(
                        selection_data.GetCellData().GetArray("vtkOriginalCellIds")
                    )
                    is_point_data = False
                else:
                    raise ValueError(
                        f"The data array '{array_name}' is not available in either point data or cell data."
                    )

            if original_ids is None or data_array is None:
                raise ValueError(
                    f"The selection or data array '{array_name}' is not available."
                )

            # Check if the array is a vector
            if vtk_data_array.GetNumberOfComponents() > 1:
                data_array[original_ids, self.__component] = self.__value
            else:
                data_array[original_ids] = self.__value

            vtk_array = numpy_to_vtk(data_array)
            vtk_array.SetName(array_name)
            vtk_array.Modified()
            if is_point_data:
                input_data.GetPointData().AddArray(vtk_array)
            else:
                input_data.GetCellData().AddArray(vtk_array)

        input_data.Modified()

        return 1

    def _get_array_selection(self):
        return self.__arrayselection

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()

    @smproperty.intvector(name="Component: ", default_values=0)
    @smdomain.intrange(min=0, max=2)
    def SetComponent(self, component):
        self.__component = component
        self.Modified()

    @smproperty.doublevector(name="Value: ", default_values=0.0)
    @smdomain.doublerange()
    def SetValue(self, value):
        self.__value = value
        self.Modified()
