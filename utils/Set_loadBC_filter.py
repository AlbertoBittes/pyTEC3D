import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

class SetLoadBCFilter(vtk.vtkAlgorithm):
    def __init__(self):
        self.SetNumberOfInputPorts(2)
        self.SetNumberOfOutputPorts(1)

    def RequestData(self, request, inInfo, outInfo):
        dof = self.dof
        value = self.value

        # Get the input data objects
        input0 = vtk.vtkDataSet.GetData(inInfo[0])
        input1 = vtk.vtkDataSet.GetData(inInfo[1])

        if input0 is None or input1 is None:
            raise ValueError("One or both input VTK objects are missing.")

        # Get the point data arrays
        point_data0 = input0.GetPointData()
        point_data1 = input1.GetPointData()

        # Get the "array one" point data array in the first input
        array_one = point_data0.GetArray(dof)
        if array_one is None:
            array_one = vtk.vtkDoubleArray()
            array_one.SetName(dof)
            array_one.SetNumberOfComponents(1)
            array_one.SetNumberOfTuples(input0.GetNumberOfPoints())
            point_data0.AddArray(array_one)

        # Get the vtkOriginalIds array from the second input
        original_ids = point_data1.GetArray('vtkOriginalIds')
        if original_ids is None:
            raise ValueError("The second input VTK object does not have a 'vtkOriginalIds' array.")

        # Convert VTK arrays to NumPy arrays
        original_ids_np = vtk_to_numpy(original_ids)
        array_one_np = vtk_to_numpy(array_one)

        # Ensure original_ids_np values are within bounds
        if np.any((original_ids_np < 0) | (original_ids_np >= input0.GetNumberOfPoints())):
            raise IndexError("One or more original IDs are out of bounds for the first input VTK object.")

        # Assign values to the dofs based on the original ids
        array_one_np[original_ids_np] = value

        # Update the VTK array with the modified NumPy array
        array_one.Modified()

        # Set the modified input0 as the output
        output = vtk.vtkDataSet.GetData(outInfo)
        output.ShallowCopy(input0)

        return 1