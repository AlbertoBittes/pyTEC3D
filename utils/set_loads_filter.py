load = "jx"  # Replace with actual load name
value = 1.0  # Replace with actual value

import vtk
import numpy as np

# Example usage
input0 = inputs[0]  # Replace with actual input0
input1 = inputs[1]  # Replace with actual input1

# Get the input data objects
if input0 is None or input1 is None:
    raise ValueError("One or both input VTK objects are missing.")

# Get the cell data arrays
cell_data0 = input0.GetCellData()
cell_data1 = input1.GetCellData()

# Get the "array one" cell data array in the first input
array_one = cell_data0.GetArray(load)
if array_one is None:
    raise ValueError("The first input VTK object does not have a " + load + " array.")

# Get the vtkOriginalCellIds array from the second input
original_ids = cell_data1.GetArray("vtkOriginalCellIds")
if original_ids is None:
    raise ValueError("The second input VTK object does not have a 'vtkOriginalCellIds' array.")

# Ensure original_ids values are within bounds
if np.any((original_ids < 0) | (original_ids >= input0.GetNumberOfCells())):
    raise IndexError("One or more original IDs are out of bounds for the first input VTK object.")

array_one[original_ids] = value

# Update the VTK array with the modified NumPy array
array_one.Modified()

# Append the modified array to cell_data0
cell_data0.append(array_one, load)
cell_data0.Modified()
