
dof = "ax_bc"  # Replace with actual dof name
value = 1.0  # Replace with actual value

import vtk
import numpy as np
# Example usage
input0 = input[0]  # Replace with actual input0
input1 = input[1]  # Replace with actual input1

# Get the input data objects
if input0 is None or input1 is None:
    raise ValueError("One or both input VTK objects are missing.")

# Get the point data arrays
point_data0 = input0.GetPointData()
point_data1 = input1.GetPointData()

# Get the "array one" point data array in the first input
array_one = point_data0.GetArray(dof)
if array_one is None:
    raise ValueError("The first input VTK object does not have a " + dof + " array.")


# Get the vtkOriginalIds array from the second input
original_ids = point_data1.GetArray('vtkOriginalIds').astype(int)
if original_ids is None:
    raise ValueError("The second input VTK object does not have a 'vtkOriginalIds' array.")

# Convert VTK arrays to NumPy arrays

# Ensure original_ids_np values are within bounds
if np.any((original_ids < 0) | (original_ids >= input0.GetNumberOfPoints())):
    raise IndexError("One or more original IDs are out of bounds for the first input VTK object.")

# Assign values to the dofs based on the original ids
array_one[original_ids] = value

point_data0.append(array_one, dof)
point_data0.Modified()

# Set the modified input0 as the output
output.PointData.append(array_one, dof)

