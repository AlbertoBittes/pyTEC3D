load = "Js"  # Replace with actual load name
component = "x"  # Replace with ('x', 'y', or 'z')
value = 1.0  # Replace with chosen load value

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

component_index = 0
# Check if the array is a 3D vector
num_components = array_one.GetNumberOfComponents()
if num_components == 1:
    component_index = 0
else:
    # Map the component selector to the corresponding index
    component_map = {"x": 0, "y": 1, "z": 2}
    component_index = component_map.get(component)
    if component_index is None:
        raise ValueError("Invalid component specified.")

# Get the vtkOriginalCellIds array from the second input
original_ids = cell_data1.GetArray("vtkOriginalCellIds")
if original_ids is None:
    raise ValueError(
        "The second input VTK object does not have a 'vtkOriginalCellIds' array."
    )

# Convert VTK arrays to NumPy arrays
original_ids_np = np.array(original_ids)
if np.any((original_ids_np < 0) | (original_ids_np >= input0.GetNumberOfCells())):
    raise IndexError(
        "One or more original IDs are out of bounds for the first input VTK object."
    )

# Assign values to the loads based on the original ids
array_one_np = np.array(
    [
        array_one.GetComponent(i, component_index)
        for i in range(array_one.GetNumberOfTuples())
    ]
)

# Update the values at the specified original IDs
array_one_np[original_ids_np] = value

# Convert the NumPy array back to a VTK array
for i, idx in enumerate(original_ids_np):
    array_one.SetComponent(idx, component_index, array_one_np[idx])

# Update the cell data
cell_data0.Modified()
