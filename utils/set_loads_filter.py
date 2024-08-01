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
cell_data0 = input0.CellData
cell_data1 = input1.CellData

# Get the "array one" cell data array in the first input
array_one = cell_data0[load]
if array_one is None:
    raise ValueError("The first input VTK object does not have a " + load + " array.")

component_index = 0
# Check if the array is a 3D vector
num_components = array_one.shape[1] if array_one.ndim > 1 else 1
if num_components == 1:
    component_index = 0
else:
    # Map the component selector to the corresponding index
    component_map = {"x": 0, "y": 1, "z": 2}
    component_index = component_map.get(component)
    if component_index is None:
        raise ValueError("Invalid component specified.")

# Get the vtkOriginalCellIds array from the second input
original_ids = cell_data1["vtkOriginalCellIds"]
if original_ids is None:
    raise ValueError(
        "The second input VTK object does not have a 'vtkOriginalCellIds' array."
    )

# Ensure original_ids are within bounds
if np.any((original_ids < 0) | (original_ids >= input0.GetNumberOfCells())):
    raise IndexError(
        "One or more original IDs are out of bounds for the first input VTK object."
    )

# Update the values at the specified original IDs
array_one[original_ids, component_index] = value

# Update the cell data
cell_data0.Modified()