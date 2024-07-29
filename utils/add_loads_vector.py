import vtk

def add_time_dependent_vector_point_data(vtu_file, output_file, time_steps):
    """
    Adds time-dependent vector point data to an existing VTK unstructured grid, initializing all arrays to zero.

    Parameters:
    - vtu_file: str, path to the input VTU file.
    - output_file: str, path to the output VTU file.
    - time_steps: list of float, the time steps for which to add vector data.
    """
    # Read the input VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()
    grid = reader.GetOutput()

    for time_step in time_steps:
        # Create a vector point data array
        vector_data = vtk.vtkDoubleArray()
        vector_data.SetNumberOfComponents(3)
        vector_data.SetName(f"VectorData_t{time_step}")

        # Initialize vector data to zero for each point
        zero_vector = [0.0, 0.0, 0.0]
        for i in range(grid.GetNumberOfPoints()):
            vector_data.InsertNextTuple(zero_vector)

        # Add the vector point data array to the grid
        grid.GetPointData().AddArray(vector_data)

        # Mark the grid as modified
        grid.Modified()

    # Write the modified grid to a file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(grid)
    writer.Write()

# Example usage
vtu_file = "input.vtu"  # Replace with your input VTU file path
output_file = "output.vtu"  # Replace with your desired output VTU file path
time_steps = [0.0, 1.0, 2.0, 3.0]  # Example time steps
add_time_dependent_vector_point_data(vtu_file, output_file, time_steps)
