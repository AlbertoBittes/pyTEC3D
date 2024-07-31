import vtk

# Create a sample VTK unstructured grid object if you don't have one
def create_sample_unstructured_grid():
    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(0.0, 1.0, 0.0)
    points.InsertNextPoint(0.0, 0.0, 1.0)

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    tetra = vtk.vtkTetra()
    tetra.GetPointIds().SetId(0, 0)
    tetra.GetPointIds().SetId(1, 1)
    tetra.GetPointIds().SetId(2, 2)
    tetra.GetPointIds().SetId(3, 3)

    grid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())
    return grid

# Function to add time-dependent vector cell data initialized to zero
def add_time_dependent_vector_cell_data(grid, num_timesteps):
    num_cells = grid.GetNumberOfCells()
    
    # Create a vtkTemporalDataSet to store each timestep
    temporal_data_set = vtk.vtkTemporalDataSet()
    
    for t in range(num_timesteps):
        # Create a copy of the grid for each timestep
        grid_copy = vtk.vtkUnstructuredGrid()
        grid_copy.DeepCopy(grid)
        
        # Create a vtkFloatArray to store the vector data
        vector_data = vtk.vtkFloatArray()
        vector_data.SetNumberOfComponents(3)  # Vector has 3 components (x, y, z)
        vector_data.SetName("TimeDependentVectorData")
        
        # Initialize the vector data to zero for each cell
        for i in range(num_cells):
            vector_data.InsertNextTuple3(0.0, 0.0, 0.0)
        
        # Add the vector data to the cell data of the unstructured grid
        grid_copy.GetCellData().AddArray(vector_data)
        
        # Create a vtkTemporalDataSetElement to store the grid and the time step
        temporal_element = vtk.vtkTemporalDataSetElement()
        temporal_element.SetTimeStep(t)
        temporal_element.SetDataObject(grid_copy)
        
        # Add the temporal element to the temporal data set
        temporal_data_set.SetTimeStep(t, temporal_element)
    
    return temporal_data_set

# Example usage
unstructured_grid = create_sample_unstructured_grid()
num_timesteps = 10  # Define the number of timesteps
temporal_data_set = add_time_dependent_vector_cell_data(unstructured_grid, num_timesteps)

# Print the cell data to verify
for t in range(num_timesteps):
    temporal_element = temporal_data_set.GetTimeStep(t)
    grid = vtk.vtkUnstructuredGrid.SafeDownCast(temporal_element.GetDataObject())
    cell_data = grid.GetCellData().GetArray("TimeDependentVectorData")
    print(f"Time step {t}:")
    for i in range(cell_data.GetNumberOfTuples()):
        print(cell_data.GetTuple3(i))