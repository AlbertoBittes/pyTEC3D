import vtk

class TimeStepGenerator(vtk.vtkPythonAlgorithmBase):
    def __init__(self, time_steps, input_grid):
        super(TimeStepGenerator, self).__init__()
        self.SetNumberOfInputPorts(0)
        self.SetNumberOfOutputPorts(1)
        self.time_steps = time_steps
        self.input_grid = input_grid

    def RequestInformation(self, request, inInfo, outInfo):
        outInfo = outInfo.GetInformationObject(0)
        
        # Specify available time steps
        outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(), self.time_steps, len(self.time_steps))
        
        # Specify time range
        timeRange = [self.time_steps[0], self.time_steps[-1]]
        outInfo.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), timeRange, 2)
        
        return 1

    def RequestData(self, request, inInfo, outInfo):
        outInfo = outInfo.GetInformationObject(0)
        
        # Get the requested time step
        timeStep = 0.0
        if outInfo.Has(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()):
            timeStep = outInfo.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())
        
        # Create a new unstructured grid for the output
        outputGrid = vtk.vtkUnstructuredGrid()
        outputGrid.DeepCopy(self.input_grid)
        
        # Create and fill PointData
        pointDataArray = vtk.vtkDoubleArray()
        pointDataArray.SetName("PointData")
        pointDataArray.SetNumberOfComponents(1)
        pointDataArray.SetNumberOfTuples(self.input_grid.GetNumberOfPoints())
        
        for i in range(self.input_grid.GetNumberOfPoints()):
            value = float(i) + timeStep  # Example: value changes with time
            pointDataArray.SetValue(i, value)
        outputGrid.GetPointData().AddArray(pointDataArray)
        
        # Create and fill CellData
        cellDataArray = vtk.vtkDoubleArray()
        cellDataArray.SetName("CellData")
        cellDataArray.SetNumberOfComponents(1)
        cellDataArray.SetNumberOfTuples(self.input_grid.GetNumberOfCells())
        
        for i in range(self.input_grid.GetNumberOfCells()):
            value = float(i) + timeStep  # Example: value changes with time
            cellDataArray.SetValue(i, value)
        outputGrid.GetCellData().AddArray(cellDataArray)
        
        # Set the output data
        outInfo.Set(vtk.vtkDataObject.DATA_OBJECT(), outputGrid)
        
        return 1

class TimeStepModifier(vtk.vtkPythonAlgorithmBase):
    def __init__(self):
        super(TimeStepModifier, self).__init__()
        self.SetNumberOfInputPorts(1)
        self.SetNumberOfOutputPorts(1)

    def RequestData(self, request, inInfo, outInfo):
        inInfo = inInfo[0].GetInformationObject(0)
        outInfo = outInfo.GetInformationObject(0)
        
        # Get the input unstructured grid
        inputGrid = vtk.vtkUnstructuredGrid.SafeDownCast(inInfo.Get(vtk.vtkDataObject.DATA_OBJECT()))
        
        # Create a new unstructured grid for the output
        outputGrid = vtk.vtkUnstructuredGrid()
        outputGrid.DeepCopy(inputGrid)
        
        # Modify PointData
        pointDataArray = inputGrid.GetPointData().GetArray("PointData")
        for i in range(pointDataArray.GetNumberOfTuples()):
            value = pointDataArray.GetValue(i)
            pointDataArray.SetValue(i, value * 2)  # Example modification: double the value
        
        # Modify CellData
        cellDataArray = inputGrid.GetCellData().GetArray("CellData")
        for i in range(cellDataArray.GetNumberOfTuples()):
            value = cellDataArray.GetValue(i)
            cellDataArray.SetValue(i, value * 2)  # Example modification: double the value
        
        # Set the output data
        outInfo.Set(vtk.vtkDataObject.DATA_OBJECT(), outputGrid)
        
        return 1

# User input for time steps
time_steps = [0.0, 0.5, 1.0, 1.5, 2.0]

# User-provided unstructured grid
input_grid = vtk.vtkUnstructuredGrid()

# Create points
points = vtk.vtkPoints()
points.InsertNextPoint(0.0, 0.0, 0.0)
points.InsertNextPoint(1.0, 0.0, 0.0)
points.InsertNextPoint(1.0, 1.0, 0.0)
points.InsertNextPoint(0.0, 1.0, 0.0)

# Create a cell
cellArray = vtk.vtkCellArray()
cell = vtk.vtkQuad()
for i in range(4):
    cell.GetPointIds().SetId(i, i)
cellArray.InsertNextCell(cell)

input_grid.SetPoints(points)
input_grid.SetCells(cell.GetCellType(), cellArray)

# Create instances of the custom algorithms
timeStepGenerator = TimeStepGenerator(time_steps, input_grid)
timeStepModifier = TimeStepModifier()

# Set up the pipeline
timeStepModifier.SetInputConnection(timeStepGenerator.GetOutputPort())

# Loop over all time steps and execute the pipeline
for time_step in time_steps:
    executive = vtk.vtkStreamingDemandDrivenPipeline.SafeDownCast(timeStepModifier.GetExecutive())
    executive.UpdateTimeStep(time_step)
    timeStepModifier.Update()
    
    # Get the output
    outputGrid = vtk.vtkUnstructuredGrid.SafeDownCast(timeStepModifier.GetOutput())
    
    # Write the output to a file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(f"output_time_step_{time_step}.vtu")
    writer.SetInputData(outputGrid)
    writer.Write()
    
    # Print the modified point and cell data for each time step
    pointDataArray = outputGrid.GetPointData().GetArray("PointData")
    cellDataArray = outputGrid.GetCellData().GetArray("CellData")
    
    print(f"Time Step: {time_step}")
    print("Point Data:", [pointDataArray.GetValue(i) for i in range(pointDataArray.GetNumberOfTuples())])
    print("Cell Data:", [cellDataArray.GetValue(i) for i in range(cellDataArray.GetNumberOfTuples())])
    print()