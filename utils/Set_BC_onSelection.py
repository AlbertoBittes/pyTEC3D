from paraview.simple import *
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# Load the custom filter class
exec(open("D:\\PlasmaArc\\SetLoadBCFilter.py").read())

# Get the selected sources
selected_sources = GetSelectedSources()

if len(selected_sources) != 2:
    raise ValueError("Please select exactly two sources.")

input0 = selected_sources[0]
input1 = selected_sources[1]

# Ask the user for the dof name and value
dof = input("Enter the dof name: ")
value = float(input("Enter the value: "))

# Create an instance of the custom filter
filter = SetLoadBCFilter()
filter.dof = dof
filter.value = value

# Set the inputs
filter.SetInputDataObject(0, input0.GetClientSideObject().GetOutputDataObject(0))
filter.SetInputDataObject(1, input1.GetClientSideObject().GetOutputDataObject(0))

# Update the filter
filter.Update()

# Get the output
output = filter.GetOutput()

# Create a new ParaView source from the output
output_source = PVTrivialProducer()
output_source.GetClientSideObject().SetOutput(output)

# Add the output to the pipeline browser
Show(output_source)

# Remove the original inputs from the pipeline browser
Delete(input0)
Delete(input1)

# Render the view
Render()