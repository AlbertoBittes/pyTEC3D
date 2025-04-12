import os
import numpy as np
import pyvista as pv
import json
from ansys.mapdl import reader as pymapdl_reader
import vtk

# Get the mesh file path from the configuration
fileName = "3D_loop_cond_period.cdb"

# Check if the mesh file exists
if not os.path.exists(fileName):
    raise FileNotFoundError(f"Mesh file '{fileName}' not found.")

# Get the file extension
file_extension = os.path.splitext(fileName)[1]

# Initialize the archive variable
archive = None

# Check the file extension and read the file accordingly
if file_extension == ".cdb":
    # Read ansys archive file
    archive = pymapdl_reader.Archive(fileName)
elif file_extension == ".vtu":
    # Read VTK grid file using vtkUnstructuredGridReader
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    archive = reader.GetOutput()
else:
    raise ValueError(f"Unsupported file extension: {file_extension}")

# Check if the archive is successfully read
if archive is None:
    raise ValueError("Failed to read the mesh file.")


# Example usage:
# Assuming `unstructured_grid` is your vtkUnstructuredGrid object
dimension = get_spatial_dimension(archive)

num_points = archive.grid.n_points
num_cells = archive.grid.n_cells
# PVTU format doesn't support time varying data, only implemented for VTKHDF5
# This is a workaround to store time varying data in a single VTK file
# One day I will try to fix it, but I will need to rewrite also the libmesh code
for i in range(1, number_time_steps + 1):
    if dimension == 3:
        A_i = np.full((num_points, 3), np.nan)
        V_i = np.full(num_points, np.nan)
    elif dimension == 2:
        A_i = np.full((num_points, 1), np.nan)
        # V_i = np.full(num_points, np.nan)
    archive.grid.point_data[f"A_{i}"] = A_i
    if dimension == 3:
        archive.grid.point_data[f"V_{i}"] = V_i

for i in range(1, number_time_steps + 1):
    if dimension == 3:
        Js_i = np.zeros((num_cells, 3))
    elif dimension == 2:
        Js_i = np.zeros((num_cells, 1))
    archive.grid.cell_data[f"Js_{i}"] = Js_i

# Add a point data array called boundary_id to the grid
boundary_id = np.full(num_points, 0)
archive.grid.point_data["boundary_id"] = boundary_id
# Check if PERIODIC_BOUNDARY_1 is present in the point data
if "PERIODIC_BOUNDARY_1" in archive.grid.point_data:
    periodic_boundary_1 = archive.grid.point_data["PERIODIC_BOUNDARY_1"]
    boundary_id = archive.grid.point_data["boundary_id"]
    
    # Assign the value 7 to all points where PERIODIC_BOUNDARY_1 is 1
    boundary_id[periodic_boundary_1 == 1] = 7
    
    # Update the boundary_id in the point data
    archive.grid.point_data["boundary_id"] = boundary_id

# Check if PERIODIC_BOUNDARY_2 is present in the point data
if "PERIODIC_BOUNDARY_2" in archive.grid.point_data:
    periodic_boundary_2 = archive.grid.point_data["PERIODIC_BOUNDARY_2"]
    boundary_id = archive.grid.point_data["boundary_id"]
    
    # Assign the value 8 to all points where PERIODIC_BOUNDARY_2 is 1
    boundary_id[periodic_boundary_2 == 1] = 8
    
    # Update the boundary_id in the point data
    archive.grid.point_data["boundary_id"] = boundary_id
# else:
#     # Initialize and add point datasets
#     num_points = archive.grid.n_points
#     if dimension == 3:
#         A = np.full((num_points, 3), np.nan)
#     elif dimension == 2:
#         A = np.full((num_points, 1), np.nan)
#     archive.grid.point_data["A"] = A
#     # Initialize and add cell datasets
#     num_cells = archive.grid.n_cells
#     if dimension == 3:
#         Js = np.zeros((num_cells, 3))
#     elif dimension == 2:
#         Js = np.zeros((num_cells, 1))
#     archive.grid.cell_data["Js"] = Js

# Print overview of data read from cdb
print(archive)
archive.grid.save("3D_loop_cond_period.vtu", binary=False)
