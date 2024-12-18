import os
import numpy as np
import pyvista as pv
import json
from ansys.mapdl import reader as pymapdl_reader
import vtk

# Get the mesh file path from the configuration
fileName = "3D_loop_cond.cdb"
config_file = "config_3D_loop.json"

# Check if the configuration file exists
if not os.path.exists(config_file):
    raise FileNotFoundError(f"Configuration file '{config_file}' not found.")

with open(config_file, "r") as f:
    config = json.load(f)

# Check and warn if any of the settings in the JSON are not defined
required_settings = {
    "Transient Analysis Settings": [
        "transient_analysis",
        "time_steps",
    ],
}

for section, keys in required_settings.items():
    if section not in config:
        raise ValueError(f"'{section}' section is missing in the configuration file.")
    else:
        for key in keys:
            if key not in config[section]:
                raise ValueError(
                    f"'{key}' is missing in the '{section}' section of the configuration file."
                )

# Get the transient_analysis setting
transient_analysis = config.get("Transient Analysis Settings", {}).get(
    "transient_analysis", False
)

# Get the time_steps setting
time_steps = config.get("Transient Analysis Settings", {}).get("time_steps", [])

# Get the number of substeps
if transient_analysis:
    number_time_steps = len(
        config.get("Transient Analysis Settings", {}).get("time_steps", 0)
    )
else:
    number_time_steps = 0

# Check if the mesh file path is provided
if not fileName:
    raise ValueError("Mesh file path is not provided in the configuration file.")

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
elif file_extension == ".vtk":
    # Read VTK grid file using vtkUnstructuredGridReader
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    archive = reader.GetOutput()
else:
    raise ValueError(f"Unsupported file extension: {file_extension}")

# Check if the archive is successfully read
if archive is None:
    raise ValueError("Failed to read the mesh file.")


def get_spatial_dimension(unstructured_grid):
    if isinstance(unstructured_grid.grid, pv.UnstructuredGrid):
        bounds = unstructured_grid.grid.bounds
        if bounds[4] == bounds[5]:  # Check if Z bounds are the same
            return 2
        else:
            return 3
    elif isinstance(unstructured_grid, vtk.vtkUnstructuredGrid):
        unstructured_grid.GetBounds()
        if bounds[4] == bounds[5]:  # Check if Z bounds are the same
            return "2D"
        else:
            return "3D"
    else:
        raise TypeError("The provided object is not a pyvista.UnstructuredGrid.")


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
        V_i = np.full(num_points, np.nan)
    archive.grid.point_data[f"A_{i}"] = A_i
    archive.grid.point_data[f"V_{i}"] = V_i

for i in range(1, number_time_steps + 1):
    if dimension == 3:
        Js_i = np.zeros((num_cells, 3))
    elif dimension == 2:
        Js_i = np.zeros((num_cells, 1))
    archive.grid.cell_data[f"Js_{i}"] = Js_i
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
archive.grid.save("3D_loop_cond.vtu", binary=False)
