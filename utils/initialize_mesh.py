import os
import numpy as np
import json
from ansys.mapdl import reader as pymapdl_reader
import vtk

# Get the mesh file path from the configuration
fileName = r"E:\Archive\PlasmaArc_Meshes\3D_ITER_noPlasma\3D_ITER_noPlasma.cdb"
config_file = "E:\\Archive\\PlasmaArc_Meshes\\3D_ITER_noPlasma\\conf.json"

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
        "number_substeps",
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
number_substeps = len(
    config.get("Transient Analysis Settings", {}).get("time_steps", 0)
)

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

if transient_analysis:
    num_points = archive.grid.n_points
    num_cells = archive.grid.n_cells
    # PVTU format doesn't support time varying data, only implemented for VTKHDF5
    # This is a workaround to store time varying data in a single VTK file
    # One day I will try to fix it, but I will need to rewrite also the libmesh code
    for i in range(1, number_substeps + 1):
        A_i = np.full((num_points, 3), np.nan)
        V_i = np.full(num_points, np.nan)
        archive.grid.point_data[f"A_{i}"] = A_i
        archive.grid.point_data[f"V_{i}"] = V_i

    for i in range(1, number_substeps + 1):
        Js_i = np.zeros((num_cells, 3))
        archive.grid.cell_data[f"Js_{i}"] = Js_i
else:
    # Initialize and add point datasets
    num_points = archive.grid.n_points
    A = np.full((num_points, 3), np.nan)  # Create a 2D array for the vector data
    archive.grid.point_data["A"] = A
    # Initialize and add cell datasets
    num_cells = archive.grid.n_cells
    Js = np.zeros((num_cells, 3))  # Create a 2D array for the vector data
    archive.grid.cell_data["Js"] = Js

# Print overview of data read from cdb
print(archive)
archive.grid.save(
    "E:\\Archive\\PlasmaArc_Meshes\\3D_ITER_noPlasma\\3D_ITER_noPlasma.vtu", binary=True
)
