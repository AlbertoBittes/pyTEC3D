import vtk
import numpy as np
import struct
import os


# Example usage
if __name__ == "__main__":
    # Specify the VTK file, point data array names, and optional integer
    vtk_file_path = "your_file.vtk"  # Replace with the path to your VTK file
    array_names = ["Array1", "Array2"]  # Replace with the names of your arrays
    num = 1  # Optional integer to append to the output file name

    # Process the boundary conditions
    process_boundary_conditions(vtk_file_path, array_names, num)