import vtk
import struct
import os


# Example usage
if __name__ == "__main__":
    vtk_file_path = "example.vtk"  # Path to the VTK file
    array_names = ["LoadArray1", "LoadArray2"]  # List of cell data array names to write
    num = 1  # Optional integer to append to the output file name

    process_vtk_grid(vtk_file_path, array_names, num)