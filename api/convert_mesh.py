import os
import sys
import numpy as np
from ansys.mapdl import reader as pymapdl_reader


def process_mesh(file_path):
    """
    Processes the mesh file, modifies the grid, and saves it as a .vtu file.

    Parameters:
        file_path (str): Path to the input mesh file.
    """
    # Check if the mesh file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Mesh file '{file_path}' not found.")

    # Get the file extension
    file_extension = os.path.splitext(file_path)[1]

    # Check the file extension and read the file accordingly
    if file_extension == ".cdb":
        # Read ansys archive file
        archive = pymapdl_reader.Archive(file_path)
    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")

    # Check if the archive is successfully read
    if archive is None:
        raise ValueError("Failed to read the mesh file.")

    # Get the number of points and cells
    num_points = archive.grid.n_points
    num_cells = archive.grid.n_cells

    # Remove specified arrays if they exist in the grid
    # Check and remove from cell data
    for array_name in ["ansys_elem_type_num"]:
        if array_name in archive.grid.cell_data:
            del archive.grid.cell_data[array_name]

    # Check and remove from point data
    for array_name in ["angles", "angles_magnitude", "origid", "VTKorigID"]:
        if array_name in archive.grid.point_data:
            del archive.grid.point_data[array_name]

    # Add a point data array called boundary_id to the grid
    boundary_id = np.full(num_points, 0)
    archive.grid.point_data["boundary_id"] = boundary_id

    # Check for arrays with "PERIODIC_BOUNDARY_" or "INFINITE_BOUNDARY_" in their names
    # for array_name in archive.grid.point_data.keys():
    #     if array_name.startswith("PERIODIC_BOUNDARY_") or array_name.startswith(
    #         "INFINITE_BOUNDARY_"
    #     ):
    #         # Extract the number after the last underscore
    #         boundary_number = int(array_name.split("_")[-1])
    #         periodic_boundary = archive.grid.point_data[array_name]
    #         boundary_id = archive.grid.point_data["boundary_id"]

    #         # Assign the boundary number to points where the array value is 1
    #         boundary_id[periodic_boundary == 1] = boundary_number

    #         # Update the boundary_id in the point data
    #         archive.grid.point_data["boundary_id"] = boundary_id

    # Save the modified grid to a .vtu file
    output_file = os.path.join(
        os.path.dirname(file_path),
        os.path.splitext(os.path.basename(file_path))[0] + ".vtu",
    )
    archive.grid.save(output_file, binary=True)

    # Print overview of data read from cdb
    print(f"Processed mesh saved to '{output_file}'.")


# Allow the script to be run from the command line
if __name__ == "__main__":
    print("APEC-M Mesh Converter")
    # input_file = input("Please enter the path to the mesh file: ")

    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    elif len(sys.argv) == 1:
        # Prompt the user for the file path if no argument is provided
        input_file = input("Enter the path of the file to be converted: ").strip()
    else:
        print(
            "Usage: python -m apec_m.convert_mesh <path_to_mesh_file> or python -m apec_m.convert_mesh"
        )
        sys.exit(1)

    try:
        process_mesh(input_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
