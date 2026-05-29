import os
import sys
import numpy as np
from ansys.mapdl import reader as pymapdl_reader

from .logging_utils import get_logger
from .validation_utils import ensure_pathlike


LOGGER = get_logger("apec_m.convert_mesh")


def convert_mesh_from_cdb(file_path):
    """
    Processes the mesh file, modifies the grid, and saves it as a .vtu file.

    Args:
        file_path: Path to the input mesh file.

    Returns:
        None
    """
    file_path = ensure_pathlike(file_path, "file_path")
    LOGGER.info("Starting CDB->VTU conversion for '%s'.", file_path)

    # Check if the mesh file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Mesh file '{file_path}' not found.")

    # Get the file extension
    file_extension = os.path.splitext(file_path)[1].lower()

    # Check the file extension and read the file accordingly
    if file_extension == ".cdb":
        # Read ansys archive file
        archive = pymapdl_reader.Archive(file_path)
    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")

    # Check if the archive is successfully read
    if archive is None:
        raise ValueError("Failed to read the mesh file.")

    if getattr(archive, "grid", None) is None:
        raise ValueError("Converted archive does not contain a valid grid.")

    if not hasattr(archive.grid, "cell_data") or not hasattr(archive.grid, "point_data"):
        raise ValueError("Converted grid is missing required cell_data/point_data attributes.")

    # Normalize selected ANSYS cell-data names to the names used elsewhere by APEC-M.
    for old_name, new_name in {
        "ansys_etype": "TYPE",
        "ansys_material_type": "MATERIAL",
        "ansys_real_constant": "CONSTANT",
    }.items():
        if old_name in archive.grid.cell_data:
            archive.grid.cell_data[new_name] = archive.grid.cell_data[old_name]
            del archive.grid.cell_data[old_name]

    # Remove specified arrays if they exist in the grid
    # Check and remove from cell data
    for array_name in ["ansys_elem_type_num"]:
        if array_name in archive.grid.cell_data:
            del archive.grid.cell_data[array_name]

    # Check and remove from point data
    for array_name in ["VTKorigID"]:
        if array_name in archive.grid.point_data:
            del archive.grid.point_data[array_name]

    # Save the modified grid to a .vtu file
    output_file = os.path.join(
        os.path.dirname(file_path),
        os.path.splitext(os.path.basename(file_path))[0] + ".vtu",
    )
    archive.grid.save(output_file, binary=True)

    # Print overview of data read from cdb
    LOGGER.info(
        "Processed mesh saved to '%s' (points=%s cells=%s).",
        output_file,
        archive.grid.n_points,
        archive.grid.n_cells,
    )


# Allow the script to be run from the command line
if __name__ == "__main__":
    LOGGER.info("APEC-M Mesh Converter")
    # input_file = input("Please enter the path to the mesh file: ")

    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    elif len(sys.argv) == 1:
        # Prompt the user for the file path if no argument is provided
        input_file = input("Enter the path of the file to be converted: ").strip()
    else:
        LOGGER.error(
            "Usage: python -m apec_m.convert_mesh <path_to_mesh_file> or python -m apec_m.convert_mesh"
        )
        sys.exit(1)

    try:
        convert_mesh_from_cdb(input_file)
    except Exception as e:
        LOGGER.exception("Mesh conversion failed: %s", e)
        sys.exit(1)
