import os
import sys

import numpy as np
import vtk


def _normalize_path(path_value):
    if not isinstance(path_value, (str, os.PathLike)):
        raise TypeError("file_path must be a string or os.PathLike")
    normalized = os.fspath(path_value)
    if not normalized:
        raise ValueError("file_path must not be empty")
    return normalized


def _remove_l2_norm_range_info_keys(grid):
    key = vtk.vtkDataArray.L2_NORM_RANGE()
    for data_getter in (grid.GetPointData, grid.GetCellData, grid.GetFieldData):
        data = data_getter()
        if data is None:
            continue
        for idx in range(data.GetNumberOfArrays()):
            arr = data.GetArray(idx)
            if arr is None:
                continue
            info = arr.GetInformation()
            if info is not None and info.Has(key):
                info.Remove(key)


def convert_mesh_from_cdb(file_path):
    file_path = _normalize_path(file_path)
    print(f"Starting CDB->VTU conversion for '{file_path}'")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Mesh file '{file_path}' not found.")

    file_extension = os.path.splitext(file_path)[1].lower()
    if file_extension != ".cdb":
        raise ValueError(f"Unsupported file extension: {file_extension}")

    try:
        from ansys.mapdl import reader as pymapdl_reader
    except ImportError as exc:
        raise ImportError(
            "Converting .cdb files requires 'ansys-mapdl-reader'. Install it with pip."
        ) from exc

    archive = pymapdl_reader.Archive(file_path)
    if archive is None:
        raise ValueError("Failed to read the mesh file.")

    grid = archive.grid
    if grid is None:
        raise ValueError("Converted archive does not contain a valid grid.")

    if not hasattr(grid, "cell_data") or not hasattr(grid, "point_data"):
        raise ValueError("Converted grid is missing required cell_data/point_data attributes.")

    for old_name, new_name in {
        "ansys_etype": "TYPE",
        "ansys_material_type": "MATERIAL",
        "ansys_real_constant": "CONSTANT",
    }.items():
        if old_name in grid.cell_data:
            grid.cell_data[new_name] = grid.cell_data[old_name]
            del grid.cell_data[old_name]

    if "COORDINATE_SYSTEM" not in grid.cell_data:
        print("Adding missing cell data array 'COORDINATE_SYSTEM' with default value 0.")
        grid.cell_data["COORDINATE_SYSTEM"] = np.zeros(grid.n_cells, dtype=np.int32)

    if "CP_0" not in grid.point_data:
        print("Adding missing point data array 'CP_0' with default value 0.")
        grid.point_data["CP_0"] = np.zeros(grid.n_points, dtype=np.int32)

    if "ansys_elem_type_num" in grid.cell_data:
        del grid.cell_data["ansys_elem_type_num"]

    if "VTKorigID" in grid.point_data:
        del grid.point_data["VTKorigID"]

    output_file = os.path.join(
        os.path.dirname(file_path),
        os.path.splitext(os.path.basename(file_path))[0] + ".vtu",
    )
    _remove_l2_norm_range_info_keys(grid)
    grid.save(output_file, binary=True)

    print(f"Processed mesh saved to '{output_file}' (points={grid.n_points} cells={grid.n_cells}).")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    elif len(sys.argv) == 1:
        input_file = input("Enter the path of the file to be converted: ").strip()
    else:
        print("Usage: python -m apec_m_mesh_converter.convert_mesh <path_to_mesh_file>")
        sys.exit(1)

    try:
        convert_mesh_from_cdb(input_file)
    except Exception as exc:
        print(f"Mesh conversion failed: {exc}")
        sys.exit(1)
