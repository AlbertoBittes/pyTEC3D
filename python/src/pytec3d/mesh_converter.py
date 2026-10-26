"""Conversion of ANSYS ``.cdb`` meshes to pytec3d ``.vtu`` meshes.

Needs ``ansys-mapdl-reader``, installed by ``build.sh`` (by hand: ``pip install --no-deps ansys-mapdl-reader appdirs
tqdm``: its declared requirements conflict with the conda packages of the solver). Command line::

    pytec3d-convert-cdb mesh.cdb      # writes mesh.vtu next to it
"""

from __future__ import annotations

import argparse
import os
import sys

from vtkmodules.vtkCommonCore import vtkDataArray

from . import mesh as meshes
from ._logging import get_logger
from ._typing import PathType

LOGGER = get_logger("pytec3d.mesh_converter")

#: ANSYS cell arrays and their pytec3d names.
_RENAMED_CELL_ARRAYS = {"ansys_etype": "TYPE", "ansys_material_type": "MATERIAL", "ansys_real_constant": "CONSTANT"}
#: ANSYS arrays pytec3d does not use.
_DROPPED_CELL_ARRAYS = ("ansys_elem_type_num",)
_DROPPED_POINT_ARRAYS = ("VTKorigID",)


def _remove_norm_range_keys(grid: object) -> None:
    """Remove the cached L2 norm ranges, which VTK readers otherwise complain about."""
    key = vtkDataArray.L2_NORM_RANGE()
    for data in (grid.GetPointData(), grid.GetCellData(), grid.GetFieldData()):  # type: ignore[attr-defined]
        for i in range(data.GetNumberOfArrays()):
            array = data.GetArray(i)
            info = array.GetInformation() if array is not None else None
            if info is not None and info.Has(key):
                info.Remove(key)


def convert_cdb(path: PathType) -> str:
    """Convert an ANSYS ``.cdb`` mesh to a ``.vtu`` file next to it.

    The ANSYS element type, material and real constant ids become the TYPE, MATERIAL and CONSTANT cell arrays, and
    the grid gets the arrays of a model (:func:`pytec3d.mesh.update_grid`, for every element type of the library).

    Args:
        path: The ``.cdb`` file.

    Returns:
        Path of the ``.vtu`` file written.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If it is not a ``.cdb`` file or holds no grid.
        ImportError: If ``ansys-mapdl-reader`` is not installed.
    """
    file_path = os.fspath(path)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Mesh file '{file_path}' not found.")
    if os.path.splitext(file_path)[1].lower() != ".cdb":
        raise ValueError(f"'{file_path}' is not a .cdb file.")
    try:
        from ansys.mapdl import reader as mapdl_reader
    except ImportError as exc:
        raise ImportError(
            "Converting .cdb files needs 'ansys-mapdl-reader': pip install --no-deps ansys-mapdl-reader appdirs tqdm."
        ) from exc

    LOGGER.info("Converting '%s'.", file_path)
    grid = mapdl_reader.Archive(file_path).grid
    if grid is None:
        raise ValueError(f"'{file_path}' holds no grid.")
    for old, new in _RENAMED_CELL_ARRAYS.items():
        if old in grid.cell_data:
            grid.cell_data[new] = grid.cell_data.pop(old)
    for name in _DROPPED_CELL_ARRAYS:
        grid.cell_data.pop(name, None)
    for name in _DROPPED_POINT_ARRAYS:
        grid.point_data.pop(name, None)
    meshes.update_grid(grid)

    output = os.path.splitext(file_path)[0] + ".vtu"
    _remove_norm_range_keys(grid)
    grid.save(output, binary=True)
    LOGGER.info("Mesh written to '%s' (%d points, %d cells).", output, grid.n_points, grid.n_cells)
    return output


def main(argv: list[str] | None = None) -> int:
    """Command line entry point: ``pytec3d-convert-cdb FILE.cdb``."""
    parser = argparse.ArgumentParser(description="Convert an ANSYS .cdb mesh to a pytec3d .vtu mesh.")
    parser.add_argument("cdb_file", help="The .cdb file; the .vtu file is written next to it.")
    arguments = parser.parse_args(argv)
    try:
        convert_cdb(arguments.cdb_file)
    except (OSError, ValueError, ImportError) as exc:
        LOGGER.error("Conversion failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
