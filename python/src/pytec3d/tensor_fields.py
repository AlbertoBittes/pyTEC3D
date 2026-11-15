"""Tensor fields of the solution as 9-component VTK arrays, assembled from the rows the solver computes.

The solver evaluates a tensor field as its three rows, vector fields (the stress sigma as SX = sigma . e_x, SY, SZ;
the strain eps as EX, EY, EZ: element and nodal fields of the elastic element types). After each solve,
:func:`compose_tensor_fields` writes the tensors as one array each, cell and point data, 9 components in row-major
order (xx, xy, xz, yx, yy, yz, zx, zy, zz: the layout of VTK tensors, which ParaView shows as a tensor, its
components and its eigenvalues): ``SIGMA`` and ``EPS``. A tensor whose rows are not all in the mesh is skipped.
"""

from __future__ import annotations

import numpy as np
from vtkmodules.util import numpy_support as vtknp

#: Tensor array -> its row arrays.
TENSOR_FIELDS: dict[str, tuple[str, str, str]] = {"SIGMA": ("SX", "SY", "SZ"), "EPS": ("EX", "EY", "EZ")}


def _compose(data, name: str, rows: tuple[str, str, str]) -> bool:
    arrays = [data.GetArray(row) for row in rows]
    if any(a is None or a.GetNumberOfComponents() != 3 for a in arrays):
        return False
    values = np.concatenate([vtknp.vtk_to_numpy(a).reshape(-1, 3) for a in arrays], axis=1)
    tensor = vtknp.numpy_to_vtk(np.ascontiguousarray(values), deep=True)
    tensor.SetName(name)
    data.RemoveArray(name)
    data.AddArray(tensor)
    return True


def compose_tensor_fields(mesh) -> list[str]:
    """Write the tensor arrays of :data:`TENSOR_FIELDS` (cell and point data) from their rows; returns those written
    (as "<name> (cell)", "<name> (point)")."""
    written = []
    for name, rows in TENSOR_FIELDS.items():
        if _compose(mesh.GetCellData(), name, rows):
            written.append(f"{name} (cell)")
        if _compose(mesh.GetPointData(), name, rows):
            written.append(f"{name} (point)")
    return written
