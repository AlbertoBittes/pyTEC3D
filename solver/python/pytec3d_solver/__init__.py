"""The compiled C++ solver of pyTEC3D (``pytec3d_solver.pytec3d_solver``), driven by :class:`pytec3d.Solver`.

Distribution ``pytec3d-solver`` (``solver/`` of the repository): the solver core, built against libMesh, PETSc and
MPI, with its pybind11 module. It does not depend on the ``pytec3d`` package, nor on VTK: it receives NumPy views of
the model arrays.
"""

import os
import sys
from importlib.metadata import PackageNotFoundError, version

# Windows: the DLLs the module needs (the MinGW runtime, and MS-MPI from the conda environment) are next to it or in
# the Library\bin of the environment; Python only searches the directory of the module and the system ones.
if sys.platform == "win32":
    for _directory in (os.path.dirname(__file__), os.path.join(sys.prefix, "Library", "bin")):
        if os.path.isdir(_directory):
            os.add_dll_directory(_directory)

try:
    __version__ = version("pytec3d-solver")
except PackageNotFoundError:  # Not installed: a build tree.
    __version__ = "0.0.0"
