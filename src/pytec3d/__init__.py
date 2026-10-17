"""pyTEC3D: low-frequency electromagnetic finite elements, a Python model driving a C++ (libMesh / PETSc) solver.

Main objects::

    from pytec3d import Solver, selection

    with Solver() as solver:
        solver.read_model("coil.vtu")              # the mesh, with the model data stored in it
        solver.model.read_json("model_data.json")  # or the model data from a JSON file
        selection.select(solver.model.mesh, "MATERIAL", 2, location="elem")
        selection.nsle(solver.model.mesh)
        selection.set_value_on_selection(solver.model.mesh, "CP_0", 1)
        solver.load_model_into_solver()
        solver.solve()
        solver.write_output()

Modules:

* :mod:`pytec3d.model` - :class:`Model`, the mesh and the model data;
* :mod:`pytec3d.solver` - :class:`Solver`, the solution of a model;
* :mod:`pytec3d.selection` - selecting points and cells, setting values on them;
* :mod:`pytec3d.mesh` - reading, writing and editing the arrays of a mesh;
* :mod:`pytec3d.model_edit` - coordinate mappings, periodic couplings, derived fields;
* :mod:`pytec3d.coordinate_system` - Cartesian, cylindrical and spherical coordinate systems;
* :mod:`pytec3d.tree_gauge` - tree-cotree gauge of edge fields;
* :mod:`pytec3d.algorithm` - :class:`Algorithm`, plasma equilibrium and evolution;
* :mod:`pytec3d.element_library` - the element types and the arrays they need;
* :mod:`pytec3d.validation` - the checks run before loading and solving;
* :mod:`pytec3d.errors` - the exceptions.
"""

from importlib.metadata import PackageNotFoundError, version

from . import element_library, mesh, model_edit, selection, tree_gauge
from ._logging import set_log_level
from .algorithm import Algorithm
from .coordinate_system import (
    CartesianCoordinateSystem,
    CoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
    XAxisCylindricalCoordinateSystem,
    YAxisCylindricalCoordinateSystem,
    ZAxisCylindricalCoordinateSystem,
)
from .errors import (
    ConvergenceError,
    MissingArrayError,
    ModelLockedError,
    PyTEC3DError,
    SolverStateError,
    ValidationError,
)
from .model import Model
from .solver import Solver

try:
    __version__ = version("pytec3d")
except PackageNotFoundError:  # Not installed: running from the source tree.
    __version__ = "0.0.0"

__all__ = [
    "Algorithm",
    "CartesianCoordinateSystem",
    "ConvergenceError",
    "CoordinateSystem",
    "CylindricalCoordinateSystem",
    "MissingArrayError",
    "Model",
    "ModelLockedError",
    "PyTEC3DError",
    "Solver",
    "SolverStateError",
    "SphericalCoordinateSystem",
    "ValidationError",
    "XAxisCylindricalCoordinateSystem",
    "YAxisCylindricalCoordinateSystem",
    "ZAxisCylindricalCoordinateSystem",
    "__version__",
    "element_library",
    "mesh",
    "model_edit",
    "selection",
    "set_log_level",
    "tree_gauge",
]
