"""The pyTEC3D element library: the JSON definitions of the domains, fields, element types and nonlinear laws.

The single source of the element library, shared by the two other parts of pyTEC3D:

* the C++ solver (``pytec3d-solver``) compiles it into its module when it is built
  (``solver/tools/generate_registries.py``): the kernels and expressions of the element types become C++ code;
* the Python API (``pytec3d``) reads it at run time (:mod:`pytec3d.element_library`) to know the element types, their
  fields and the arrays a model needs, with or without a solver installed.

Files (in :func:`data_dir`): ``domains.json``, ``fields.json``, ``element_types.json``, ``nonlinear_laws.json``. A
solver and an API of the same element library agree on every name; the API checks it against the registries compiled
into the solver (:func:`pytec3d.validation.solver_registries`).
"""

import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

#: The files of the element library, in :func:`data_dir`.
FILES = ("domains.json", "fields.json", "element_types.json", "nonlinear_laws.json")

try:
    __version__ = version("pytec3d-element-library")
except PackageNotFoundError:  # Not installed: the source tree.
    __version__ = "0.0.0"


def data_dir() -> Path:
    """Directory of the JSON files of the element library."""
    return Path(__file__).resolve().parent / "data"


def load(name: str) -> dict[str, Any]:
    """The content of the file ``name`` of the element library (one of :data:`FILES`).

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    path = data_dir() / name
    with path.open(encoding="utf-8") as handle:
        data: dict[str, Any] = json.load(handle)
    return data


def main() -> None:
    """``python -m pytec3d_element_library``: print the directory of the files (used by the build of the solver)."""
    print(data_dir())
