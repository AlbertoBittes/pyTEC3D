# pytec3d

The Python API of pyTEC3D: the model (a VTK unstructured grid and its model data), the selections, the meshes, the
analysis driver, the plasma algorithm and the post-processing (`from pytec3d import Solver, Model, selection`).

Pure Python, Python 3.12. The solver is the separate distribution `pytec3d-solver` (`../solver`); the element library
is `pytec3d-element-library` (`../element_library`). With conda, `conda install -c <channel> -c conda-forge pytec3d`
installs everything; with pip, from the root of the repository: `pip install ./element_library ./python`.

See the documentation of the repository (`docs/`).
