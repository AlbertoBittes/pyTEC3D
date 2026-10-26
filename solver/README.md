# pytec3d-solver

The C++ solver of pyTEC3D: a core library (libMesh, PETSc / MUMPS, MPI; no Python, no VTK) that works on views of the
caller's arrays (`include/model_input.h`), and its pybind11 module `pytec3d_solver.pytec3d_solver`, driven by the
`pytec3d` Python package.

- `CMakeLists.txt`, `pyproject.toml`: the build (CMake, through scikit-build-core for the Python distribution).
- `include/`, `src/`: the core and the bindings (`src/python_bindings.cpp`).
- `tools/`: the compiler of the JSON element library (`../element_library`, package `pytec3d-element-library`) into the module.
- `libmesh/`: the libMesh build (pinned commit, patches) the solver links.
- `examples/`: a C++ program that drives the core without Python.

See `docs/solver_interface.rst` and `docs/installation.rst` of the repository.
