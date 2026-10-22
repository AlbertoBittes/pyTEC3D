# pyTEC3D

python Tokamak Eddy-currents Code 3D: A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platform

- Linux
- Python 3.12

## Project structure

- `src/pytec3d`: the Python package (`from pytec3d import Solver, Model, selection`)

  - `model`, `solver`: the model (mesh + model data) and its solution;
  - `selection`, `mesh`, `model_edit`, `coordinate_system`, `tree_gauge`: model preparation;
  - `algorithm`: plasma equilibrium and evolution;
  - `pytec3d_paraview_plugin`: ParaView filters and console helpers;
  - `mesh_converter`: ANSYS `.cdb` to `.vtu` conversion (`pytec3d-convert-cdb`);
  - `resources`: the JSON element library compiled into the solver.

- `src/pytec3d_solver`: the C++ solver; the compiled module `pytec3d_solver.pytec3d_solver` is installed there.
- `examples/`: commented example analyses (2D and 3D coils, an electromagnet, ITER equilibria and vertical
  displacements in 2D and 3D), each with its mesh, model data and `macro.py`.
- `docs/`: the documentation (Sphinx), with the tutorial notebooks in `docs/tutorials`.
- `tests/test_benchmarks.py`: the validation of the installation, which runs the benchmark repository.
- `pyTEC3D_APDL_benchmarks/` (git submodule): the benchmarks against ANSYS MAPDL and APEC, with the stored reference
  results (no ANSYS needed to validate pyTEC3D) and the tools for parametric studies.

## Quick start (local workflow)

1. Clone the repository with its benchmark submodule.

  ```bash
  git clone https://github.com/AlbertoBittes/pyTEC3D.git
  cd pyTEC3D
  git submodule update --init pyTEC3D_APDL_benchmarks
  ```

2. Create the conda environment of `environment.yml` (conda-forge; Python 3.12, PETSc 3.26, VTK 9, Open MPI 5) and
  build.

  ```bash
  conda env create -f environment.yml
  conda activate pytec3d
  ./build.sh
  ```

  `build.sh` builds libMesh and the solver module, installs the package with `pip install .` (with the test,
  docs and notebook dependencies), installs `ansys-mapdl-reader` with `--no-deps` (its pins conflict with the conda
  packages), then validates the installation against the stored benchmark references. Flags:

  ```bash
  ./build.sh                           # pip install .     (regular, for users; the default)
  ./build.sh --pip=editable            # pip install -e .  (editable, for developers)
  ./build.sh --validate=no             # without the validation against the benchmarks
  ./build.sh --help
  ```
