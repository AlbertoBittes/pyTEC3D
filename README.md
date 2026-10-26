# pyTEC3D

python Tokamak Eddy-currents Code 3D: A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platform

- Linux and Windows (x86_64), through the conda packages; building from source: Linux
- Python 3.12 (the Python of ParaView)

## Project structure

- `src/pytec3d`: the Python package (`from pytec3d import Solver, Model, selection`), pure Python (distribution
  `pytec3d`: `pip install ./element_library .` needs no compiler)
- `element_library/`: the JSON element library (distribution `pytec3d-element-library`), shared by the API, which reads
  it, and the solver, which compiles it into its module

  - `model`, `solver`: the model (mesh + model data) and its solution;
  - `selection`, `mesh`, `model_edit`, `coordinate_system`, `tree_gauge`: model preparation;
  - `algorithm`: plasma equilibrium and evolution;
  - `solver_input`: the NumPy views of the mesh arrays handed to the solver;
  - `pytec3d_paraview_plugin`: ParaView filters and console helpers;
  - `mesh_converter`: ANSYS `.cdb` to `.vtu` conversion (`pytec3d-convert-cdb`);
  - `resources`: an example of model data.

- `solver/`: the C++ solver, its own project (distribution `pytec3d-solver`, `solver/pyproject.toml`): a core library
  that depends neither on Python nor on VTK (it works on views of the caller's arrays, `solver/include/model_input.h`),
  its pybind11 module `pytec3d_solver.pytec3d_solver`, the libMesh build it links (`solver/libmesh/`) and a C++ example.
  See `docs/solver_interface.rst`.
- `conda/`: the conda packages (`recipe.yaml`, rattler-build): `pytec3d` = `pytec3d-base` (the API) + `pytec3d-solver`;
  `.github/workflows/conda.yml` builds them for Linux and Windows (`ci/windows/`: the MSYS2 build).
- `examples/`: commented example analyses (2D and 3D coils, an electromagnet, ITER equilibria and vertical
  displacements in 2D and 3D), each with its mesh, model data and `macro.py`.
- `docs/`: the documentation (Sphinx), with the tutorial notebooks in `docs/tutorials`.
- `tests/test_benchmarks.py`: the validation of the installation, which runs the benchmark repository.
- `pyTEC3D_APDL_benchmarks/` (git submodule): the benchmarks against ANSYS MAPDL and APEC, with the stored reference
  results (no ANSYS needed to validate pyTEC3D) and the tools for parametric studies.

## Installation with conda (Linux, Windows)

```bash
conda create -n pytec3d -c <channel> -c conda-forge python=3.12 pytec3d   # API + solver + MPI
conda install -c <channel> -c conda-forge pytec3d-base                     # or the Python API alone
```

The Python API alone also installs with pip: `pip install ./element_library .` (see `docs/installation.rst`).

## Quick start (building from source, Linux)

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

  `build.sh` builds libMesh and the solver, installs the package with `pip install .` (with the MPI, test,
  docs and notebook dependencies) and the solver distribution `pytec3d-solver`, installs `ansys-mapdl-reader` with
  `--no-deps` (its pins conflict with the conda packages), then validates the installation against the stored
  benchmark references. Flags:

  ```bash
  ./build.sh                           # pip install .     (regular, for users; the default)
  ./build.sh --pip=editable            # pip install -e .  (editable, for developers)
  ./build.sh --validate=no             # without the validation against the benchmarks
  ./build.sh --help
  ```
