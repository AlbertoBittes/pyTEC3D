# pyTEC3D

python Tokamak Eddy-currents Code 3D: A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platforms

- Linux and Windows (x86_64), through the conda packages
- Python 3.12 (the Python of ParaView)

## Project structure

Three independent distributions, built from this repository:

- `python/`: the Python API (distribution `pytec3d`, `python/pyproject.toml`), pure Python
  (`from pytec3d import Solver, Model, selection`). In `python/src/pytec3d`:

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
- `element_library/`: the JSON element library (distribution `pytec3d-element-library`), shared by the API, which reads
  it, and the solver, which compiles it into its module.

And around them:

- `conda/`: the conda packages (`recipe.yaml`, rattler-build): `pytec3d` = `pytec3d-base` (the API) + `pytec3d-solver`
  (+ `pytec3d-element-library`, `pytec3d-libmesh`), with the build scripts of the solver (`build_windows.sh`: the MSYS2
  build of the Windows solver); `.github/workflows/conda.yml` builds them for Linux and Windows.
- `examples/`: commented example analyses (2D and 3D coils, an electromagnet, ITER equilibria and vertical
  displacements in 2D and 3D), each with its mesh, model data and `macro.py`.
- `docs/`: the documentation (Sphinx), with the tutorial notebooks in `docs/tutorials`.
- `tests/test_benchmarks.py`: the validation of an installation, which runs the benchmark repository.
- `pyTEC3D_APDL_benchmarks/` (git submodule): the benchmarks against ANSYS MAPDL and APEC, with the stored reference
  results (no ANSYS needed to validate pyTEC3D) and the tools for parametric studies.
- `environment.yml`: the conda environment to develop pyTEC3D (`docs/installation.rst`, "Development").

## Installation (conda, Linux and Windows)

```bash
conda create -n pytec3d -c <channel> -c conda-forge python=3.12 pytec3d   # API + solver + MPI
conda install -c <channel> -c conda-forge pytec3d-base                     # or the Python API alone
```

The Python API alone also installs with pip, from the root of the repository: `pip install ./element_library ./python`.
The `.cdb` converter needs `pip install --no-deps ansys-mapdl-reader appdirs tqdm` after pytec3d.

## Development

```bash
git clone https://github.com/AlbertoBittes/pyTEC3D.git && cd pyTEC3D
git submodule update --init pyTEC3D_APDL_benchmarks
conda env create -f environment.yml && conda activate pytec3d-dev
pip install --no-deps -e ./element_library -e ./python       # the API from the sources
rattler-build build -r conda/recipe.yaml -c conda-forge       # the packages, solver included, into output/
conda install -c ./output -c conda-forge pytec3d-solver       # the solver just built
python pyTEC3D_APDL_benchmarks/bench.py validate              # the validation against the benchmarks
```
