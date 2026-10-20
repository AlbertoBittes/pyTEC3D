# pyTEC3D

A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

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

## Quick start (local workflow)

1. Clone and enter the repository.

  ```bash
  git clone https://github.com/AlbertoBittes/pyTEC3D.git
  cd pyTEC3D
  ```

2. Create a Python 3.12 development environment.

  ```bash
 conda create -n pytec3d -c conda-forge \
                        python=3.12 \
                        compilers \
                        git \
                        cmake \
                        ninja \
                        pkg-config \
                        openmpi \
                        petsc \
                        vtk \
                        pyvista \
                        pybind11 \
                        spdlog \
                        numpy \
                        pandas \
                        mpi4py \
                        pip \
                        -y
  conda activate pytec3d
  ./build.sh
  ```

  `build.sh` builds libMesh and the solver module, installs the package with `pip install -e .`, then installs
  `ansys-mapdl-reader` with `--no-deps` (its pins conflict with the conda packages). Flags:

  ```bash
  ./build.sh --pip=regular             # pip install .   instead of   pip install -e .
  ./build.sh --pip=none                # no pip install (run it by hand afterwards)
  ./build.sh --extras=dev,notebooks    # also the optional dependencies
  ./build.sh --help
  ```
