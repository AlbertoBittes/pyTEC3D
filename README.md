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
  python -m pip install .
  ```

  Editable install:

  ```bash
  python -m pip install -e .
  ```

## Fast rebuild loop for solver-only C++ changes

After a full initial run of `build.sh`, rebuild with CMake:

```bash
cmake -S src -B build/pytec3d -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DLIBMESH_DIR="$CONDA_PREFIX" \
  -DLIBMESH_METHOD="${pyTEC3D_LIBMESH_METHOD:-opt}"

cmake --build build/pytec3d -j
```

## Extension output location

CMake builds the solver module in `build/pytec3d-<method>/`; `cmake --install build/pytec3d-<method>` copies it to
`src/pytec3d_solver/pytec3d_solver*.so`, where both regular and editable installs (`pip install -e .`) import it.
If needed, override the package root during build with `PYTEC3D_PYTHON_PACKAGE_ROOT=/path/to/python_packages`.