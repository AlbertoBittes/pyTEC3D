# pyTEC3D

A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platform

- Linux
- Python 3.12

## Project structure

- `python/pytec3d`: main Python API
- `python/pytec3d_coordinate_system`: coordinate_system helpers
- `python/pytec3d_solver`: compiled solver extension location
- `python/pytec3d_paraview`: ParaView support module
- `python/pytec3d_algorithm`: pyTEC3D algorith module
- `python/pytec3d_mesh_converter`: mesh conversion package

## Quick start (local workflow)

1. Clone and enter the repository.

  ```bash
  git clone <your-repo-url>
  cd PlasmaArc
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

The compiled extension modules are written to the Python package tree (default: `python/`), not kept as import targets in `build/`.

- `pytec3d_solver` extension output: `python/pytec3d_solver`
- `pytec3d_coordinate_system` extension output: `python/pytec3d_coordinate_system`

This keeps imports consistent for both regular installs and editable installs (`pip install -e .`).
If needed, override the package root during build with `PYTEC3D_PYTHON_PACKAGE_ROOT=/path/to/python_packages`.