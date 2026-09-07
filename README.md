# APEC-M

A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platform

- Linux
- Python 3.12

## Project structure

- `python/apec_m`: main Python API
- `python/apec_m_coordinate_system`: coordinate_system helpers
- `python/apec_m_solver`: compiled solver extension location
- `python/apec_m_paraview`: ParaView support module
- `python/apec_m_algorithm`: APEC algorith module
- `python/apec_m_mesh_converter`: mesh conversion package

## Quick start (local workflow)

1. Clone and enter the repository.

  ```bash
  git clone <your-repo-url>
  cd PlasmaArc
  ```

2. Create a Python 3.12 development environment.

  ```bash
 conda create -n apec-m -c conda-forge \
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
                        openpyxl \
                        mpi4py \
                        pip \
                        -y
  conda activate apec-m
  ./scripts/build.sh
  python -m pip install .
  ```

  Editable install:

  ```bash
  python -m pip install -e .
  ```

## Fast rebuild loop for solver-only C++ changes

After a full initial run of `scripts/build.sh`, rebuild with CMake:

```bash
cmake -S src -B build/apec_m -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DLIBMESH_DIR="$CONDA_PREFIX" \
  -DLIBMESH_METHOD="${APEC_M_LIBMESH_METHOD:-opt}"

cmake --build build/apec_m -j
```