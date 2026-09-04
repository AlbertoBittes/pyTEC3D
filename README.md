# APEC-M

A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported platform

- Linux
- Python 3.12 only

## Project structure

- `python/apec_m`: main Python API
- `python/apec_m_coordinate_system`: coordinate helpers
- `python/apec_m_solver`: compiled solver extension location
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
  ```

2. Install the Python packages from this repository.

  Editable install (recommended for development):

  ```bash
  python -m pip install -e .
  ```

  Fixed install:

  ```bash
  python -m pip install .
  ```

  You can also do this during build by setting `APEC_M_PIP_INSTALL_MODE`:

  ```bash
  APEC_M_PIP_INSTALL_MODE=editable ./scripts/build.sh
  # or
  APEC_M_PIP_INSTALL_MODE=fixed ./scripts/build.sh
  ```

The root pip install includes `apec_m`, `apec_m_coordinate_system`, and `apec_m_mesh_converter`.
Core solver is built locally into `python/apec_m_solver` and imported as `apec_m_solver.apec_m_solver`.

## Sanity checks

```bash
python -c "import apec_m_solver.apec_m_solver"
python -c "from apec_m import APEC_M; APEC_M()"
```

## No conda-build packaging path

This repository now uses only the local development/build flow above.
Conda recipe build/distribution files are intentionally removed from the active workflow.
