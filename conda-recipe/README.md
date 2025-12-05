# libmesh-conda-recipe

This project provides a conda recipe for building and installing LibMesh, a finite element library designed for the solution of partial differential equations.

## Overview

LibMesh is a powerful library that facilitates the development of finite element applications. This conda recipe simplifies the installation process by managing dependencies and automating the build steps.

## Structure

The project contains the following key components:

- **recipe/**: This directory includes the conda recipe files.
  - **meta.yaml**: Contains metadata for the conda package, including name, version, and dependencies (VTK and PETSc).
  - **build.sh**: A shell script that automates the build process for LibMesh.
  - **README.md**: Documentation for the conda recipe.

- **README.md**: General information about the project and instructions for setting up the conda environment.

## Conda-driven installation workflow

The steps below assume you are working inside a conda environment (Python >= 3.12) and want to build everything from source with proper linkage.

1. Create and activate a fresh environment with the needed toolchain and libraries:
   ```bash
   conda create -n apec-m-env python=3.12 cmake ninja pkg-config git openmpi petsc vtk boost -y
   conda activate apec-m-env
   ```
2. Build LibMesh from source using the recipe:
   ```bash
   conda build recipe
   conda install --use-local libmesh
   ```
3. Clone the PlasmaArc repository (this project) if you have not already, so you can build the `apec_m` solver linked to the LibMesh you just installed:
   ```bash
   git clone https://github.com/AlbertoBittes/PlasmaArc.git
   cd PlasmaArc
   ```
4. Build and install the `apec_m` solver with your conda toolchain and LibMesh installation:
   ```bash
   mkdir -p build && cd build
   cmake ../src -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DLIBMESH_DIR=$CONDA_PREFIX -DPETSC_DIR=$CONDA_PREFIX -DVTK_DIR=$CONDA_PREFIX/lib/cmake/vtk -DBOOST_ROOT=$CONDA_PREFIX -DBOOST_INCLUDEDIR=$CONDA_PREFIX/include -DBOOST_LIBRARYDIR=$CONDA_PREFIX/lib
   ninja
   install -Dm755 ../bin/apec_m $CONDA_PREFIX/bin/apec_m
   ```
5. Install the Python API (requires local `openpyxl` vendored in `api/local_packages`):
   ```bash
   cd ..
   python -m pip install --no-build-isolation -e .
   ```
   Ensure `api/local_packages` is on `PYTHONPATH` so `openpyxl` can be imported when running the API.

## Additional Information

For more detailed instructions on building LibMesh, refer to the `recipe/README.md` file.
