# libmesh-conda-recipe

This directory contains the conda recipe for building and installing LibMesh, a finite element library designed for the solution of partial differential equations.

## Overview

LibMesh is a powerful library that provides a flexible framework for implementing finite element methods. This conda recipe automates the process of building LibMesh with its dependencies, specifically VTK and PETSc.

## Installation Instructions (Conda workflow)

1. Start in a clean conda environment with Python >= 3.12 and the required toolchain/libs:

   ```bash
   conda create -n apec-m-env python=3.12 cmake ninja pkg-config git openmpi petsc vtk boost -y
   conda activate apec-m-env
   ```

2. From this directory (`libmesh-conda-recipe/recipe`), build LibMesh with the bundled script:

   ```bash
   conda build .
   conda install --use-local libmesh
   ```

3. To compile `apec_m` linked to the newly built LibMesh, clone the PlasmaArc repository and run CMake/Ninja with the conda paths:

   ```bash
   git clone https://github.com/AlbertoBittes/PlasmaArc.git
   cd PlasmaArc
   mkdir -p build && cd build
   cmake ../src -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DLIBMESH_DIR=$CONDA_PREFIX -DPETSC_DIR=$CONDA_PREFIX -DVTK_DIR=$CONDA_PREFIX/lib/cmake/vtk -DBOOST_ROOT=$CONDA_PREFIX -DBOOST_INCLUDEDIR=$CONDA_PREFIX/include -DBOOST_LIBRARYDIR=$CONDA_PREFIX/lib
   ninja
   install -Dm755 ../bin/apec_m $CONDA_PREFIX/bin/apec_m
   ```

4. Install the Python API and ensure the vendored `openpyxl` in `api/local_packages` is discoverable:

   ```bash
   cd ..
   python -m pip install --no-build-isolation -e .
   # Add api/local_packages to PYTHONPATH if needed so openpyxl can be imported
   ```

## Build Process

The build process is automated through the `build.sh` script, which performs the following steps:

- Initializes git submodules if the source is a git clone.
- Creates a separate build directory.
- Configures the build with the specified installation prefix.
- Compiles the LibMesh library.
- Optionally runs tests and installs the library to the specified directory.

## Dependencies

This recipe requires the following dependencies:

- VTK
- PETSc

Ensure that these packages are available in your conda environment before building LibMesh.

## Additional Information

For more details on LibMesh and its capabilities, please refer to the official [LibMesh documentation](http://libmesh.github.io/).
