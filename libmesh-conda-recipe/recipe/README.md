# libmesh-conda-recipe

This directory contains the conda recipe for building and installing LibMesh, a finite element library designed for the solution of partial differential equations.

## Overview

LibMesh is a powerful library that provides a flexible framework for implementing finite element methods. This conda recipe automates the process of building LibMesh with its dependencies, specifically VTK and PETSc.

## Installation Instructions

To build and install LibMesh using this conda recipe, follow these steps:

1. Clone the LibMesh repository or download the release tarball.
2. Navigate to the `libmesh-conda-recipe/recipe` directory.
3. Use the conda build command to create the package:

   ```
   conda build .
   ```

4. After the build process completes, you can install the package using:

   ```
   conda install --use-local libmesh
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