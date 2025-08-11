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

## Usage

To build and install LibMesh using this conda recipe, follow these steps:

1. Clone the repository or download the release tarball of LibMesh.
2. Navigate to the `libmesh-conda-recipe` directory.
3. Use conda to build the recipe:
   ```
   conda build recipe
   ```
4. Install the built package using:
   ```
   conda install --use-local libmesh
   ```

## Additional Information

For more detailed instructions on building LibMesh, refer to the `recipe/README.md` file.