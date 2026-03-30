# APEC-M local conda build

This repository uses a local conda recipe to build `apec_m` from source.

## Intended workflow

From a cloned repository:

```bash
conda create -n apec-m-build python=3.12 conda-build -y
conda activate apec-m-build
conda build conda-recipe/recipe && conda install --use-local -y apec_m
```

## What the recipe does

The recipe build script:

1. clones and builds `libmesh` from source into the conda build prefix
2. compiles the `apec_m_solver` pybind module against that libmesh/PETSc/VTK stack
3. installs the Python API and the compiled `apec_m_solver` extension into site-packages
4. installs `openpyxl` in the environment
5. vendors `openpyxl` into `python/apec_m/local_packages` during recipe build for ParaView-side script usage

After `conda install --use-local apec_m`, the package is ready to use immediately.
