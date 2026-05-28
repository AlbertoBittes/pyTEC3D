# APEC-M local conda build

This repository uses a local conda recipe to build `apec_m` from source.

## Supported modes

There are two install modes:

- `copy` (default): package files are installed into the environment (normal conda package flow).
- `editable`: package points to your working tree (`pip install -e .`) for local debugging.

There are three libMesh link methods:

- `opt`
- `dbg`
- `devel`


## Difference between 1), 2), and 3)

1. `conda build` package build (`copy` mode):
- conda creates an isolated build environment and installs dependencies from `meta.yaml`
- runs `conda-recipe/recipe/build.sh` inside that build environment
- produces a conda package (`apec_m`) to install with `conda install --use-local`
- install is non-editable (files are installed into the target environment)

2. local `build.sh` non-editable (`copy` mode):
- you run `conda-recipe/recipe/build.sh` directly in your currently activated environment
- you must pre-install build dependencies in that environment
- does not produce a conda package; it installs directly into that environment
- install is non-editable (copied into environment)

3. local `build.sh` editable (`editable` mode):
- same as 2) for dependency responsibility (you pre-install them)
- does not produce a conda package
- uses `pip install -e .`, so Python imports point to your repo working tree
- intended for active development/debugging

In short:
- use 1) for reproducible packaging/distribution
- use 2) for quick local non-editable installation
- use 3) for local development with live code edits


## 1) Conda package build (`copy` mode, recommended)

If you are starting from a clean machine, create/activate an environment with `git` first, then clone APEC-M and build from the repository root:

```bash
conda create -n apec-m -c conda-forge python=3.12 conda-build compilers git -y
conda activate apec-m

git clone https://github.com/AlbertoBittes/APEC-M.git
cd APEC-M

export APEC_M_LIBMESH_METHOD=opt   # or dbg / devel
export APEC_M_INSTALL_MODE=copy
conda build conda-recipe/recipe
conda install --use-local -y apec_m
```

Notes:
- this recipe currently uses a local source path (`source: path: ../..`), so `conda build` must be run from inside a local APEC-M clone.
- `APEC_M_INSTALL_MODE` must stay `copy` for `conda build`.
- `editable` is intentionally blocked during `conda build`.

## 2) Local non-editable build via `build.sh` (`copy` mode)

Use this if you want to run the same script outside `conda build`.

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 compilers cmake ninja pkg-config git petsc "vtk>=9.5" pybind11 mpi4py pip -y
conda activate apec-m-dev

git clone https://github.com/AlbertoBittes/APEC-M.git
cd APEC-M

export APEC_M_LIBMESH_METHOD=opt   # or dbg / devel
export APEC_M_INSTALL_MODE=copy
bash conda-recipe/recipe/build.sh
```

## 3) Local editable debug build via `build.sh` (`editable` mode)

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 compilers cmake ninja pkg-config git petsc "vtk>=9.5" pybind11 mpi4py pip -y
conda activate apec-m-dev

git clone https://github.com/AlbertoBittes/APEC-M.git
cd APEC-M

export APEC_M_LIBMESH_METHOD=dbg   # typical for debugging
export APEC_M_INSTALL_MODE=editable
bash conda-recipe/recipe/build.sh
```

This installs `apec_m` as editable (`pip install -e .`), so Python imports follow your local repo changes.

## What `build.sh` does

1. clones and builds `libmesh` from source into the target prefix
2. compiles `apec_m_solver` with CMake/Ninja, linking selected `libmesh_<method>` + `timpi_<method>`
3. installs Python API (`copy` or `editable` mode)
4. installs `ansys-mapdl-reader==0.55.2` via pip
5. vendors `openpyxl` into `python/apec_m/local_packages` for ParaView-side usage
