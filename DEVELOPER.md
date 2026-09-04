# APEC-M Developer Settings

This guide is for local development (editable workflow), not package publishing.

## Supported Platform
- Linux

## 1) Create Developer Environment

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 conda-build boa compilers git cmake ninja pkg-config openmpi petsc vtk pybind11 spdlog -y
conda activate apec-m-dev
```

## 2) Core Developer Variables

```bash
export APEC_M_LIBMESH_METHOD=opt
export PREFIX="$CONDA_PREFIX"
export LIBMESH_DIR="$CONDA_PREFIX"
```

Optional for debug builds:

```bash
export APEC_M_LIBMESH_METHOD=dbg
```

## 3) Initial Local Build (Conda Recipe Flow)

Run from repository root:

```bash
conda mambabuild conda-recipe/recipes/apec-m-coordinate -c conda-forge
conda mambabuild conda-recipe/recipes/apec-m-python -c conda-forge --use-local
conda mambabuild conda-recipe/recipes/apec-m-core -c conda-forge --use-local
conda mambabuild conda-recipe/recipes/apec-m -c conda-forge --use-local
conda install -y -c conda-forge --use-local apec-m
```

If `mambabuild` is unavailable, replace `conda mambabuild` with `conda build`.

## 4) Editable Python Development

Use source directly:

```bash
export PYTHONPATH="$PWD/python:$PYTHONPATH"
```

Or install editable Python API package:

```bash
python -m pip install -e python/apec_m_python_project
```

## 5) Fast C++ Rebuild Loop (No Full Recipe Rebuild)

After editing C++ sources, rebuild solver module in place:

```bash
cmake -S src -B build/apec_m -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DLIBMESH_DIR="$CONDA_PREFIX" \
  -DLIBMESH_METHOD="$APEC_M_LIBMESH_METHOD" \
  -DAPEC_M_BUILD_PYBIND_MODULE=ON
cmake --build build/apec_m -j
```

If needed, install the rebuilt module into the active environment:

```bash
cmake --install build/apec_m
```

## 6) Quick Sanity Checks

```bash
python -c "import apec_m_core.apec_m_solver"
python -c "from apec_m import APEC_M; APEC_M()"
```

## 7) When To Rebuild Full Core Recipe

Run full `apec-m-core` recipe build again if you changed:
- libMesh patch logic in `conda-recipe/recipes/apec-m-core/build.sh`
- core packaging metadata
- toolchain/dependency assumptions
