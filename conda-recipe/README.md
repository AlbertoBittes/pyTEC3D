# APEC-M Conda Packaging

This repository uses `conda-build` as the primary distribution flow.

## Current Scope
- Linux only
- Single solver-focused conda recipe:
	- `apec-m` (builds and installs `apec_m_core.apec_m_solver`, `apec_m_coordinate`, `apec_m`, and `apec_m_mesh_converter`)
- Python suite is pip-distributed:
	- `apec-m-coordinate`
	- `apec-m-mesh-converter`
	- `apec-m-python` (depends on both)

## Source Of Truth
- Full recipe: `conda-recipe/recipes/apec-m`
- Solver build implementation reused by full recipe: `conda-recipe/recipes/apec-m-core/build.sh`

## Build Release Artifacts Once (Maintainer)
From repository root:

```bash
conda create -n apec-m-build -c conda-forge python=3.12 conda-build boa compilers git conda-index -y
conda activate apec-m-build
export APEC_M_LIBMESH_METHOD=opt
conda mambabuild conda-recipe/recipes/apec-m -c conda-forge --python 3.12
```

## Experimental: Build Solver With MSYS2 + pacman (Windows)

This is an experimental native build flow for the solver on Windows using MSYS2.
It does not produce a conda package.

Open the MSYS2 UCRT64 shell and run:

```bash
pacman -Syu --noconfirm
pacman -S --needed --noconfirm \
	git make pkgconf autoconf automake libtool \
	mingw-w64-ucrt-x86_64-toolchain \
	mingw-w64-ucrt-x86_64-cmake \
	mingw-w64-ucrt-x86_64-ninja \
	mingw-w64-ucrt-x86_64-openmpi \
	mingw-w64-ucrt-x86_64-petsc \
	mingw-w64-ucrt-x86_64-vtk \
	mingw-w64-ucrt-x86_64-boost \
	mingw-w64-ucrt-x86_64-spdlog \
	mingw-w64-ucrt-x86_64-python \
	mingw-w64-ucrt-x86_64-python-pip \
	mingw-w64-ucrt-x86_64-pybind11

cd /path/to/APEC-M
export PREFIX=/ucrt64
export PATH=/ucrt64/bin:$PATH
export APEC_M_LIBMESH_METHOD=opt

bash conda-recipe/recipes/apec-m-core/build.sh
```

Quick check:

```bash
python -c "import apec_m_core.apec_m_solver"
```

Collect generated artifacts into a distributable folder and index it:

```bash
mkdir -p conda-release/linux-64 conda-release/noarch
cp -a "$CONDA_PREFIX/conda-bld/linux-64"/apec-m*.tar.bz2 conda-release/linux-64/
cp -a "$CONDA_PREFIX/conda-bld/linux-64"/apec-m*.conda conda-release/linux-64/ 2>/dev/null || true
cp -a "$CONDA_PREFIX/conda-bld/noarch"/apec-m*.tar.bz2 conda-release/noarch/ 2>/dev/null || true
cp -a "$CONDA_PREFIX/conda-bld/noarch"/apec-m*.conda conda-release/noarch/ 2>/dev/null || true
conda index conda-release
```

## Install From Release Artifacts (Users)

Users do not build anything. They only install from the provided artifact folder:

```bash
conda create -n apec-m -c conda-forge python=3.12 -y
conda activate apec-m
conda install -y -c conda-forge -c file://$PWD/conda-release apec-m
```

## End-User Guide: From GitHub Tarball To Install

### One-command install (copy/paste)

```bash
set -euo pipefail
ARCHIVE=/tmp/conda-release.tar.gz
DEST=$HOME/apec-release

curl -fL https://github.com/AlbertoBittes/APEC-M/releases/latest/download/conda-release.tar.gz -o "$ARCHIVE"
rm -rf "$DEST/conda-release"
mkdir -p "$DEST"
tar -xzf "$ARCHIVE" -C "$DEST"

conda create -n apec-m -c conda-forge python=3.12 -y
conda activate apec-m
conda install -y -c conda-forge -c "file://$DEST/conda-release" apec-m
```

This requires the GitHub release asset to be named exactly `conda-release.tar.gz`.

### 1. Download the release artifact from GitHub

- Open the project Releases page.
- Download the packaged local channel archive (example name: `conda-release.tar.gz`).

### 2. Extract the archive

```bash
mkdir -p "$HOME/apec-release"
tar -xzf /path/to/conda-release.tar.gz -C "$HOME/apec-release"
```

After extraction, you should have a local channel root like:

- `$HOME/apec-release/conda-release/linux-64`
- `$HOME/apec-release/conda-release/noarch`

### 3. Create a clean environment

```bash
conda create -n apec-m -c conda-forge python=3.12 -y
conda activate apec-m
```

### 4. Install APEC-M from the extracted artifacts

```bash
conda install -y -c conda-forge -c file://$HOME/apec-release/conda-release apec-m
```

### 5. Optional: pip-only Python suite (no solver)

```bash
python -m pip install apec-m-python
```

Use this only if you want Python modules without the conda solver package.

### 6. Quick verification

```bash
python -c "import apec_m"
python -c "import apec_m_coordinate"
python -c "import apec_m_mesh_converter"
python -c "import apec_m_core.apec_m_solver"
```

### Notes
- If `mambabuild` is not available, replace each `conda mambabuild` with `conda build`.
- Keep all produced `apec-m*` artifacts together in `conda-release` before running `conda index`.
- If the artifact folder is shared on another machine, replace `$PWD/conda-release` with that absolute path.
