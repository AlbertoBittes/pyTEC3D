# Recipe details

This recipe builds `apec_m` locally from source. It is designed for the case where `libmesh` is not distributed as a prebuilt conda package.

## Build sequence (`build.sh`)

`build.sh` performs:

1. clone `libmesh` into the build workspace
2. configure and compile `libmesh` against conda-provided MPI/PETSc/Boost with VTK disabled
3. install `libmesh` into the conda build prefix for the selected method (`opt` or `dbg` or `devel`)
4. configure and compile `apec_m_solver` (pybind module) with CMake/Ninja, linking VTK directly and selecting one libMesh method
5. install Python API modules into the conda package prefix with pip (`--no-deps`)
6. install `ansys-mapdl-reader==0.55.2` via pip (PyPI)
7. vendor `openpyxl` into `python/apec_m/local_packages` during recipe build for ParaView-side script usage

## Environment flags

- `APEC_M_LIBMESH_METHOD=opt|dbg|devel`  
  selects which `libmesh_<method>` / `timpi_<method>` libraries are linked by `apec_m_solver`.
- `APEC_M_INSTALL_MODE=copy|editable`  
  `copy` is the default packaging mode. `editable` is for local debug runs and is blocked in `conda build`.
- optional compile flag overrides:
  - `APEC_M_METHOD_CFLAGS`
  - `APEC_M_METHOD_CXXFLAGS`
  - `APEC_M_EXTRA_CFLAGS`
  - `APEC_M_EXTRA_CXXFLAGS`
  - `APEC_M_EXTRA_CPPFLAGS`
  - `APEC_M_EXTRA_LDFLAGS`

Default method-oriented flags in `build.sh`:
- `opt`: `-O3 -g0`
- `devel`: `-O2 -g`
- `dbg`: `-O0 -g3` and `-DNDEBUG` is removed (`-UNDEBUG` enforced)

## Usage matrix

- `conda build`:
  - allowed: `APEC_M_INSTALL_MODE=copy`
  - blocked: `APEC_M_INSTALL_MODE=editable`
- local `bash conda-recipe/recipe/build.sh`:
  - allowed: both `copy` and `editable`

## Commands

From repository root.

### A) Package build (non-editable, recommended)

```bash
conda create -n apec-m -c conda-forge python=3.12 conda-build compilers -y
conda activate apec-m
export APEC_M_LIBMESH_METHOD=opt   # or dbg / devel
conda build conda-recipe/recipe
conda install --use-local -y apec_m
```

### B) Local script run (non-editable copy mode)

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 compilers cmake ninja pkg-config git petsc "vtk>=9.5" pybind11 mpi4py pip -y
conda activate apec-m-dev
export APEC_M_LIBMESH_METHOD=opt   # or dbg / devel
export APEC_M_INSTALL_MODE=copy
bash conda-recipe/recipe/build.sh
```

### C) Local script run (editable debug mode)

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 compilers cmake ninja pkg-config git petsc "vtk>=9.5" pybind11 mpi4py pip -y
conda activate apec-m-dev
export APEC_M_LIBMESH_METHOD=dbg   # or opt / devel
export APEC_M_INSTALL_MODE=editable
bash conda-recipe/recipe/build.sh
```
