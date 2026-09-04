# APEC-M Developer Guide

This guide is the only supported build workflow: local development with Python 3.12 and local solver compilation.

## 1) Create environment (Python 3.12)

```bash
./scripts/create_dev_env.sh apec-m-dev
conda activate apec-m-dev
```

## 2) Enable local modules (editable-by-source)

```bash
source ./scripts/activate_local_python.sh
```

This exposes local modules directly from `python/`:
- `apec_m`
- `apec_m_coordinate_system`
- `apec_m_solver`

## 3) Build native dependencies and solver

```bash
./scripts/build_solver_local.sh
```

Default build method is optimized (`opt`).

For debug symbols:

```bash
export APEC_M_LIBMESH_METHOD=dbg
./scripts/build_solver_local.sh
```

## 4) Mesh converter package (separate pip package)

Mesh conversion remains separate because of `ansys-mapdl-reader`.

Install editable:

```bash
python -m pip install -e ./python/apec_m_mesh_converter_project
```

## 5) Fast rebuild loop for solver-only C++ changes

After a full initial run of `scripts/build_solver_local.sh`, rebuild with CMake:

```bash
cmake -S src -B build/apec_m -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DLIBMESH_DIR="$CONDA_PREFIX" \
  -DLIBMESH_METHOD="${APEC_M_LIBMESH_METHOD:-dbg}" \
  -DAPEC_M_PYTHON_SITEARCH="$PWD/python" \
  -DAPEC_M_BUILD_PYBIND_MODULE=ON
cmake --build build/apec_m -j
cmake --install build/apec_m
```

## 6) Sanity checks

```bash
python -c "import apec_m_solver.apec_m_solver"
python -c "from apec_m import APEC_M; APEC_M()"
```

## 7) Notes

- Use Python 3.12 across this project.
- Conda-build recipe packaging is removed from the active developer workflow.
