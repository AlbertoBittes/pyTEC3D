# Recipe details

This recipe builds `apec_m` locally from source. It is designed for the case where `libmesh` is not distributed as a prebuilt conda package.

## Build sequence

`build.sh` performs:

1. clone `libmesh` into the build workspace
2. configure and compile `libmesh` against conda-provided MPI/PETSc/VTK/Boost
3. install `libmesh` into the conda build prefix
4. configure and compile `apec_m_solver` (pybind module) with CMake/Ninja
5. install Python API modules into the conda package prefix with pip (`--no-deps`)
6. ensure `openpyxl` is available in the conda environment
7. vendor `openpyxl` into `python/apec_m/local_packages` during recipe build for ParaView-side script usage

## Local usage

From the repository root:

```bash
conda create -n apec-m-build python=3.12 conda-build -y
conda activate apec-m-build
conda build conda-recipe/recipe && conda install --use-local -y apec_m
```

The resulting environment contains `apec_m` + `apec_m_solver` ready to import.
