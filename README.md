# APEC-M
A new code for electromagnetic simulation of plasma equilibria and disruptions

## Local conda recipe install (from source)

```bash
conda create -n apec-m -c conda-forge python=3.12 conda-build compilers -y
conda activate apec-m
git clone https://github.com/AlbertoBittes/PlasmaArc.git
cd PlasmaArc
export APEC_M_LIBMESH_METHOD=opt   # opt|dbg|devel
conda build conda-recipe/recipe
conda install --use-local -y apec_m
```

The conda recipe performs all required steps:

1. installs the build/runtime dependencies in the conda build environment
2. clones and builds `libmesh` with the required configuration
3. compiles the `apec_m_solver` pybind module
4. installs Python API modules plus the compiled solver module
5. installs `openpyxl` in the conda environment
6. vendors `openpyxl` into `python/apec_m/local_packages` during recipe build for ParaView-side script usage

After install, `apec_m` and `apec_m.apec_m_solver` are immediately importable in the target conda environment.

## Editable vs non-editable

- `conda build` flow is non-editable (`copy` mode) and produces a relocatable package.
- For local editable debugging, run:

```bash
conda create -n apec-m-dev -c conda-forge python=3.12 compilers cmake ninja pkg-config git petsc "vtk>=9.5" pybind11 mpi4py pip -y
conda activate apec-m-dev
cd PlasmaArc
export APEC_M_LIBMESH_METHOD=dbg
export APEC_M_INSTALL_MODE=editable
bash conda-recipe/recipe/build.sh
```

This uses `pip install -e .` so imports track your working tree.

## Notes

```bash
conda activate apec-m
```
Use your installed environment (for example `apec-m`) to run scripts/tests.
