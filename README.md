# APEC-M
A new code for electromagnetic simulation of plasma equilibria and disruptions

## Local conda recipe install (from source)

```bash
conda create -n apec-m-build python=3.12 conda-build -y
conda activate apec-m-build
git clone https://github.com/AlbertoBittes/PlasmaArc.git
cd PlasmaArc
conda build conda-recipe/recipe && conda install --use-local -y apec_m
```

The conda recipe performs all required steps:

1. installs the build/runtime dependencies in the conda build environment
2. clones and builds `libmesh` with the required configuration
3. compiles the `apec_m_solver` pybind module
4. installs Python API modules plus the compiled solver module
5. installs `openpyxl` in the conda environment
6. vendors `openpyxl` into `python/apec_m/local_packages` during recipe build for ParaView-side script usage

After install, `apec_m` and `apec_m.apec_m_solver` are immediately importable in the target conda environment.

## Notes

```bash
conda activate apec-m-build
```
Use your installed environment (for example `apec-m-build`) to run scripts/tests.
