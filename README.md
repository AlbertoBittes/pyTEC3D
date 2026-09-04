# APEC-M
A finite-element code for electromagnetic simulation of plasma equilibria and disruptions.

## Supported Platform
- Linux only (current packaging target)

## Installation Model
APEC-M uses two distribution paths:
- conda package `apec-m`: full package with solver (`apec_m_core.apec_m_solver`) + Python modules (`apec_m`, `apec_m_coordinate`, `apec_m_mesh_converter`)
- pip Python suite:
  - `apec-m-mesh-converter`
  - `apec-m-python` (depends on `apec-m-coordinate` and `apec-m-mesh-converter`)

## Install Solver (Conda)

```bash
conda create -n apec-m -c conda-forge python=3.12 -y
conda activate apec-m
conda install -y -c <your-channel> -c conda-forge apec-m
```

If using locally built conda artifacts:

```bash
conda create -n apec-m -c conda-forge python=3.12 -y
conda activate apec-m
conda install --use-local -y apec-m
```

## Install Python Suite With Pip (No Solver)
For pip-only usage (without conda solver), install the full Python suite:

```bash
python -m pip install apec-m-python
```

Install components independently:

```bash
python -m pip install apec-m-coordinate
python -m pip install apec-m-mesh-converter
```

Local source installs for development:

```bash
python -m pip install ./python/apec_m_coordinate_project
python -m pip install ./python/apec_m_mesh_converter_project
python -m pip install ./python/apec_m_python_project
```

## What Gets Installed
Conda `apec-m` installs:
- `apec_m_core.apec_m_solver`
- `apec_m_coordinate`
- `apec_m`
- `apec_m_mesh_converter`
- solver-native runtime dependencies from conda recipe

Pip `apec-m-python` installs:
- `apec_m`
- `apec-m-coordinate`
- `apec-m-mesh-converter`

## Build And Release Model
- Conda recipe source of truth: `conda-recipe/recipes/apec-m/meta.yaml`
- Solver build script: `conda-recipe/recipes/apec-m-core/build.sh`
- CI workflow: `.github/workflows/conda-package.yml`

## Notes
- Solver distribution is conda-first.
- Python suite distribution is pip-first.
- Windows packaging is currently disabled.
