# APEC-M Recipe Layout

Conda packaging now uses one full recipe:

- `apec-m`: `conda-recipe/recipes/apec-m`

The full recipe builds and installs:
- native solver (`apec_m_core.apec_m_solver`)
- coordinate system package (`apec_m_coordinate`)
- main Python API package (`apec_m`)
- mesh converter package (`apec_m_mesh_converter`)

The pip distribution path remains solver-free (`apec-m-python`, `apec-m-coordinate`, `apec-m-mesh-converter`).

Other folders under `conda-recipe/recipes/` are retained only as build assets used by the full recipe (for example, solver/libMesh patch files).

The standalone pip packages remain available for Python-only usage:
- `pip install apec-m-coordinate`
- `pip install apec-m-mesh-converter`
- `pip install apec-m-python`
