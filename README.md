# APEC-M
A new code for electromagnetic simulation of plasma equilibria and disruptions

## Conda-based setup (Python >= 3.12)

```bash
conda create -n apec-m-env python=3.12 cmake ninja pkg-config git openmpi petsc vtk boost -y
conda activate apec-m-env
```

Build LibMesh via the bundled recipe (clones libMesh and links to the conda MPI/PETSc/VTK/Boost):

```bash
cd conda-recipe/recipe
conda build .
conda install --use-local libmesh
```

Compile `apec_m` against that LibMesh:

```bash
cd ../../
mkdir -p build && cd build
cmake ../src -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DLIBMESH_DIR=$CONDA_PREFIX -DPETSC_DIR=$CONDA_PREFIX -DVTK_DIR=$CONDA_PREFIX/lib/cmake/vtk -DBOOST_ROOT=$CONDA_PREFIX -DBOOST_INCLUDEDIR=$CONDA_PREFIX/include -DBOOST_LIBRARYDIR=$CONDA_PREFIX/lib
ninja
install -Dm755 ../bin/apec_m $CONDA_PREFIX/bin/apec_m
```

Install the Python API (requires the vendored `openpyxl` from `api/local_packages` to be on `PYTHONPATH`):

```bash
cd ..
python -m pip install --no-build-isolation -e .
# export PYTHONPATH=$PWD/api/local_packages:$PYTHONPATH if openpyxl is not otherwise available
```
