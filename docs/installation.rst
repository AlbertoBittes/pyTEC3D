Installation
============

Requirements
------------
- Linux
- conda
- Python 3.12

Clone Repository
----------------

.. code-block:: bash

   git clone <your-repo-url>
   cd PlasmaArc

Project Structure
-----------------

- ``python/pytec3d``: main Python API
- ``python/pytec3d_coordinate_system``: coordinate_system helpers
- ``python/pytec3d_solver``: compiled solver extension location
- ``python/pytec3d_paraview``: ParaView support module
- ``python/pytec3d_algorithm``: pyTEC3D algorithm module
- ``python/pytec3d_mesh_converter``: mesh conversion package

Quick Start (Local Workflow)
----------------------------

1. Clone and enter the repository.

.. code-block:: bash

   git clone <your-repo-url>
   cd PlasmaArc

2. Create a Python 3.12 development environment.

.. code-block:: bash

   conda create -n pytec3d -c conda-forge \
     python=3.12 \
     compilers \
     git \
     cmake \
     ninja \
     pkg-config \
     openmpi \
     petsc \
     vtk \
     pybind11 \
     spdlog \
     numpy \
     pandas \
     mpi4py \
     pip \
     -y
   conda activate pytec3d
   ./build.sh
   python -m pip install .

Editable install:

.. code-block:: bash

   python -m pip install -e .

Fast Rebuild Loop for Solver-Only C++ Changes
----------------------------------------------

After a full initial run of ``./build.sh``, rebuild with CMake:

.. code-block:: bash

   cmake -S src -B build/pytec3d -G Ninja \
     -DCMAKE_BUILD_TYPE=Debug \
     -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
     -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
     -DLIBMESH_DIR="$CONDA_PREFIX" \
     -DLIBMESH_METHOD="${pyTEC3D_LIBMESH_METHOD:-dbg}" \
     -DpyTEC3D_PYTHON_SITEARCH="$PWD/python" \
     -DpyTEC3D_BUILD_PYBIND_MODULE=ON
   cmake --build build/pytec3d -j
   cmake --install build/pytec3d

Sanity checks
-------------

.. code-block:: bash

   python -c "import pytec3d_solver.pytec3d_solver"
   python -c "from pytec3d_solver import pyTEC3D; pyTEC3D()"

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.

When building on Read the Docs, heavy runtime dependencies are mocked so docs build quickly. They are still required for full runtime functionality.
