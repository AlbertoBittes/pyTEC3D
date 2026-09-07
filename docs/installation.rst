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

- ``python/apec_m``: main Python API
- ``python/apec_m_coordinate_system``: coordinate_system helpers
- ``python/apec_m_solver``: compiled solver extension location
- ``python/apec_m_paraview``: ParaView support module
- ``python/apec_m_algorithm``: APEC algorithm module
- ``python/apec_m_mesh_converter``: mesh conversion package

Quick Start (Local Workflow)
----------------------------

1. Clone and enter the repository.

.. code-block:: bash

   git clone <your-repo-url>
   cd PlasmaArc

2. Create a Python 3.12 development environment.

.. code-block:: bash

   conda create -n apec-m -c conda-forge \
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
     openpyxl \
     mpi4py \
     pip \
     -y
   conda activate apec-m
   ./scripts/build.sh
   python -m pip install .

Editable install:

.. code-block:: bash

   python -m pip install -e .

Fast Rebuild Loop for Solver-Only C++ Changes
----------------------------------------------

After a full initial run of ``scripts/build.sh``, rebuild with CMake:

.. code-block:: bash

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

Sanity checks
-------------

.. code-block:: bash

   python -c "import apec_m_solver.apec_m_solver"
   python -c "from apec_m_solver import APEC_M; APEC_M()"

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.

When building on Read the Docs, heavy runtime dependencies are mocked so docs build quickly. They are still required for full runtime functionality.
