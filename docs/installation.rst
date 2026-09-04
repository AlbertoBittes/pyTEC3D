Installation
============

Requirements
------------
- Linux
- conda
- Python 3.12

Clone repository
----------------

.. code-block:: bash

   git clone <your-repo-url>
   cd PlasmaArc

Create environment
------------------

Use the helper script:

.. code-block:: bash

   ./scripts/create_dev_env.sh apec-m-dev
   conda activate apec-m-dev

or create it manually:

.. code-block:: bash

   conda create -n apec-m-dev -c conda-forge python=3.12 compilers git cmake ninja pkg-config openmpi petsc vtk pybind11 spdlog numpy pandas openpyxl mpi4py pip -y
   conda activate apec-m-dev

Enable local Python modules
---------------------------

.. code-block:: bash

   source ./scripts/activate_local_python.sh

Build libMesh and solver
------------------------

Run from repository root:

.. code-block:: bash

   ./scripts/build_solver_local.sh

Optional debug build:

.. code-block:: bash

   export APEC_M_LIBMESH_METHOD=dbg
   ./scripts/build_solver_local.sh

Install mesh converter (separate pip package)
---------------------------------------------

The mesh converter depends on ``ansys-mapdl-reader`` and remains separate.

.. code-block:: bash

   python -m pip install -e ./python/apec_m_mesh_converter_project

Sanity checks
-------------

.. code-block:: bash

   python -c "import apec_m_solver.apec_m_solver"
   python -c "from apec_m import APEC_M; APEC_M()"

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.

When building on Read the Docs, heavy runtime dependencies are mocked so docs build quickly. They are still required for full runtime functionality.
