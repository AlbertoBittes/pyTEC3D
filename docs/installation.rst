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

- ``src/pytec3d``: the Python package (``from pytec3d import Solver, Model, selection``)

  - ``model``, ``solver``: the model (mesh + model data) and its solution;
  - ``selection``, ``mesh``, ``model_edit``, ``coordinate_system``, ``tree_gauge``: model preparation;
  - ``algorithm``: plasma equilibrium and evolution;
  - ``pytec3d_paraview_plugin``: ParaView filters and console helpers;
  - ``mesh_converter``: ANSYS ``.cdb`` to ``.vtu`` conversion (``pytec3d-convert-cdb``);
  - ``resources``: the JSON element library compiled into the solver.

- ``src/pytec3d_solver``: the C++ solver; the compiled module ``pytec3d_solver.pytec3d_solver`` is installed there.

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

After a full initial run of ``./build.sh``, rebuild the solver module incrementally and install it into
``src/pytec3d_solver`` (``<method>`` is ``opt``, ``dbg`` or ``devel``):

.. code-block:: bash

   cmake --build build/pytec3d-<method>
   cmake --install build/pytec3d-<method>

Debug builds (``pyTEC3D_LIBMESH_METHOD=dbg``) treat compiler warnings as errors
(``-DPYTEC3D_WARNINGS_AS_ERRORS=OFF`` to disable).

Sanity checks
-------------

.. code-block:: bash

   python -c "import pytec3d_solver.pytec3d_solver"
   python -c "from pytec3d import Solver; Solver().close()"

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.

When building on Read the Docs, heavy runtime dependencies are mocked so docs build quickly. They are still required for full runtime functionality.
