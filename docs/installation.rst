Installation
============

pyTEC3D has two parts that are installed together:

* the **Python package** ``pytec3d`` (``src/pytec3d``): model preparation, analysis driver, plasma algorithm,
  post-processing helpers and the JSON element library;
* the **C++ solver module** ``pytec3d_solver.pytec3d_solver``: a Python extension (pybind11) linked against libMesh,
  PETSc, MUMPS and VTK, built from ``src/pytec3d_solver`` by ``build.sh``.

The build compiles a patched copy of libMesh once, then the solver module; ``pip`` then installs both packages in the
active environment. Every dependency comes from **conda-forge**, so that the compilers, MPI, PETSc, VTK and Python
are binary compatible.

Requirements
------------

* Linux, x86_64 (the only platform tested).
* A conda installation (Miniforge or Miniconda) able to use the ``conda-forge`` channel.
* Internet access during the first build: ``build.sh`` clones libMesh from GitHub.
* About 5 GB of disk space and 30 to 60 minutes for the first build (libMesh); later solver rebuilds take a few
  minutes.
* Python 3.12 (the package requires ``>=3.12,<3.13``).

.. note::

   Do not install ``vtk`` or ``mpi4py`` with ``pip``. The solver module is linked against the VTK and the MPI of the
   conda environment; a second copy installed by pip is a different binary and makes the module fail to import or
   crash. ``pyproject.toml`` deliberately does not list ``vtk``, and lists ``mpi4py`` only so that pip finds the conda
   one already installed.

1. Get the sources
------------------

.. code-block:: bash

   git clone https://github.com/AlbertoBittes/pyTEC3D.git
   cd pyTEC3D

2. Create the conda environment
-------------------------------

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
     pyvista \
     pybind11 \
     spdlog \
     nlohmann_json \
     numpy \
     pandas \
     mpi4py \
     pip \
     -y
   conda activate pytec3d

What each package is for:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Package
     - Role
   * - ``compilers``, ``cmake``, ``ninja``, ``pkg-config``
     - C, C++ and Fortran compilers (with the MPI wrappers ``mpicc``, ``mpicxx``, ``mpifort``) and the build tools.
   * - ``openmpi``, ``mpi4py``
     - MPI, used by PETSc, libMesh and the solver; ``mpi4py`` gives the MPI rank to the Python side.
   * - ``petsc``
     - Linear algebra: the KSP solvers and preconditioners, MUMPS (direct solver), METIS / ParMETIS (partitioning).
   * - ``vtk``
     - The mesh format (``.vtu``) and all the mesh operations of the Python package; also linked by the solver.
   * - ``pyvista``
     - Not used by the package: a convenient layer on VTK for plots and post-processing in scripts and notebooks
       (the tutorials use it). ``pyvista.wrap(solver.model.mesh)`` gives a pyvista view of the model, sharing its
       arrays. Installed from conda-forge, it uses the VTK of the environment.
   * - ``pybind11``, ``spdlog``, ``nlohmann_json``
     - The Python binding, the logging and the JSON parsing of the C++ solver.
   * - ``numpy``, ``pandas``
     - Numerical work on the mesh arrays (``pandas`` for tabular post-processing).

3. Build libMesh and the solver
-------------------------------

.. code-block:: bash

   ./build.sh

``build.sh`` runs in two stages:

#. **libMesh**: clones it into ``.deps/src/libmesh``, applies the patches of ``libmesh_patches/`` (Nedelec shape
   functions on prisms and pyramids, quadratic pyramids) and installs it in ``.deps/libmesh-<method>``. This stage is
   skipped when that installation already exists.
#. **Solver module**: configures ``src/CMakeLists.txt`` in ``build/pytec3d-<method>`` (always from scratch), builds
   the module and installs it as ``src/pytec3d_solver/pytec3d_solver*.so``. The JSON files of the element library
   (``src/pytec3d/resources``) are compiled into the module at this point by ``src/tools/generate_registries.py``.
   The script finally checks that the module resolves all its shared libraries without ``LD_LIBRARY_PATH``.

Options are given as environment variables:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Variable
     - Effect
   * - ``pyTEC3D_LIBMESH_METHOD=opt|dbg|devel``
     - Build type of libMesh and of the solver: ``opt`` (default, ``-O3``), ``devel`` (``-O2 -g``, assertions on),
       ``dbg`` (``-O0 -g3``, all checks; compiler warnings are errors, ``-DPYTEC3D_WARNINGS_AS_ERRORS=OFF`` to
       disable). Each method has its own libMesh installation and build directory.
   * - ``FORCE_LIBMESH=1``
     - Delete and rebuild libMesh (clone, patch, configure, build).
   * - ``CPU_COUNT=<n>``
     - Number of parallel build jobs (default: all the cores).
   * - ``PREFIX``, ``PYTHON``
     - Installation prefix of the dependencies (default ``$CONDA_PREFIX``) and Python interpreter.
   * - ``DEPS_DIR``, ``PATCH_DIR``
     - Location of the libMesh checkout / installation (``.deps``) and of the patches.
   * - ``pyTEC3D_METHOD_CFLAGS``, ``pyTEC3D_METHOD_CXXFLAGS``, ``pyTEC3D_EXTRA_{CFLAGS,CXXFLAGS,CPPFLAGS,LDFLAGS}``
     - Override the optimization flags of the method, or add flags.

Example: a debug build on 8 cores, ``pyTEC3D_LIBMESH_METHOD=dbg CPU_COUNT=8 ./build.sh``.

4. Install the Python packages
------------------------------

.. code-block:: bash

   python -m pip install -e .        # editable: the environment uses the sources of src/ directly
   # or
   python -m pip install .           # regular: copies the packages into site-packages

The **editable** install is recommended while developing or when the sources are updated often: the environment
imports ``src/pytec3d`` and ``src/pytec3d_solver`` directly, so a ``git pull`` or a rebuild of the module takes effect
at once. With a **regular** install, ``pip install .`` must be run again after every change of the sources and after
every rebuild of the solver, otherwise the environment keeps using the old copy. ``python -c "import pytec3d;
print(pytec3d.__file__)"`` shows which copy is imported.

Optional extras:

.. code-block:: bash

   python -m pip install -e ".[cdb]"   # ansys-mapdl-reader: conversion of ANSYS .cdb meshes (pytec3d-convert-cdb)
   python -m pip install -e ".[dev]"   # pytest, ruff, mypy
   python -m pip install -e ".[notebooks]"   # matplotlib, ipykernel: the tutorial notebooks

To run the notebooks, also install Jupyter (``conda install -c conda-forge jupyterlab``) and start it from
``docs/tutorials``; they also use pyvista, part of the conda environment above.

5. Check the installation
-------------------------

.. code-block:: bash

   python -c "import pytec3d_solver.pytec3d_solver"          # the compiled module imports
   python -c "from pytec3d import Solver; Solver().close()"   # the package and the solver runtime start
   python -c "import pytec3d; print(pytec3d.__version__)"

Then run a first analysis (about 3 seconds):

.. code-block:: bash

   mkdir -p /tmp/first_run && cp tests/coil_2D_avolt/{coil_2D_avolt.vtu,model_data.json,macro.py} /tmp/first_run
   cd /tmp/first_run && python macro.py
   ls output          # coil_2D_avolt_out_1.vtu  coil_2D_avolt_out_2.vtu  coil_2D_avolt_steps.csv
   less coil_2D_avolt.log

Working on a copy keeps the files of ``tests/`` unchanged: the macros save the solved model over their input
``.vtu`` (see :doc:`examples/index`).

Running in parallel
-------------------

The solver is MPI-parallel; the script is the same in serial and in parallel:

.. code-block:: bash

   mpirun -np 4 python macro.py

Every rank runs the script; the mesh is distributed by libMesh when the model is loaded, and only rank 0 writes the
output files and, by default, prints (``Solver(root_only_output=True)``). Threads can be added per MPI process with
``Solver(num_threads=n)``. See :doc:`user_manual/solver_settings` for the solver options of large models.

Rebuilding after a change
-------------------------

* **C++ sources** (``src/pytec3d_solver/core``, ``include``): incremental rebuild and install of the module,

  .. code-block:: bash

     cmake --build build/pytec3d-<method>
     cmake --install build/pytec3d-<method>

* **Element library** (``src/pytec3d/resources/*.json``: element types, fields, domains, nonlinear laws): these files
  are compiled into the module, so the same incremental rebuild is needed. CMake regenerates the registries
  automatically when a JSON file changes. ``PYTEC3D_RESOURCES`` points the Python side to another resources
  directory, but the solver always uses the library it was built with.
* **Python sources**: nothing to do with an editable install; ``pip install .`` again with a regular one.

Building this documentation
---------------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt    # sphinx, sphinx-rtd-theme, myst-nb
   sphinx-build -b html docs docs/_build/html

The API pages are generated from the docstrings (``sphinx.ext.autosummary``); the heavy dependencies (VTK, the
solver module, MPI) are mocked, so the documentation builds without the compiled solver. The tutorial notebooks are
shown with the outputs stored in it, not executed.

Troubleshooting
---------------

``ImportError: Could not import the C++ solver module 'pytec3d_solver.pytec3d_solver'``
   The module is not built or not installed in the active environment: run ``./build.sh``, then ``pip install``
   (step 4), with the ``pytec3d`` environment active.

``unresolved shared libraries (RPATH problem)`` at the end of ``build.sh``
   The module was linked against libraries outside the environment. Rebuild from a clean shell with only the conda
   environment active (no ``LD_LIBRARY_PATH`` pointing to system MPI, PETSc or VTK).

The module imports but crashes, or ``vtkUnstructuredGrid`` type errors
   A second VTK or MPI is installed by pip in the environment: ``pip uninstall vtk mpi4py`` and reinstall them with
   conda.

Changes of the sources have no effect
   The environment imports a regular (copied) install: reinstall with ``pip install -e .``.

The element library changed but the solver rejects the new names
   The JSON files are compiled into the module: rebuild it (see above).

More output, or less
   ``pytec3d.set_log_level("DEBUG")`` or the environment variable ``PYTEC3D_LOG_LEVEL=DEBUG``; ``NO_COLOR=1``
   disables the colored labels. PETSc options such as ``-ksp_monitor`` are set in ``solver_settings["petsc_options"]``
   (see :doc:`user_manual/solver_settings`).

The whole default installation in one block
-------------------------------------------

Steps 1 to 5 with the default options (optimized build, editable install, the notebook extras), to copy and paste
into a shell where ``conda`` is available. The first run takes 30 to 60 minutes (libMesh).

.. code-block:: bash

   git clone https://github.com/AlbertoBittes/pyTEC3D.git
   cd pyTEC3D
   conda create -n pytec3d -c conda-forge -y \
     python=3.12 compilers git cmake ninja pkg-config openmpi petsc vtk pyvista pybind11 spdlog nlohmann_json \
     numpy pandas mpi4py pip jupyterlab
   conda activate pytec3d
   ./build.sh
   python -m pip install -e ".[notebooks]"
   python -c "from pytec3d import Solver; Solver().close()" && echo "pyTEC3D is installed"
