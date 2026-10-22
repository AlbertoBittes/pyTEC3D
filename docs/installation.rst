Installation
============

pyTEC3D has two parts that are installed together:

* the **Python package** ``pytec3d`` (``src/pytec3d``): model preparation, analysis driver, plasma algorithm,
  post-processing helpers and the JSON element library;
* the **C++ solver module** ``pytec3d_solver.pytec3d_solver``: a Python extension (`pybind11
  <https://pybind11.readthedocs.io/>`_) linked against `libMesh <https://libmesh.github.io/>`_, `PETSc
  <https://petsc.org/>`_, `MUMPS <https://mumps-solver.org/>`_ and `VTK <https://vtk.org/>`_.

One script, ``build.sh``, compiles a patched copy of libMesh (once), builds the solver module, installs both packages
with ``pip`` in the active conda environment, checks the installation and validates it against the benchmark
references. Every dependency comes from `conda-forge <https://conda-forge.org/>`_, so that the compilers, MPI,
PETSc, VTK and Python are binary compatible.

Requirements
------------

* Linux, x86_64 (the only platform tested).
* A conda installation (`Miniforge <https://github.com/conda-forge/miniforge>`_ or Miniconda) able to use the
  ``conda-forge`` channel.
* Internet access during the first build: ``build.sh`` clones libMesh from GitHub.
* About 5 GB of disk space and 30 to 60 minutes for the first build (libMesh); later builds take a few minutes.
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
   git submodule update --init pyTEC3D_APDL_benchmarks    # the benchmarks that validate the build

2. Create the conda environment
-------------------------------

The environment is defined by ``environment.yml`` at the root of the repository (the conda-forge packages, with
Python 3.12, PETSc 3.26, VTK 9 and Open MPI 5 pinned):

.. code-block:: bash

   conda env create -f environment.yml
   conda activate pytec3d

After an update of the sources that changes ``environment.yml``, update the environment with ``conda env update -f
environment.yml`` before building again.

What each package is for:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Package
     - Role
   * - ``compilers``, ``cmake``, ``ninja``, ``pkg-config``
     - C, C++ and Fortran compilers (with the MPI wrappers ``mpicc``, ``mpicxx``, ``mpifort``) and the build tools.
   * - ``openmpi``, ``mpi4py``
     - MPI (`Open MPI <https://www.open-mpi.org/>`_), used by PETSc, libMesh and the solver; `mpi4py
       <https://mpi4py.readthedocs.io/>`_ gives the MPI rank to the Python side.
   * - ``petsc``
     - Linear algebra: the KSP solvers and preconditioners, MUMPS (direct solver), METIS / ParMETIS (partitioning).
   * - ``vtk``
     - The mesh format (``.vtu``) and all the mesh operations of the Python package; also linked by the solver.
   * - ``pyvista``
     - Not used by the package: a convenient layer on VTK for plots and post-processing in scripts and notebooks
       (the tutorials use it). ``pyvista.wrap(solver.model.mesh)`` gives a `pyvista <https://docs.pyvista.org/>`_
       view of the model, sharing its arrays. Installed from conda-forge, it uses the VTK of the environment.
   * - ``pybind11``, ``spdlog``, ``nlohmann_json``
     - The Python binding, the logging and the JSON parsing of the C++ solver.
   * - ``numpy``, ``pandas``
     - Numerical work on the mesh arrays (``pandas`` for tabular post-processing).

3. Build and install
--------------------

.. code-block:: bash

   ./build.sh                    # regular install, for users (default)
   ./build.sh --pip=editable     # editable install, for developers

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Install
     - Effect
   * - ``--pip=regular`` (default)
     - ``pip install .``: the packages are copied into ``site-packages``. For users; run ``./build.sh`` again after
       updating the sources (``git pull``).
   * - ``--pip=editable``
     - ``pip install -e .``: the environment imports ``src/pytec3d`` and ``src/pytec3d_solver`` directly, so a change
       of the Python sources or a rebuild of the module takes effect at once. For developers.

``build.sh`` runs these stages, each of which stops the build with an error when it fails:

#. **libMesh**: clones it, applies the patches of ``libmesh_patches/`` (Nedelec shape functions on prisms and
   pyramids, quadratic pyramids) and installs it in ``.deps/``. Skipped when that installation already exists.
#. **Solver module**: builds ``src/pytec3d_solver/pytec3d_solver*.so``, with the JSON files of the element library
   (``src/pytec3d/resources``) compiled into it, and checks that it resolves all its shared libraries and links the
   MPI of PETSc only.
#. **Python packages**: ``pip install`` of ``pytec3d`` with the test, documentation and notebook dependencies
   (pytest, Sphinx, JupyterLab, matplotlib...), then a check that the installed copy imports.
#. **ansys-mapdl-reader**: the reader of ANSYS ``.cdb`` meshes (``pytec3d-convert-cdb``,
   :func:`pytec3d.mesh_converter.convert_cdb`), installed with ``--no-deps``: its declared requirements conflict with
   the conda packages the solver is built against, while the package works with them.
#. **Validation**: the installed pyTEC3D runs the benchmark cases and is compared with the stored ANSYS MAPDL / APEC
   results, about 10 minutes, no ANSYS needed (:doc:`testing`). ``--validate=long`` also runs the 3D ITER cases (about
   1 h more), ``--validate=no`` skips the validation.

All the cores are used, and everything is built in and for the active conda environment.

For developers
~~~~~~~~~~~~~~

* **Install** with ``./build.sh --pip=editable``.
* **After a change** of the C++ sources (``src/pytec3d_solver/core``, ``include``) or of the element library
  (``src/pytec3d/resources/*.json``, compiled into the module), run ``./build.sh --pip=editable`` again: libMesh is
  reused and only the module is rebuilt (``--validate=no`` to skip the validation while iterating). Changes of the
  Python sources need nothing with an editable install.
* **Build type**: ``pyTEC3D_LIBMESH_METHOD=opt|devel|dbg ./build.sh``. ``opt`` (default, ``-O3``), ``devel``
  (``-O2 -g``, assertions on), ``dbg`` (``-O0 -g3``, all checks, compiler warnings are errors). Each method has its
  own libMesh installation, built the first time it is used.
* **Rebuild libMesh** from scratch (new clone, patches, configure, build): ``FORCE_LIBMESH=1 ./build.sh``.

A first run
-----------

A first analysis (about 3 seconds), on a copy of an example so that the repository stays unchanged (the macro saves
the solved model over its input ``.vtu``, see :doc:`examples/index`):

.. code-block:: bash

   cp -r examples/coil_2D_avolt /tmp/first_run
   cd /tmp/first_run && python macro.py
   ls output          # coil_2D_avolt_out_1.vtu  coil_2D_avolt_out_2.vtu  coil_2D_avolt_steps.csv
   less coil_2D_avolt.log

The tutorial notebooks run with ``jupyter lab`` started from ``docs/tutorials`` (:doc:`tutorials/index`).

Running in parallel
-------------------

The solver is MPI-parallel; the script is the same in serial and in parallel:

.. code-block:: bash

   mpirun -np 4 python macro.py

Every rank runs the script; the mesh is distributed by libMesh when the model is loaded, and only rank 0 writes the
output files and, by default, prints (``Solver(root_only_output=True)``). Threads can be added per MPI process with
``Solver(num_threads=n)``. See :doc:`user_manual/solver_settings` for the solver options of large models, and the
`Open MPI documentation of mpirun <https://docs.open-mpi.org/en/main/man-openmpi/man1/mpirun.1.html>`_ for process
placement.

Building this documentation
---------------------------

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

The API pages are generated from the docstrings (`autosummary
<https://www.sphinx-doc.org/en/master/usage/extensions/autosummary.html>`_); the heavy dependencies (VTK, the solver
module, MPI) are mocked, so the documentation builds without the compiled solver (``docs/requirements.txt`` lists what
Read the Docs installs). The tutorial notebooks are shown with the outputs stored in them, not executed (`MyST-NB
<https://myst-nb.readthedocs.io/>`_).

Troubleshooting
---------------

``ImportError: Could not import the C++ solver module 'pytec3d_solver.pytec3d_solver'``
   The module is not built or not installed in the active environment: run ``./build.sh`` with the ``pytec3d``
   environment active.

``unresolved shared libraries (RPATH problem)`` at the end of the module build
   The module was linked against libraries outside the environment. Rebuild from a clean shell with only the conda
   environment active (no ``LD_LIBRARY_PATH`` pointing to system MPI, PETSc or VTK).

The module imports but crashes, or ``vtkUnstructuredGrid`` type errors
   A second VTK or MPI is installed by pip in the environment: ``pip uninstall vtk mpi4py`` and reinstall them with
   conda.

Changes of the sources have no effect
   The environment imports a regular (copied) install: run ``./build.sh`` again, or ``./build.sh --pip=editable``
   while developing.

The element library changed but the solver rejects the new names
   The JSON files are compiled into the module: run ``./build.sh`` again.

The validation fails
   The summary printed by ``build.sh`` gives the cases and quantities out of tolerance; see :doc:`testing`.

More output, or less
   ``pytec3d.set_log_level("DEBUG")`` or the environment variable ``PYTEC3D_LOG_LEVEL=DEBUG``; ``NO_COLOR=1``
   disables the colored labels. PETSc options such as ``-ksp_monitor`` are set in ``solver_settings["petsc_options"]``
   (see :doc:`user_manual/solver_settings`).

The whole installation in one block
-----------------------------------

To copy and paste into a shell where ``conda`` is available (``./build.sh --pip=editable`` for developers). The
first run takes 30 to 60 minutes (libMesh).

.. code-block:: bash

   git clone https://github.com/AlbertoBittes/pyTEC3D.git
   cd pyTEC3D
   git submodule update --init pyTEC3D_APDL_benchmarks
   conda env create -f environment.yml
   conda activate pytec3d
   ./build.sh
