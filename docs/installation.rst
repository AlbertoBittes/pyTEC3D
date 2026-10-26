Installation
============

pyTEC3D has three parts, independent distributions built from the same repository:

* the **Python API** ``pytec3d`` (``src/pytec3d``, ``pyproject.toml``): model preparation, analysis driver, plasma
  algorithm, post-processing helpers and the JSON element library. Pure Python: it installs anywhere, with ``numpy``
  and ``vtk``.
* the **C++ solver** ``pytec3d-solver`` (``solver/``, ``solver/pyproject.toml``), the module
  ``pytec3d_solver.pytec3d_solver``: a solver core in C++ linked against `libMesh <https://libmesh.github.io/>`_,
  `PETSc <https://petsc.org/>`_, `MUMPS <https://mumps-solver.org/>`_ and MPI, with a `pybind11
  <https://pybind11.readthedocs.io/>`_ binding. It does not use VTK nor the Python package: it receives NumPy views of
  the arrays of the mesh (:doc:`solver_interface`).
* the **element library** ``pytec3d-element-library`` (``element_library/``): the JSON definitions of the domains,
  fields, element types and nonlinear laws (:doc:`user_manual/element_library`), data only. The solver compiles it
  into its module when it is built, the API reads it at run time: both come from this single source.

Conda packages
--------------

The simplest installation, on Linux (x86_64) and Windows (x86_64): the conda packages, built by the GitHub workflow
``.github/workflows/conda.yml`` from ``conda/recipe.yaml``, and published on the channel of the project with
`conda-forge <https://conda-forge.org/>`_ for the dependencies:

.. code-block:: bash

   conda create -n pytec3d -c <channel> -c conda-forge python=3.12 pytec3d
   conda activate pytec3d

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Package
     - Platforms
     - Content
   * - ``pytec3d``
     - all
     - Everything: ``pytec3d-base`` and ``pytec3d-solver``, with their dependencies.
   * - ``pytec3d-base``
     - all (noarch)
     - The Python API (the pip distribution ``pytec3d``), with the element library, NumPy, VTK, pyvista and mpi4py.
       Without the solver, everything works except loading a model into a solver.
   * - ``pytec3d-element-library``
     - all (noarch)
     - The JSON element library. ``pytec3d-solver`` only installs next to the element library it was built with.
   * - ``pytec3d-solver``
     - linux-64, win-64
     - The solver, for Python 3.12 (the Python of ParaView). Linux: linked against the PETSc, MUMPS and Open MPI of
       conda-forge and ``pytec3d-libmesh``. Windows: PETSc, MUMPS, libMesh and the GCC runtime are linked statically
       into the module (built in MSYS2, ``ci/windows/build.sh``), shipped with the Fortran interface of MS-MPI; MS-MPI
       itself is the conda package ``msmpi``.
   * - ``pytec3d-libmesh``
     - linux-64
     - The patched libMesh the solver links: the libMesh development branch at the time of the build
       (``solver/libmesh/build_libmesh.sh``; the commit is in ``share/pytec3d/libmesh-commit``).

**MPI.** The solver depends on the MPI it is linked with (Open MPI on Linux, MS-MPI on Windows), and ``mpi4py``, a
dependency of ``pytec3d-base``, has a conda build per MPI: in an environment with the solver, conda always installs the
mpi4py of that MPI (a conda-forge environment holds only one MPI). ``conda install pytec3d-base`` alone installs the
mpi4py of any MPI; ``conda install pytec3d-solver`` installs the solver alone. In parallel:
``mpiexec -n 4 python macro.py``.

**The Python API with pip.** ``pip install ./element_library .`` (from the repository; the element library first, the
API depends on it) installs ``pytec3d`` alone. ``pip install ".[mpi]"`` adds mpi4py,
which then builds against the MPI it finds (``mpicc`` on ``PATH``), or uses an MPI installed with pip next to it, e.g.
``pip install ".[mpi]" openmpi`` (or ``mpich``, ``impi-rt``). With a solver of a conda environment, install mpi4py with
conda instead, so that it uses the MPI of the solver.

**Publishing.** The workflow builds the packages at every push and pull request (the artifacts ``conda-linux-64`` and
``conda-win-64``, a channel layout), installs them on both systems and solves a model in serial and on 2 MPI ranks
(``conda/tests/test_square_conductor.py``). On a tag ``v*`` it uploads them to the anaconda.org channel of the
repository variable ``ANACONDA_OWNER`` (with the secret ``ANACONDA_API_KEY``) and attaches them to the GitHub release.
To use the packages of a run without a channel: download the two artifacts into one directory, index it with
``rattler-index fs <dir>`` and install with ``-c file://<dir> -c conda-forge``.

Building from source (Linux)
----------------------------

For development: ``build.sh`` compiles the patched libMesh (once), builds the solver (``pip install ./solver``),
installs both distributions with ``pip`` in the active conda environment, checks the installation and validates it
against the benchmark references. Every dependency of the solver comes from conda-forge (``environment.yml``), so that
the compilers, MPI, PETSc and Python are binary compatible.

Requirements
~~~~~~~~~~~~

* Linux, x86_64.
* A conda installation (`Miniforge <https://github.com/conda-forge/miniforge>`_ or Miniconda) able to use the
  ``conda-forge`` channel.
* Internet access during the first build: ``build.sh`` clones libMesh from GitHub.
* About 5 GB of disk space and 30 to 60 minutes for the first build (libMesh); later builds take a few minutes.
* Python 3.12 or newer.

.. note::

   Do not install ``mpi4py`` with ``pip``: it must use the MPI the solver is linked against, the one of the conda
   environment (the ``mpi`` extra of ``pyproject.toml`` lists it only so that pip finds the conda one already
   installed). VTK is not linked by the solver: the ``vtk`` of conda or of PyPI both work.

1. Get the sources
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/AlbertoBittes/pyTEC3D.git
   cd pyTEC3D
   git submodule update --init pyTEC3D_APDL_benchmarks    # the benchmarks that validate the build

2. Create the conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
     - The mesh format (``.vtu``) and all the mesh operations of the Python package (not used by the solver).
   * - ``pyvista``
     - A requirement of the package: the reader of ANSYS ``.cdb`` meshes (``ansys-mapdl-reader``) needs it, and it is
       a convenient layer on VTK for plots and post-processing in scripts and notebooks (the tutorials use it).
       ``pyvista.wrap(solver.model.mesh)`` gives a `pyvista <https://docs.pyvista.org/>`_
       view of the model, sharing its arrays. Installed from conda-forge, it uses the VTK of the environment.
   * - ``pybind11``, ``spdlog``, ``nlohmann_json``
     - The Python binding, the logging and the JSON parsing of the C++ solver.
   * - ``numpy``, ``pandas``
     - Numerical work on the mesh arrays (``pandas`` for tabular post-processing).

3. Build and install
~~~~~~~~~~~~~~~~~~~~

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
     - ``pip install -e .``: the environment imports ``src/pytec3d`` directly, so a change of the Python sources takes
       effect at once. For developers.

The solver distribution ``pytec3d-solver`` is always installed as a copy (``pip install ./solver``, built in
``build/pytec3d-<method>``): run ``./build.sh`` again after a change of the C++ sources.

``build.sh`` runs these stages, each of which stops the build with an error when it fails:

#. **libMesh** (``solver/libmesh/build_libmesh.sh``, shared with the conda recipe and the Windows build): clones the
   libMesh development branch (``LIBMESH_REF=<branch, tag or commit>`` to choose another), applies the patches of ``solver/libmesh/patches`` (Nedelec shape functions on prisms and pyramids,
   quadratic pyramids) and installs it in ``.deps/``. Skipped when that installation already exists.
#. **Python API**: ``pip install`` of ``pytec3d`` with the MPI, test, documentation and notebook dependencies
   (mpi4py, pytest, Sphinx, JupyterLab, matplotlib...).
#. **Solver**: ``pip install ./solver`` (scikit-build-core and CMake, in ``build/pytec3d-<method>``) builds the solver
   core (``libpytec3d_core.a``, no Python, no VTK) with the JSON files of the element library
   (``element_library/``) compiled into it, and the module ``pytec3d_solver.cpython-3XX-*.so``; checks that the
   module resolves all its shared libraries, links the MPI of PETSc only and has the interface version of the package.
   ``PYTEC3D_CMAKE_ARGS`` passes more CMake options (e.g. ``-DPYTEC3D_BUILD_CPP_EXAMPLES=ON``,
   :doc:`solver_interface`).
#. **ansys-mapdl-reader**: the reader of ANSYS ``.cdb`` meshes (``pytec3d-convert-cdb``,
   :func:`pytec3d.mesh_converter.convert_cdb`), installed with ``--no-deps``: its declared requirements conflict with
   the conda packages the solver is built against, while the package works with them.
#. **Validation**: the installed pyTEC3D runs the benchmark cases and is compared with the stored ANSYS MAPDL / APEC
   results, about 10 minutes, no ANSYS needed (:doc:`testing`). ``--validate=long`` also runs the 3D ITER cases (about
   1 h more), ``--validate=no`` skips the validation.

All the cores are used, and everything is built in and for the active conda environment.

For developers
^^^^^^^^^^^^^^

* **Install** with ``./build.sh --pip=editable``.
* **After a change** of the C++ sources (``solver/``) or of the element library
  (``element_library/``, the package ``pytec3d-element-library``, compiled into the module), run ``./build.sh --pip=editable`` again: libMesh is
  reused and only the solver is rebuilt and reinstalled (``--validate=no`` to skip the validation while iterating).
  Changes of the Python sources need nothing with an editable install.
* **The C++ interface** of the solver (the model input, the session) is described in :doc:`solver_interface`.
* **Build type**: ``pyTEC3D_LIBMESH_METHOD=opt|devel|dbg ./build.sh``. ``opt`` (default, ``-O3``), ``devel``
  (``-O2 -g``, assertions on), ``dbg`` (``-O0 -g3``, all checks, compiler warnings are errors). Each method has its
  own libMesh installation, built the first time it is used.
* **Rebuild libMesh** from scratch (new clone, patches, configure, build): ``FORCE_LIBMESH=1 ./build.sh``.
* **The conda packages** locally (Linux): ``rattler-build build -r conda/recipe.yaml -c conda-forge`` (in an
  environment with ``rattler-build``); the packages are in ``output/``.

The Python package alone
------------------------

``pytec3d`` installs without a compiler, MPI or the solver, in any Python 3.12 environment:

.. code-block:: bash

   pip install .                                        # from the repository
   pip install "git+https://github.com/AlbertoBittes/pyTEC3D.git"

Everything works except loading a model into a solver (:meth:`pytec3d.Solver.load_model_into_solver`), which needs
the ``pytec3d-solver`` distribution of a ``build.sh`` run, and the check of the model data against the registries of
the solver (:func:`pytec3d.validation.check_model`).

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
   The solver (distribution ``pytec3d-solver``) is not built or not installed in the active environment: run
   ``./build.sh`` with the ``pytec3d`` environment active.

``ImportError: ... has interface version N, this pytec3d needs M``
   The installed solver was built from other sources than the installed ``pytec3d``: run ``./build.sh`` again.

``unresolved shared libraries (RPATH problem)`` at the end of the module build
   The module was linked against libraries outside the environment. Rebuild from a clean shell with only the conda
   environment active (no ``LD_LIBRARY_PATH`` pointing to system MPI or PETSc).

The module imports but crashes in the first parallel call
   A second MPI is installed by pip in the environment: ``pip uninstall mpi4py`` and reinstall it with conda.

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
