Installation
============

pyTEC3D has three parts, independent distributions built from the same repository:

* the **Python API** ``pytec3d`` (``python/``, ``python/pyproject.toml``): model preparation, analysis driver, plasma
  algorithm, post-processing helpers. Pure Python: it installs anywhere, with ``numpy`` and ``vtk``.
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
       conda-forge and ``pytec3d-libmesh``. Windows: built in MSYS2 (``ci/windows/build.sh``), PETSc and libMesh
       linked statically into the module, which is shipped with the DLLs of MUMPS, ScaLAPACK, OpenBLAS and the
       Fortran runtime of MSYS2; MS-MPI itself is the conda package ``msmpi``.
   * - ``pytec3d-libmesh``
     - linux-64
     - The patched libMesh the solver links: the libMesh development branch at the time of the build
       (``solver/libmesh/build_libmesh.sh``; the commit is in ``share/pytec3d/libmesh-commit``).

**MPI.** The solver depends on the MPI it is linked with (Open MPI on Linux, MS-MPI on Windows), and ``mpi4py``, a
dependency of ``pytec3d-base``, has a conda build per MPI: in an environment with the solver, conda always installs the
mpi4py of that MPI (a conda-forge environment holds only one MPI). ``conda install pytec3d-base`` alone installs the
mpi4py of any MPI; ``conda install pytec3d-solver`` installs the solver alone. In parallel:
``mpiexec -n 4 python macro.py``.

**The Python API with pip.** ``pip install ./element_library ./python`` (from the root of the repository; the element
library first, the API depends on it) installs ``pytec3d`` alone. ``pip install "./python[mpi]"`` adds mpi4py,
which then builds against the MPI it finds (``mpicc`` on ``PATH``), or uses an MPI installed with pip next to it, e.g.
``pip install "./python[mpi]" openmpi`` (or ``mpich``, ``impi-rt``). With a solver of a conda environment, install mpi4py with
conda instead, so that it uses the MPI of the solver.

**Publishing.** The workflow builds the packages at every push and pull request (the artifacts ``conda-linux-64`` and
``conda-win-64``, a channel layout), installs them on both systems and solves a model in serial and on 2 MPI ranks
(``conda/tests/test_square_conductor.py``). On a tag ``v*`` it uploads them to the anaconda.org channel of the
repository variable ``ANACONDA_OWNER`` (with the secret ``ANACONDA_API_KEY``) and attaches them to the GitHub release.
To use the packages of a run without a channel: download the two artifacts into one directory, index it with
``rattler-index fs <dir>`` and install with ``-c file://<dir> -c conda-forge``.

**Python.** Every package is for Python 3.12, the Python of ParaView, which runs the API in its plugin and console
(:doc:`user_manual/paraview`).

**The .cdb converter.** ``pytec3d-convert-cdb`` (:func:`pytec3d.mesh_converter.convert_cdb`) needs
``ansys-mapdl-reader``, which is not a dependency: its declared requirements conflict with recent NumPy, VTK and
pyvista, while it works with them. Install it after pytec3d, without its dependencies:

.. code-block:: bash

   pip install --no-deps ansys-mapdl-reader appdirs tqdm

Development
-----------

pyTEC3D is only built as packages: the solver is never compiled into an environment by hand.

* **The environment**: ``environment.yml`` (conda-forge, Python 3.12) has the dependencies of the API, the tools of the
  repository (pytest, ruff, mypy, Sphinx, JupyterLab) and rattler-build:

  .. code-block:: bash

     git clone https://github.com/AlbertoBittes/pyTEC3D.git && cd pyTEC3D
     git submodule update --init pyTEC3D_APDL_benchmarks     # the benchmarks of the validation
     conda env create -f environment.yml
     conda activate pytec3d-dev

* **The Python API** from the sources, editable (a change of ``python/`` or ``element_library/`` takes effect at once):

  .. code-block:: bash

     pip install --no-deps -e ./element_library -e ./python

* **The packages**, solver included (``conda/recipe.yaml``; the first build compiles libMesh, 30 to 60 minutes):

  .. code-block:: bash

     rattler-build build -r conda/recipe.yaml -c conda-forge                 # into output/
     conda install -c ./output -c conda-forge pytec3d-solver                 # the solver into this environment

  ``LIBMESH_REF=<branch, tag or commit>`` builds another libMesh than ``master``. After a change of the solver
  (``solver/``) or of the element library, build the packages again: the element library is compiled into the
  solver. The C++ interface of the solver is described in :doc:`solver_interface`; its C++ example
  (``solver/examples``) is built with ``-DPYTEC3D_BUILD_CPP_EXAMPLES=ON`` when the solver is built with CMake.

* **The validation** against the benchmark references (:doc:`testing`), with the installed pyTEC3D:

  .. code-block:: bash

     python pyTEC3D_APDL_benchmarks/bench.py validate        # or: pytest tests/test_benchmarks.py -s

* **Checks**: ``ruff check .`` and ``mypy --config-file python/pyproject.toml`` (from the root of the repository).

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

   mpiexec -n 4 python macro.py

Every rank runs the script; the mesh is distributed by libMesh when the model is loaded, and only rank 0 writes the
output files and, by default, prints (``Solver(root_only_output=True)``). Threads can be added per MPI process with
``Solver(num_threads=n)``. See :doc:`user_manual/solver_settings` for the solver options of large models, and the
`Open MPI documentation of mpirun <https://docs.open-mpi.org/en/main/man-openmpi/man1/mpirun.1.html>`_ (Linux) or of
`MS-MPI <https://learn.microsoft.com/en-us/message-passing-interface/mpiexec>`_ (Windows) for process placement.

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
   The solver (package ``pytec3d-solver``) is not installed in the active environment:
   ``conda install -c <channel> -c conda-forge pytec3d-solver``, or ``pytec3d`` for everything. Without it, everything
   works except loading a model into a solver.

``ImportError: ... has interface version N, this pytec3d needs M``
   The installed solver and API are of different releases: install both of the same release (``conda install
   pytec3d`` keeps them together).

The module imports but crashes in the first parallel call
   A second MPI is installed by pip in the environment: ``pip uninstall mpi4py`` and reinstall it with conda.

A parallel run stops with ``UCX ERROR ... uct_iface_open`` or ``Failed to create UCP worker`` (Linux)
   Open MPI tried a network adapter it cannot use (e.g. on Azure virtual machines). On one machine, shared memory
   and TCP are enough: ``export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm,tcp OMPI_MCA_osc=^ucx`` before ``mpiexec``.

The element library changed but the solver rejects the new names
   The element library is compiled into the solver: build the packages again (``rattler-build build -r
   conda/recipe.yaml -c conda-forge``) and install the new solver.

The validation fails
   The summary printed by ``bench.py validate`` gives the cases and quantities out of tolerance; see :doc:`testing`.

More output, or less
   ``pytec3d.set_log_level("DEBUG")`` or the environment variable ``PYTEC3D_LOG_LEVEL=DEBUG``; ``NO_COLOR=1``
   disables the colored labels. PETSc options such as ``-ksp_monitor`` are set in ``solver_settings["petsc_options"]``
   (see :doc:`user_manual/solver_settings`).
