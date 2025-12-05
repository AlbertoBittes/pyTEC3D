Installation
============

Requirements
------------
- Python 3.9 or newer
- C++ toolchain for building/running the solver (e.g., GCC/Clang and CMake)
- Recommended: a virtual environment for the Python tools

Install the Python API
----------------------

From a cloned repository:

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate  # on Windows use .venv\\Scripts\\activate
   pip install --upgrade pip
   pip install -e .

This installs the command-line entry points defined in ``pyproject.toml`` (e.g., ``apec`` and ``convert_mesh``).

Install the C++ solver
----------------------

The solver sources live under ``src/`` with accompanying headers in ``include/``.

Basic build outline (fill in with your project-specific steps):

.. code-block:: bash

   mkdir -p build
   cd build
   cmake ..
   cmake --build .

If you have specific compiler flags, third-party libraries, or environment modules, document them here. Replace this section with the authoritative build/run recipe your team uses.

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.
- ``ansys-mapdl-reader`` enables ANSYS interoperability.

When building on Read the Docs, these heavy dependencies are mocked so documentation builds quickly. They are still required for full runtime functionality.
