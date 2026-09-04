Installation
============

Requirements
------------
- Linux
- conda

Install from conda channel
--------------------------

Create a clean environment and install ``apec-m``:

.. code-block:: bash

   conda create -n apec-m -c <your-channel> -c conda-forge python=3.12 apec-m -y
   conda activate apec-m

Install from local build artifacts
----------------------------------

If you built the package locally with ``conda mambabuild``:

.. code-block:: bash

   conda create -n apec-m -c conda-forge python=3.12 -y
   conda activate apec-m
   conda install --use-local -y apec-m

Install partial stacks from conda
---------------------------------

Coordinate package only:

.. code-block:: bash

   conda install -c <your-channel> -c conda-forge apec-m-coordinate

Python API package (depends on ``apec-m-coordinate``):

.. code-block:: bash

   conda install -c <your-channel> -c conda-forge apec-m-python

Conda package build (developers)
--------------------------------

From repository root:

.. code-block:: bash

   conda create -n apec-m-build -c conda-forge python=3.12 conda-build boa compilers git -y
   conda activate apec-m-build
   export APEC_M_LIBMESH_METHOD=opt
   conda mambabuild conda-recipe/recipes/apec-m-coordinate -c conda-forge --output-folder conda-dist
   conda index conda-dist
   conda mambabuild conda-recipe/recipes/apec-m-python -c conda-forge -c file://$PWD/conda-dist --output-folder conda-dist
   conda index conda-dist
   conda mambabuild conda-recipe/recipes/apec-m-core -c conda-forge -c file://$PWD/conda-dist --output-folder conda-dist
   conda index conda-dist
   conda mambabuild conda-recipe/recipes/apec-m -c conda-forge -c file://$PWD/conda-dist --output-folder conda-dist

Optional dependencies
---------------------
- ``vtk`` is required for mesh IO and visualization helpers.
- ``pandas`` and ``numpy`` power tabular and numeric processing.

When building on Read the Docs, these heavy dependencies are mocked so documentation builds quickly. They are still required for full runtime functionality.
