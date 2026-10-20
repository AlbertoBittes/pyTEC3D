pyTEC3D Documentation
=====================

pyTEC3D is a finite-element code for low-frequency (eddy-current) electromagnetic analyses, written for the
simulation of plasma equilibria, plasma evolutions and disruptions and of their effects on the conducting
structures of a tokamak. A Python package prepares the model, drives the analysis and post-processes the results; a
parallel C++ solver (libMesh, PETSc, MUMPS) assembles and solves the finite-element systems.

* New users: start from :doc:`getting_started`, install the code with :doc:`installation`, then follow the
  :doc:`tutorials/index` (Jupyter notebooks built on the test cases).
* Setting up an analysis: the :doc:`user_manual/index` describes every step and option, and the
  :doc:`user_manual/element_library` every element type, its fields, materials, loads and outputs.
* The equations behind the element types are in :doc:`theory/index`; complete annotated analyses are in
  :doc:`examples/index`; the Python reference is in :doc:`python_api`.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   installation
   testing
   tutorials/index
   user_manual/index
   theory/index
   examples/index
   python_api
