User Manual
===========

This manual explains how pyTEC3D works and how to set up an analysis, one step at a time. Read
:doc:`workflow` first: it describes the objects of the package and the life cycle of an analysis. The
:doc:`tutorial notebooks <../tutorials/index>` then run complete analyses interactively. The other pages are a
reference for each step; the :doc:`element_library` documents every formulation, its fields, materials, real
constants, loads and outputs.

All quantities are in SI units (m, kg, s, A, V, T, Ohm m). The mesh coordinates must be in metres.

.. toctree::
   :maxdepth: 2
   :caption: Using pyTEC3D

   workflow
   mesh
   model_data
   element_library
   selection
   boundary_conditions
   loads
   couplings
   coordinate_systems
   nonlinear
   analysis
   solver_settings
   outputs
   plasma_algorithm
   paraview
   validation
