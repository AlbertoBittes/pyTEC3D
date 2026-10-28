Theory
======

This part derives the equations that pyTEC3D solves, from Maxwell's equations to the discrete system, and explains
how the axisymmetric plasma equilibrium is coupled with the eddy currents of the structures. It follows the element
library: every term of an element type in :doc:`../user_manual/element_library` (and in
``element_library/pytec3d_element_library/data/element_types.json``) is derived here, and every section points to the element type and the
domain that implement it.

**Foundations.**

* :doc:`domains`: the seven domains of pyTEC3D (``NON_CONDUCTING`` ... ``STRANDED_COIL``), their electric and magnetic
  character and the unknowns each one has in 2D and 3D.
* :doc:`notation`: the symbols, with the names of the corresponding fields, parameters and arrays of the code.
* :doc:`maxwell`: the magneto-quasi-static (eddy-current) approximation, the A-V potentials, the gauge and the
  boundary conditions.
* :doc:`constitutive`: reluctivity, conductivity, magnetization and the B-H law, in the element coordinate system of
  a cell.

**Formulations.**

* :doc:`formulation_2d`: ``AV2D_nodal``, planar and axisymmetric: the out-of-plane vector potential, the out-of-plane
  voltage of a conductor with a terminal, the element matrices of each domain.
* :doc:`formulation_3d`: ``AV3D_nodal`` (Lagrange vector potential with a Coulomb gauge term) and ``AV3D_edge``
  (Nedelec elements with a tree-cotree gauge).
* :doc:`circuits`: the domains driven through terminals: stranded coils and their circuit equation, conductors with
  a terminal in 2D, injected currents in 3D, and the time-integrated variants (``_timint``).
* :doc:`nonlinear`: saturating materials: residual, Jacobian and Newton-Raphson iteration in block form.

**The code.** :doc:`implementation`: how the kernels of the element types become the element matrices, the time
scheme, the constraints (Dirichlet values, nodal frames, couplings) and the output fields.

**Plasma.** :doc:`plasma`: the Grad-Shafranov equilibrium, the free-boundary problem and its Picard solution with
position control, the plasma evolution coupled with the eddy currents of 2D and 3D structures, the halo currents.

.. toctree::
   :maxdepth: 2

   domains
   notation
   maxwell
   constitutive
   formulation_2d
   formulation_3d
   circuits
   nonlinear
   implementation
   plasma
   references
