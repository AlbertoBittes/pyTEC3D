References
==========

Finite elements for eddy currents
---------------------------------

.. [BiroPreis1989] O. Biro, K. Preis, "On the use of the magnetic vector potential in the finite element analysis of
   three-dimensional eddy currents", IEEE Transactions on Magnetics 25(4), 3145-3159 (1989).
.. [Biro1999] O. Biro, "Edge element formulations of eddy current problems", Computer Methods in Applied Mechanics and
   Engineering 169, 391-405 (1999).
.. [Nedelec1980] J.-C. Nedelec, "Mixed finite elements in R3", Numerische Mathematik 35, 315-341 (1980): the edge
   (curl-conforming) elements of ``AV3D_edge``.
.. [Bossavit1998] A. Bossavit, *Computational Electromagnetism: Variational Formulations, Complementarity, Edge
   Elements*, Academic Press (1998).
.. [AlbaneseRubinacci1988] R. Albanese, G. Rubinacci, "Integral formulation for 3D eddy-current computation using
   edge elements", IEE Proceedings A 135(7), 457-462 (1988): the tree-cotree gauge of edge elements.
.. [MangesCendes1995] J. B. Manges, Z. J. Cendes, "A generalized tree-cotree gauge for magnetic field computation",
   IEEE Transactions on Magnetics 31(3), 1342-1347 (1995).
.. [Meunier2008] G. Meunier (ed.), *The Finite Element Method for Electromagnetic Modeling*, ISTE / Wiley (2008):
   potential formulations, stranded and massive conductors, circuit coupling.
.. [MAPDLTheory2024] *Mechanical APDL Theory Reference*, Release 2024 R2, section 5.4 "Stranded Coil Analyses"
   (A-VOLT-EMF formulation: coil volume, symmetry factor, time-integrated variant).

Time integration and nonlinear solution
---------------------------------------

.. [Hughes2000] T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*,
   Dover (2000): the theta method and the Newmark family for first- and second-order systems.
.. [CrankNicolson1947] J. Crank, P. Nicolson, "A practical method for numerical evaluation of solutions of partial
   differential equations of the heat-conduction type", Proceedings of the Cambridge Philosophical Society 43(1),
   50-67 (1947).
.. [Newmark1959] N. M. Newmark, "A method of computation for structural dynamics", Journal of the Engineering
   Mechanics Division (ASCE) 85(EM3), 67-94 (1959).
.. [DennisSchnabel1996] J. E. Dennis, R. B. Schnabel, *Numerical Methods for Unconstrained Optimization and Nonlinear
   Equations*, SIAM, Classics in Applied Mathematics 16 (1996): Newton's method with line search.
.. [Brent1973] R. P. Brent, *Algorithms for Minimization without Derivatives*, Prentice-Hall (1973).

Linear solvers and software
---------------------------

.. [Kirk2006] B. S. Kirk, J. W. Peterson, R. H. Stogner, G. F. Carey, "libMesh: a C++ library for parallel adaptive
   mesh refinement/coarsening simulations", Engineering with Computers 22(3-4), 237-254 (2006);
   https://libmesh.github.io/.
.. [PETScManual] S. Balay et al., *PETSc/TAO Users Manual*, Argonne National Laboratory, ANL-21/39;
   https://petsc.org/release/manual/.
.. [Amestoy2001] P. R. Amestoy, I. S. Duff, J. Koster, J.-Y. L'Excellent, "A fully asynchronous multifrontal solver
   using distributed dynamic scheduling", SIAM Journal on Matrix Analysis and Applications 23(1), 15-41 (2001):
   MUMPS, https://mumps-solver.org/.
.. [KarypisKumar1998] G. Karypis, V. Kumar, "A fast and high quality multilevel scheme for partitioning irregular
   graphs", SIAM Journal on Scientific Computing 20(1), 359-392 (1998): METIS.
.. [SaadSchultz1986] Y. Saad, M. H. Schultz, "GMRES: a generalized minimal residual algorithm for solving
   nonsymmetric linear systems", SIAM Journal on Scientific and Statistical Computing 7(3), 856-869 (1986).
.. [HiptmairXu2007] R. Hiptmair, J. Xu, "Nodal auxiliary space preconditioning in H(curl) and H(div) spaces", SIAM
   Journal on Numerical Analysis 45(6), 2483-2509 (2007): the preconditioner of curl-curl systems (hypre AMS).
.. [Schroeder2006] W. Schroeder, K. Martin, B. Lorensen, *The Visualization Toolkit: An Object-Oriented Approach to
   3D Graphics*, 4th ed., Kitware (2006); https://vtk.org/.
.. [Ahrens2005] J. Ahrens, B. Geveci, C. Law, "ParaView: an end-user tool for large data visualization", in *The
   Visualization Handbook*, Elsevier (2005); https://www.paraview.org/.

Plasma equilibrium and evolution
--------------------------------

.. [APEC2018] P. Testoni, M. Cavinato, A. Portone, "A novel code for the simulation of plasma equilibrium and
   evolution", Fusion Engineering and Design (2019), https://doi.org/10.1016/j.fusengdes.2019.02.019 (SOFT 2018).
.. [APEC2019] P. Testoni, A. Portone, "Introduction and validation of the APEC code for the simulation of plasma
   equilibrium and evolution in presence of 3D passive structures", IEEE Transactions on Magnetics 56(4), 2020
   (COMPUMAG 2019).
.. [APEC2020] P. Testoni, A. Portone, "Upgrade of the APEC code: implementation of halo currents and multi-turn coils
   analyses", SOFT 2020.
.. [Shafranov1958] V. D. Shafranov, Sov. Phys. JETP 6, 545 (1958); Zh. Eksp. Teor. Fiz. 33, 710 (1957).
.. [GradRubin1958] H. Grad, H. Rubin, "Hydromagnetic equilibria and force-free fields", Proceedings of the 2nd United
   Nations International Conference on the Peaceful Uses of Atomic Energy, Geneva, vol. 31, 190 (1958).
.. [Freidberg2014] J. P. Freidberg, *Ideal MHD*, Cambridge University Press (2014): the equilibrium of toroidal
   plasmas, the Grad-Shafranov equation, poloidal beta and internal inductance.
.. [Jardin2010] S. C. Jardin, *Computational Methods in Plasma Physics*, CRC Press (2010): the numerical solution of
   the free-boundary equilibrium problem, Picard iterations and their vertical instability.
.. [Wesson2011] J. Wesson, *Tokamaks*, 4th ed., Oxford University Press (2011).
.. [Hender2007] T. C. Hender et al., "Progress in the ITER Physics Basis, Chapter 3: MHD stability, operational
   limits and disruptions", Nuclear Fusion 47(6), S128-S202 (2007): disruptions, vertical displacement events and
   halo currents.
.. [HumphreysKellman1999] D. A. Humphreys, A. G. Kellman, "Analytic modeling of axisymmetric disruption halo
   currents", Physics of Plasmas 6(7), 2742-2756 (1999).
.. [Villone2013] F. Villone, L. Barbato, S. Mastrostefano, S. Ventre, "Coupling of nonlinear axisymmetric plasma
   evolution with three-dimensional volumetric conductors", Plasma Physics and Controlled Fusion 55(9), 095008
   (2013).
