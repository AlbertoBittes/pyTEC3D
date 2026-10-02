Getting Started
===============

Why --code name--?
-------------
Plasma disruptions can induce high electromagnetic loads on tokamak components,
making accurate plasma modelling crucial for design. pyTEC3D ([ref] + link to repo) is a FEM code developed at Fusion for Energy,
based on ANSYS APDL, to simulate plasma equilibrium and evolution.
The main difference with other codes, which solve completely axisymmetric problems,
is that pyTEC3D allows direct coupling of the axisymmetric plasma equilibrium equations
with the full 3D eddy currents equations in the passive structures.
While leveraging ANSYS features, pyTEC3D's dependence on ANSYS introduces performance issues,
including iteration restarts and repeated matrix factorization, limiting scalability for large models.
Additionally, reliance on costly HPC licenses restricts its usability with large models.
--code name-- aims to address these challenges as an independent, free, scalable, and customizable code
able to replicate pyTEC3D's algorithm and overcome its limitations.


What is --code name--?
-------------

--code name-- is a FEM code for time-domain EMAG eddy current problems, employing the A, V-A formulation [ref].
The code is separated into 2 distinct modules:
- a Pre/Post-Processing module, written in Python, that handles model and analysis preparation and output visualization.
- a Solver, written in C++, that performs the actual FEM calculations.

Pre/Post-Processor module includes a Graphical User Interface and/or Python API for scripting:
- ParaView (https://www.paraview.org/) to visualize and manipulate the mesh. 
- Python API based on VTK (https://vtk.org) libraries for automatizing the Pre/Post-Processing and interact with the GUI in ParaView.

--code name-- models are based on VTK format. The user can generate a mesh with any meshing tool and import it
as a VTK Unstructured Grid (link to vtk data formats page). By doing so --code name-- leverages the consolidated and powerful
capabilities of Python as a programming language and Paraview as a visualization tool to
prepare and post-process the models and setup the analyses.

---example of GUI screenshot---

The solver module is a C++ code that:
- Reads the model prepared with the Pre/Post-Processing module.
- Distributes FEM Degrees-Of-Freedom and build system matrix in parallel, by leveraging the open-source FEM library libMesh (https://libmesh.github.io/), according to A-V, A formulation of Maxwell equations, for 2D, 2D axisymmetric and 3D meshes.
- Solves in parallel the FEM system by using open-source PETSC/MUMPS solvers (https://petsc.org/release/, https://mumps-solver.org/index.php)

Project links
-------------
- Source repository: https://github.com/AlbertoBittes/PlasmaArc
- Issue tracker: use the GitHub issues tab on the repository
- Contacts: <maintainer@example.com> (replace with your team contacts)

What you get
------------
- Geometry and mesh preparation helpers for pyTEC3D/PlasmaArc workflows
- IO bridges for VTK, Excel, ANSYS, and ParaView integration
- A C++ solver configured via input files generated from the Python API

Quick start
-----------
1. Install the Python API and its dependencies (see :doc:`installation`).
2. Prepare or import a mesh using the Python helpers.
3. Configure solver inputs, run the C++ solver, and visualize results.
4. Browse :doc:`examples` for end-to-end walkthroughs once they are added.
