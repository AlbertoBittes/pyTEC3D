Tutorials
=========

Jupyter notebooks that build and run complete analyses step by step, taken from the test cases of ``tests/``, with
explanations, results and plots (matplotlib and `pyvista <https://docs.pyvista.org>`_). The pages below show them with
their stored outputs; to run them yourself:

.. code-block:: bash

   conda activate pytec3d
   python -m pip install -e ".[notebooks]"          # matplotlib, ipykernel (pyvista comes with the conda environment)
   conda install -c conda-forge jupyterlab
   cd docs/tutorials && jupyter lab

Each notebook copies its test case into a temporary folder, so the files of ``tests/`` are not modified.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Notebook
     - What it shows
   * - :doc:`tutorial_coil_2D`
     - The first analysis: model, model data, element library, boundary conditions, a DOF coupling (terminal of a
       circuit-coupled conductor), current- and voltage-driven steps, results, files, and the caveats of element
       results.
   * - :doc:`tutorial_stranded_coil_iron`
     - A stranded coil with its circuit equation (VOLT, EMF), a voltage step response, a saturating iron frame (B-H
       curve, Newton iterations, operating range of the law), with and without iron.
   * - :doc:`tutorial_coil_3d_sector`
     - A 3D coil on a 30 degree cyclic sector: cylindrical selections, element coordinate system, periodic couplings
       with rotated nodal frames, eddy currents of a passive ring, a full-ring 3D view.
   * - :doc:`tutorial_edge_elements`
     - Nedelec edge elements: midside-node Dirichlet values, the tree-cotree gauge and how to look at the tree, a
       shorted solid rod reacting to a source.
   * - :doc:`tutorial_loops_3d`
     - A quarter model with two symmetry planes: named boundary arrays, symmetry conditions, the current induced in a
       passive loop when the source is switched off, a mirrored 3D view.

.. rubric:: Plasma simulations

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Notebook
     - What it shows
   * - :doc:`tutorial_iter_equilibrium`
     - The ITER 15 MA free-boundary equilibrium (Grad-Shafranov, Picard iterations): the machine model, coil currents
       and initial guess, where the position-control parameters come from, the iterations one by one (convergence,
       control currents, plasma boundary), flux map and separatrix, current density, pressure and toroidal-field
       profiles, plasma volume, stored energy and poloidal beta, a 3D view. About 10 s.
   * - :doc:`tutorial_iter_evolution`
     - A vertical displacement event of the same plasma, coupled with the eddy currents of the vessel: the
       equilibrium as initial state, the scenario waveforms, the time loop with a Picard iteration per step, the
       trajectory of the plasma, the vessel currents, the transition from diverted to limited, snapshots and a 3D
       view. About 2 minutes.
   * - :doc:`tutorial_iter_3d_halo`
     - The same event on a 3D sector of the machine (cyclic symmetry), with the coils as voltage-driven stranded
       coils and halo currents: the 3D setup (element coordinate systems, coil circuits, periodic couplings with
       rotated frames, potential of the structures), the plasma projected on the 3D stations, core and halo
       regions, the halo current injected into the wall and its poloidal path in the vessel, coil current changes,
       vertical forces on the structures, a 3D view. About 20 minutes.

.. toctree::
   :maxdepth: 1
   :hidden:

   tutorial_coil_2D
   tutorial_stranded_coil_iron
   tutorial_coil_3d_sector
   tutorial_edge_elements
   tutorial_loops_3d
   tutorial_iter_equilibrium
   tutorial_iter_evolution
   tutorial_iter_3d_halo

Regenerating the stored outputs
-------------------------------

The notebooks are written, executed and saved by the scripts of ``docs/tutorials/tools`` (no Jupyter needed: pyvista
plots are rendered off screen):

.. code-block:: bash

   cd docs/tutorials/tools
   python build_tutorial_coil_2D.py ../tutorial_coil_2D.ipynb
   python build_tutorials.py ..                 # the other notebooks (or name one: edge_elements, loops_3d...)
   python build_plasma_tutorials.py ..          # iter_equilibrium, iter_evolution, iter_3d_halo (25 min)
