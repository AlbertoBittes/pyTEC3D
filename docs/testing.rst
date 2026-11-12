Testing and Benchmarks
======================

An installation of pyTEC3D is **validated** against the benchmark references with
``python pyTEC3D_APDL_benchmarks/bench.py validate`` (:doc:`installation`): the installed code runs the benchmark cases
and is compared with stored ANSYS Mechanical APDL results, about 10 minutes, without ANSYS. The validation fails when a
case is out of tolerance, and prints one summary table per suite (three suites):

* the **field formulations** against MAPDL on 23 analyses: 2D and 3D, nodal and edge elements on hexahedra, prisms
  and pyramids, stranded coils, massive and circuit-coupled conductors, time-integrated formulations, iron with a B-H
  curve;
* the **scalar-potential formulations** (``PHI3D_nodal``, ``TPHI3D_edge``, ``TPHI2D_edge``) against analytical
  solutions, the A, V-A types and MAPDL: magnetostatics with B-H curves and magnets, eddy currents with cuts and
  fixed-T coils, superconductor power laws (Bean's and Norris's AC losses, Jc(B));
* the **plasma algorithm** (:mod:`pytec3d.algorithm`) against `APEC <https://github.com/Fusion4Energy/APEC>`_, the
  APDL implementation of the same algorithm [APEC2018]_ [APEC2019]_ [APEC2020]_, on ITER equilibria and vertical
  displacements in 2D and 3D, with halo currents and stranded coils.

The same validation runs again with ``pytest tests/test_benchmarks.py -s``.

Everything else about the benchmarks (the cases and their inputs, the MAPDL / APEC models, the stored references, the
tolerances, the options of the validation, the parametric studies, the assumptions and the open points) lives in the
benchmark repository `pyTEC3D_APDL_benchmarks <https://github.com/AlbertoBittes/pyTEC3D_APDL_benchmarks>`_, the git
submodule ``pyTEC3D_APDL_benchmarks`` of pyTEC3D. Start from its `README
<https://github.com/AlbertoBittes/pyTEC3D_APDL_benchmarks/blob/main/README.md>`_; the comparison choices are in
`ASSUMPTIONS.md <https://github.com/AlbertoBittes/pyTEC3D_APDL_benchmarks/blob/main/ASSUMPTIONS.md>`_ and what is
still open in `OPEN_POINTS.md <https://github.com/AlbertoBittes/pyTEC3D_APDL_benchmarks/blob/main/OPEN_POINTS.md>`_.
