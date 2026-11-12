# Element library refactor: open points

State after the JSON element library refactor (branch `operator-version`, commits `df88eae`, `a84fca2`).
Reference for the physics: `main` at `3e16de0`; the comparison tests (no CREATE / APEC) match it.

## Decisions taken

- **Element field average:** JxW-weighted mean (sum JxW v / sum JxW, without the 2 pi r and thickness factor),
  so that sum(Jt * area) is the current through a section (`ElementEvaluator::mean_weights`).
- **Stranded coil volume:** as ANSYS (A-VOLT-EMF): Vc is the real coil volume (`volume`, or from `radius`), the
  modelled part is accounted for by `symmetry`; J = N i / S exactly only when symmetry * V_model = Vc.

## Deferred on purpose

- **Nonlinear nodal loads (`mn.`).** Linear for now; the generator rejects a law on a nodal parameter.
- **Law inputs other than a field at the quadrature points** (element means, circuit / global quantities).
- **Cross-coupling Jacobian blocks.** The Jacobian of a law term holds the derivative of the law's output with
  respect to its own input only. When the output also depends on other unknowns (an extra input of the law,
  `power_law_b`: rho(J, B), B from T and PHI; or a term whose trial is not the law input, sigma(T) on `dt(A)`), the
  blocks d(rho)/dB . J (of the T rows on the PHI and T columns through B) are missing: J uses the tangent at fixed B.
  A solution-dependent F term gets a residual but no Jacobian. Newton converges more slowly (linearly where the
  missing blocks matter: `kim_slab` 15-18 iterations per step against 8 without B), not to a wrong answer. Adding them
  needs the derivative of each output with respect to each extra input (generated from a "tangent_<input>" entry)
  and the kernels of the extra input's trials in the rows of the replaced coefficient.
- **History-dependent laws (hysteresis).** Laws are pointwise and stateless.
- **State variables of the laws.** Where the equation holds the time derivative of a quantity a law gives from the
  state (T-phi with a B-H law: dt(B(H)), B = mu(H) H), the residual of the term is the rate of the law's output by the
  time scheme (`Stage::RateResidual`, chosen by the generator when a kernel's trial is dt() of the law's own state
  input), not S(u) times the rate of u. The rates of the output are stored per element (any theta, Newmark): the
  storage is per element and per DOF of the output's block, not general state variables. Hysteresis (Preisach,
  Jiles-Atherton, play models) would need state variables per quadrature point, kept across steps and restarts.
- **Implicit laws.** A law's expressions are closed-form C++, evaluated in order: no local iteration. A law given as
  the inverse of another (sigma(E) of rho(J) with a floor, |E| = rho(J) |J|) cannot be written exactly; `power_law_e`
  is the closed-form inverse, exact only without floor. A generic implicit local in nonlinear_laws.json (residual,
  derivative, start, compiled into a safeguarded Newton loop) would allow it without hard-coding a law.

## Scalar potentials and superconductors: gaps (2026-11-12)

Measured on the `formulations` suite of pyTEC3D_APDL_benchmarks; user decisions in parentheses.

- **Circuits of the T-phi types.** No circuit unknowns: a conductor's net current is its cut DOF (current imposed by
  `_T`, voltage by `VLTG`, or free). Conductors in series share a coupling id; conductors in parallel at the same
  voltage each have their own DOF loaded by the same `VLTG`. Not available: a total current shared by conductors in
  parallel, a resistor or a supply impedance in series, an external network (a circuit of cut DOFs, as the coils of
  A, V-A have). Superconducting coils are usually the conductors themselves (cut + current or voltage drive), so this
  matters only for protection circuits, supplies with an impedance and current sharing between parallel tapes.
- **Cuts by hand.** The layer, the planar and crossing edges, their orientation and the couplings are made by the
  model (procedure: user manual, "Cuts"). No cut generator, no cohomology basis. One free cut DOF per independent
  hole; an extra free cut on the same loop is singular.
- **T0 of a coil given as fixed T by hand** (user: kept so). No helper building T0 from a winding on a general mesh.
- **PHI3D_nodal is static** (user: on purpose); its Hs adds an unknown per edge of the mesh (user: accepted; the
  workaround is Hs given on the edges).
- **Laws.** Up to MAX_LAW_INPUTS extra inputs (values only) and MAX_LAW_CURVES extra named curves. `power_law_b`:
  E(J, B) = E0(J / s(|B|)), |B| = mu0 |H| (non-magnetic conductor). Not available: the cross-coupling Jacobian (above);
  a curve of two arguments (Jc(B_par, B_perp), angular anisotropy of REBCO), n(B), temperature (no thermal model,
  no quench). Tabulated curves only (user: a fine table is the answer); the curve extends along its last segment, so
  E-J tables should end with a normal-state segment (`power_law_table(..., e_max=, rho_n=)`).
- **Hysteresis** of iron and magnets (user: deferred).
- **Nonlinear convergence at high n** (n = 30-100). A relative residual tolerance is relative to the first residual of
  the step, huge when the starting guess has J > Jc: it accepted E = 1e20 Ec at n = 100. An absolute tolerance must
  be scaled to the problem by the user (about Ec curl(w) V per T row): no automatic scale. A predictor (the previous
  rate extrapolated) as the first iterate would reduce the first residual; not implemented. Picard does not converge
  on the power law. `require_residual_reduction` is needed when a load changes sharply at a few cells (a transport
  current through a cut: `wire_2D_tphi` fails at its third step without it), but it stalls a step that starts from a
  current-free conductor in a changing field (`bulk_2D_tphi_axisym` after its field-cooled static step: the tangent at
  J = 0 is the floor, the first full step overshoots 1e4-1e5 times, the line search cannot reduce it): no single
  setting fits both, the user chooses (user manual, "Superconductors: the nonlinear settings"). A fallback (retry the
  step with the other choice, or a predictor) is not implemented. Adaptive time stepping (user: later).
- **Static steps and persistent currents.** A T-phi static step has T = 0 in the conductors (the model sets it):
  a static step after a transient wipes the induced currents; trapped fields and persistent currents need transient
  steps only.
- **2D.** `TPHI2D_edge` (planar, axisymmetric) is the 2D H-formulation. Not implemented: the T-A formulation of thin
  tapes (T on 1D lines, the tapes, in a 2D A mesh; 2D shells in 3D). It needs element types of a lower dimension
  than the mesh (EDGE3 cells inside a QUAD8 mesh: the mesh reader accepts mixed dimensions, the element types and the
  out-of-plane weight do not yet), the tangential derivative of A along the line (a vector coefficient, the tape
  direction) and a law on the sheet current.
- **Multidimensional curves** (proposed, not implemented). A law curve is a 1D table (piecewise linear, extended
  along its end segments); `power_law_b` combines two 1D curves, E0(J / s(|B|)), a separable form. A table on a
  rectangular grid of N axes (multilinear interpolation, the value and the slope along each axis; the 1D curve its
  one-axis case, existing laws unchanged) would give non-separable laws: Jc(B, theta) of REBCO (anisotropy), Jc(B, T),
  n(B). Keep the steep E0(J) as its own 1D table (fine log spacing in J); make the scaling s(B, theta, ...) the N-D
  table. The slopes along the axes are also what the cross-coupling Jacobian needs. Costs: the grid is rectangular,
  the data grows as the product of the axis sizes, each axis needs its own extension rule.
- **PHI2D_nodal** (in `ffb6b1dd`, not in the redesigned library). 2D coils and iron are the A_z types (`AVA2D_nodal`:
  one scalar unknown, no gauge, no cancellation); current-free 2D scalar-potential problems are `TPHI2D_edge` without
  conductors (`phi2d` benchmarks). A 2D twin of `PHI3D_nodal` (in-plane Hs, curl Hs = Js_z) would add nothing over
  A_z; to be restored only on request.
- **Edge values by the user.** Cuts, fixed-T coils and given Hs need the edge values in the solver's orientation:
  `tree_gauge.edges()` / `edge_circulations()` give them, and `set_value_on_selection` takes one value per selected
  id. No helper builds a cut layer or classifies its edges (the examples do it by position).
- **No restart from a saved state**: a new session starts from zero fields (and the stored rates of the law outputs,
  RateResidual, start at zero).
- **Integrals from cell means** (Joule power, AC loss as sum E . J V) are slightly low where the fields vary inside a
  cell (2 % on the ring of tutorial_tphi_cuts, 4 cells across); no quadrature-point integration of products is offered.
- **A, V-A with power_law_e** does not converge or oscillates at real superconductor parameters (slab benchmark):
  superconductors are a T-phi model.
- **PHI3D_nodal**: Hs carries the solenoidal part of Js, a non-divergence-free Js is used without a warning.
- **Setup errors not detected**: a missing or broken cut (the conductor is modelled cut open), wrong static-step
  gauges, two free cuts on one hole (singular).
- **Validation coverage**: no comparison with an independent H-formulation code on the published HTS benchmarks
  (HTS modelling workgroup); MAPDL references only for the reduced potential and A-V eddy currents.
- **Accuracy.** First order (no second-order Nedelec): the error decreases linearly with h (normal finite-element
  behaviour; refine). Power-law AC loss tends to the critical state as n grows (1.22 at n = 30, 1.09-1.11 at n = 100
  times Bean / Norris).

## Not verified

- **CREATE_* tests and the full APEC algorithm tests** were not run (too long). Their `model_data.json` were edited
  (`Az` -> `A`, laws in `nonlinear_laws`, 3D conductors `CIRCUIT_COUPLED_CONDUCTOR` -> `SOLID_CONDUCTOR`). Reduced
  runs of `ITER_2D`, `ITER_2D_evol` and `ITER_3D_evol_coils` (30 ms steps) are done by the plasma tutorials.
- **Second-order (Newmark) time scheme**: runs, but `main` has no reference to compare with.

## Known limits

- One FE order for all the Lagrange FEs of an element type (no mixed orders such as P2 / P1).
- No second-order Nedelec element: the edge fields are first order whatever option2. TPHI3D_edge and PHI3D_nodal are
  first order only ("second_order": false, checked by the Python validation: the Whitney pair of PHI and T, of PHI
  and Hs); AVA3D_edge* may run second-order V, VOLT, EMF (as MAPDL SOLID236).
- TPHI3D_edge has no source current and no circuit: a coil is T fixed on its edges (current drive only; the induced
  voltage from the flux afterwards). Voltage-driven coils are an A, V-A model.
- PHI3D_nodal is static: Hs adds one unknown per edge of the whole mesh, and with a B-H law Newton refactors the
  Hs block at each iteration (block-triangular, unsymmetric matrix).
- Coefficients constant over the element, or per quadrature point only through a law.
- 2D vector fields need a mesh whose dimension is that of the elements (no mixed 2D / 3D meshes).
- Coefficients of K, C, M are assumed constant in time between matrix rebuilds.

## Possible optimizations

- Two nodal fields with identical kernels each solve their columns (no cache per field term list).
- Kernels differing only by sign or a constant factor (`coil_av`, `-coil_av`) are computed twice.
- Transposed kernels (`sigma grad(V) . A` and `sigma A . grad(V)`) are computed separately.
- Algebraic relations between fields (H = nu B on linear elements) are not detected.
- The element matrix phase is about 0.4 s slower than `main` on the edge test (M^-1 K column solves done once per
  rebuild instead of a solve per output); the output phases are 3-4x faster, and the total is lower.

## Found while documenting the tests (2026-10-17)

- **Output rates for theta < 1** (resolved 2026-10-20). Derived fields used the mean rate `(u - u_n)/dt`, the
  equations the theta velocity: `coil_2D_avolt` with theta = 0.5 gave element `Jt` sums of +-2 MA for a 100 A load.
  They now use the scheme velocity (as MAPDL): 100.0 A. Theta = 1 is documented for sudden changes.

## Found while writing the plasma tutorials (2026-10-18)

- **Position control is proportional.** `Algorithm._apply_control_currents` sets `Js = base + I_corr / a` with the
  base fixed at the first call, and `I_corr` is recomputed from the residual field at the target at every iteration:
  the correction replaces the previous one instead of accumulating. The converged field at the target is not zero
  (`ITER_2D`: (1.6, 5.2) mT, axis 5 cm above the target, correction currents settled at +42 / -48 kA). An integral
  scheme (`base += I_corr / a`) would put the axis on the target.
- **`b*_value` of `ITER_2D` vs gains.** The values are the fields of unit gains (2x the response of the coils with
  the gains 0.5 used): the gains act as a 0.5 under-relaxation.
- **Halo fraction.** With `halo_delta1 = 0.04` the halo carries 20-40 kA of 15 MA (profile ~0 near psi_2): the halo
  force J_pol x B_phi is negligible against the eddy force (~30 MN at 0.63 s).
