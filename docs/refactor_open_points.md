# Element library refactor: open points

State after the JSON element library refactor (branch `operator-version`, commits `df88eae`, `a84fca2`).
Reference for the physics: `main` at `3e16de0`; the comparison tests (no CREATE / APEC) match it.

## Decisions pending

- **Element field average.** Element fields are the plain mean of the quadrature-point values, as on `main`.
  Candidates: area-weighted mean (sum JxW v / sum JxW, without the 2 pi r and thickness factor), so that
  sum(Jt * area) gives the current through a section; or volume-weighted. The difference is large where a field
  varies strongly inside an element (2D solid conductors on 8-node quads: up to 30-40 % per cell). Integral
  quantities (currents, losses, forces) are better computed by a dedicated integration than from element values.

## Deferred on purpose

- **Nonlinear nodal loads (`mn.`).** Linear for now; the generator rejects a law on a nodal parameter.
- **Law inputs other than a field at the quadrature points** (element means, circuit / global quantities).
- **Cross-coupling Jacobian blocks.** When a term's trial is not the law input (e.g. sigma(T) on `dt(A)`), J uses
  the secant only; a solution-dependent F term gets a residual but no Jacobian. Newton converges more slowly,
  not to a wrong answer.
- **History-dependent laws (hysteresis).** Laws are pointwise and stateless.

## Not verified

- **CREATE_* tests and the full APEC algorithm tests** were not run (too long). Their `model_data.json` were edited
  (`Az` -> `A`, laws in `nonlinear_laws`, 3D conductors `CIRCUIT_COUPLED_CONDUCTOR` -> `SOLID_CONDUCTOR`). Reduced
  runs of `ITER_2D`, `ITER_2D_evol` and `ITER_3D_evol_coils` (30 ms steps) are done by the plasma tutorials.
- **Second-order (Newmark) time scheme**: runs, but `main` has no reference to compare with.

## Known limits

- One FE order for all the Lagrange FEs of an element type (no mixed orders such as P2 / P1).
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

- **Output rates for theta < 1.** Derived fields use the mean rate `(u - u_n)/dt`, the equations the theta velocity.
  After a sudden change they differ a lot: `coil_2D_avolt` with theta = 0.5 gives element `Jt` sums of +-2 MA for a
  100 A load, alternating in sign (theta = 1: 155 A, the plain-mean issue above).

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
