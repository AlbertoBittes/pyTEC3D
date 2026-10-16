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

- **CREATE_* and APEC algorithm tests** were not run (too long). Their `model_data.json` were edited
  (`Az` -> `A`, `option4` -> `"law": "ferromagnetic"` on the iron material).
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
