# Code review: the libMesh and MFEM solvers (November 2026)

The two solver branches as of 2026-11-15:

- `libmesh-global-assembly`: the libMesh solver of `main`, with global assembly. K, C and M are assembled once as
  global unconstrained matrices. The constraints become a prolongation P with values g (u = P w + g), and the
  reduced matrix P^T K P is factorized again with the numeric phase only when the pattern is unchanged. The
  element-matrix files (`emat/`) and the hash table are gone, the field operators are kept in memory, and the RHS of
  constant-coefficient elements is cached. The element loops and output writing are threaded (OpenMP).
- `mfem-backend`: the port to MFEM 4.9. It uses ParMesh, batched geometric factors and DofToQuad-like reference
  tables, ParFiniteElementSpace, HypreParMatrix RAP, and PETSc only for KSP / MUMPS.

Both pass the whole validation (simple, iter, formulations: 56 passed, 2 skipped, the same values as `main`). They
also pass the new patch tests (below) with identical results.

## 1. Performance

ITER_3D_evol_coils_edge, 3 steps, 68k HEX20 cells, 256k DOFs, one MPI rank:

| solver | 1 thread | 8 threads | peak RSS (1 thread) |
|---|---:|---:|---:|
| libMesh `main` (emat files, per-element constraints) | 124.6 s | – | 2.08 GB |
| libMesh `libmesh-global-assembly` | 89.9 s | 35.6 s | 2.47 GB |
| MFEM `mfem-backend` | 33.9 s | 23.0 s | 2.68 GB |

Phase times of the 3 steps (libMesh branch, 1 / 8 threads; MFEM, 1 thread):

| phase | libMesh 1 | libMesh 8 | MFEM 1 |
|---|---:|---:|---:|
| operators (K, C, M, field operators) | 4.1 s | 0.8 s | 2.9 s |
| linear system (RHS + solve) | 10.0 s | 6.1 s | 7.6 s |
| nodal fields | 8.2 s | 1.3 s | 4.9 s |
| element fields | 13.5 s | 2.0 s | 1.9 s |

Most of the remaining serial gap on libMesh is output: the nodal values of Nedelec fields come from
`FEInterface::nodal_soln`, one FE reinit per element and variable, plus the element fields. That cost is the
recomputation libMesh does by design (no stored quadrature data); threads remove most of it. MFEM evaluates the
output fields as stored sparse operators (one matvec per field) on batched geometry.

## 2. libMesh solver: findings and possible improvements

Ordered by importance. File references are to `libmesh-global-assembly`.

### Correctness and robustness

1. **Operators are never invalidated** (`pytec3d_system.cpp`, `write_element_matrices`,
   `check_update_element_matrices() == false`). K, C, M, the field operators, the RHS bases and the law-element list
   are built at the first solve and kept for the whole session. A change of MATERIAL / CONSTANT / TYPE arrays or of
   material data between steps (`pyTEC3DModel::set_arrays`, a new `model.data`) is silently ignored. `main` behaves
   the same way (its emat files), so this is not a regression, but it is a trap. Suggestion: an `operators_valid_`
   flag reset by `set_arrays` and by a model-data reload (or a fingerprint of the coefficient arrays), and an error
   or warning when a cached quantity changes.
2. **The RHS basis cache relies on the load coefficients being constant per element.** This is decided per element
   at assembly time (`build_rhs_basis`). A coefficient that becomes position- or time-dependent through a law must
   keep returning `false` there. Worth a unit test that a quadrature-point-dependent load takes the full path.
3. Leftovers: the unused alias `ElemRecordIndex` (`pytec3d_system.cpp`), the "Picard iterations" comment of the
   constraint cache (`constraints.cpp`), and `check_update_element_matrices()`, which always returns false (remove it,
   or implement item 1).

### Scalability (MPI)

4. **`ReducedPetscLinearSolver::set_constraints` allgathers every constrained DOF to every rank**, to number the
   unconstrained DOFs. Memory and communication grow as O(N_constrained) per rank at each constraint update. Each
   rank only needs the reduced index of the masters it references: an exclusive scan of the local free counts, plus
   a PetscSF / ghost exchange for the off-rank masters.
5. **`check_update_boundary_conditions` copies the whole previous constraint map** (a `std::map` of `std::map`) at
   every step to compare it. Compare in place after `swap_dof_constraints()`, or compare a hash of the rows.
6. **`build_time_integration_constraints` loops over all active elements** (not only the local and ghost ones) and
   fills `unordered_set`s of all their DOFs on every static step. Restrict it to `active_local_element_ptr_range()`
   plus the ghosts, and cache the result (it only changes with `transient_effects` or the element types).
7. `build_dof_to_variable_index_map` builds an `unordered_map` over all local DOFs at every constraint update, only
   for the statistics table. Count per node and variable instead, or use the DOF ranges of the variables.
8. **Python-side replication**: every rank holds the whole VTK grid and the whole field arrays (`ModelInput` is
   copied into `pyTEC3DModel::arrays_`). `write_system_to_arrays` allreduces full-length arrays
   (n_points x components) on every rank. This is consistent with the replicated grids, but memory per rank is O(N).
   Reducing to rank 0 only (if only rank 0 writes VTK) would save one full copy per rank.

### Performance

9. **Field operators as dense blocks per element** (`nodal_operators_`, `element_operators_`). This is about 0.4 GB
   above `main` on the large model. Assembling them once into one sparse matrix per field system (as the MFEM branch
   does) would cut memory, and each output would become one `MatMult` with no `dof_indices` per element.
10. **Nonlinear `update_time_derivatives`**: the law-rate loop (`law_output`, one FE reinit per element and step) is
    serial and keyed by element id in a hash map. Thread it with `parallel_elements` and store the rates in a vector
    indexed like `law_elements()`.
11. `assemble_system_impl`: `jacobian.zero(); jacobian.add(1.0, keff)` calls `MatAXPY(DIFFERENT_NONZERO_PATTERN)`.
    PETSc detects the identical pattern, so the factorization stays numeric-only (checked: "1 full, 8 numeric only"
    on coil_2D_avoltemf_ferr). `MatCopy(..., SAME_NONZERO_PATTERN)` is cheaper and states the intent.
12. `assemble_rhs` allocates two `zero_clone()` vectors per history term per step: keep them as members.
13. `compute_law_ranges` runs one `comm().max` per (law, id): pack them into one vector reduction.

### Design

14. The per-element FE reinit is libMesh's model, kept deliberately: the cached RHS bases and the threaded loops
    already remove most of its cost in linear runs. In nonlinear runs the law terms are still reinitialized at every
    Newton iteration, which is the main remaining gap to MFEM; MPI ranks and threads mitigate it.
15. Constraints through P (`reduced_linear_solver`) plus global matrices is a clean replacement of the per-element
    constraining, and it is the same algebra as the MFEM branch. This makes the two code paths easier to compare.

## 3. MFEM compared with libMesh

### Advantages of MFEM

- **Speed**: 2.7x faster serial than the libMesh branch, 1.5x with 8 threads; 3.7x faster than `main`. The gain
  comes from batched geometry and shape tables (no per-element reinit), forms assembled once and RAP'd, field
  operators as sparse matrices, and no constraint work per element.
- **Assembly model**: the standard MFEM model (local forms -> P^T A P -> HypreParMatrix). It leaves room for
  partial assembly and device kernels later.
- **Arbitrary order** is native to the MFEM spaces (H1 and Nedelec of any order on tets, hexes, prisms). The
  solver's addressing of high-order DOFs through point arrays still has to be generalized before option2 > 2 is
  exposed.
- Simpler parallel layer: one owner per true DOF, `GroupCommunicator` for the shared ones, no ghosting functors.

### Disadvantages and gaps of MFEM

- **Serendipity elements are in-house** (QUAD8, HEX20 and PRISM15 by Vandermonde; PYRAMID13 the libMesh rational
  basis). MFEM has no 3D serendipity, so pyTEC3D maintains them.
- **No mixed-geometry GeometricFactors in MFEM**: the solver keeps its own batched quadrature data per geometry.
- **Two linear-algebra stacks**: hypre for assembly and RAP, PETSc for KSP / MUMPS (converted at each full
  factorization). A PETSc-only path (MFEM's `PetscParMatrix`) is possible but MFEM's PETSc support is thinner.
- **Memory**: about 9 % above the libMesh branch on the large model (2.68 vs 2.47 GB): forms, the RAP'd matrices
  and the PETSc copy.
- **Ported pieces to maintain**: the libMesh NewtonSolver algorithm (Brent line search, tolerances), and the
  Nedelec sign convention (lexicographic edge orientation, which the tree gauge and the edge Dirichlet arrays rely
  on).
- **Windows not proven**: the MSYS2 build (`conda/build_windows.sh`: PETSc with `--download-hypre`, static MFEM,
  MS-MPI) has never run. See section 5.

### Equal in both

The same validation results to printed precision. The new patch tests give identical errors, including the same
1.5e-3 consistency error of Nedelec cells on non-affine hexahedra. MPI and threads match serial to about 1e-8.

## 4. Verification added: patch tests (pyTEC3D_APDL_benchmarks, formulations suite)

Four cases, 71 variants and 141 checks (`common/patch.py`, `pytec3d_driver.patch_test`). The field is uniform in two
layers (mur 2 and 7, the interface a mesh plane, the field oblique). Every boundary value is exact and the error
tolerance is 1e-10 (relative L2 of the cell H and B). The observed errors are 1e-15 to 6e-15 on both backends.

| case | what it checks |
|---|---|
| `patch3d_scalar` | PHI3D_nodal and TPHI3D_edge on TETRA10, HEX20, WEDGE15, PYRAMID13; sheared lattice, and interior nodes perturbed |
| `patch3d_vector` | AVA3D_edge (order 1 and 2, tree gauge) and AVA3D_nodal (order 1 and 2) on the four shapes; perturbed tets (edge) and all perturbed shapes (nodal); plus AVA3D_edge on non-affine hexes, with its known 1.5e-3 consistency error checked at 3e-3 |
| `patch_quadratic` | second-order nodal A with a linear B (A_z = 2cxy) on every shape and on perturbed tets and triangles (first order fails it at about 0.17, so the case can detect errors) |
| `patch2d` | AVA2D_nodal (order 1 and 2) and TPHI2D_edge on TRIANGLE6 and QUAD8, planar (sheared, perturbed) and axisymmetric |

Suggested further cases: linear magnetic diffusion in a slab against the exact backward-Euler modal solution,
time-order checks by Richardson extrapolation (backward Euler and theta = 1/2), and h-convergence rates on the ball
cases.

## 5. Windows (MFEM)

`conda/build_windows.sh` builds PETSc (static, MUMPS of MSYS2, `--download-hypre`), then static MFEM against that
hypre and the METIS of MSYS2, then the module. CI: `.github/workflows/conda.yml`, job `windows`.

The MFEM stage now runs CMake with gcc / g++ and explicit MS-MPI hints, not with the `mpicc` / `mpicxx` wrappers.
The MSYS2 wrappers are shell scripts, which the native CMake cannot use as compilers. It also passes the paths in
Windows form, since MFEMConfig.cmake records them. The static link order (MFEM -> hypre -> OpenBLAS) is covered by
`PETSC_STATIC_LDFLAGS` in `solver/CMakeLists.txt`. The solver uses only MPI-2 calls, which MS-MPI supports.

None of this has been run on Windows yet: the first CI run will show whether MFEM 4.9 builds with MinGW / MS-MPI.

## 6. MFEM capabilities added after the review (MFEM branch)

- **Linear solvers on the hypre matrix** (`solver_settings["linear"]["backend"] = "mfem"`): MUMPS through MFEM's
  `MUMPSSolver` (MFEM built with `MFEM_USE_MUMPS`, `solver/mfem/build_mfem.sh`). The analysis is reused on the same
  pattern, so the factorization counts match the PETSc path. On ITER_3D: the same time, 0.1 GB less memory (no PETSc
  copy of the matrix). MFEM's Krylov solvers (cg, gmres, fgmres, minres, bcgs) with hypre BoomerAMG, AMS, Jacobi or
  MUMPS as preconditioner.
- **AMS for curl-curl systems.** The discrete gradient comes from MFEM's `ParDiscreteLinearOperator`, mapped to the
  reduced numbering with the DOF orientations. An ungauged (singular) system is detected by a random probe of K·G on
  the interior vertices; AMS then runs singular, and the gradient part of the right-hand side is removed
  (f − G (GᵀG)⁻¹ Gᵀ f), so the system is consistent. On 287k edges: 12 CG iterations (independent of the mesh),
  0.6× the wall time and 0.6× the memory of MUMPS. A tree gauge spoils AMS (223 iterations): leave A ungauged with
  AMS. Tests: `tests/test_linear_backends.py`.
- **Arbitrary order on VTK Lagrange cells** (types 69–73, one order p per mesh). The geometry of order p is set as
  MFEM nodes (equispaced H1). The DOFs of order-p Lagrange fields map to the cell points through MFEM's own VTK
  Lagrange ordering (`CreateVTKElementConnectivity`). `option2` can be any order dividing p, for the nodal element
  types. Not possible: pyramids above order 2 (MFEM has no nodal high-order pyramid) and Nedelec above order 1/2
  (moment DOFs, not point values). Tests: the `patch_lagrange` case (58 checks; the polynomial patch is exact only
  in the complete order-p space).
- **A latent bug fixed:** several hypre products wrote into new, uninitialized vectors (y = 0·y + A x). hypre
  computes 0·y, so an inf or NaN left in fresh memory survived as NaN. The visible effect was an intermittent NaN
  right-hand side (`sphere_coil_3D_period`, about 1 run in 8 under load). The outputs are now zeroed before every
  such product.
- **Still not used, by design:** partial assembly, element assembly and device kernels (single-geometry meshes only
  in MFEM, no serendipity cells), and `mfem::NewtonSolver` (the libMesh-identical Newton, Brent line search
  included, is kept).
