# NEDELEC_ONE on 3D cells: the shape function code of the libMesh patch

`../patch_nedelec_one_3D.py` inserts these fragments into libMesh's `src/fe/fe_nedelec_one_shape_3D.C`
(and edits `fe_nedelec_one.C`, `fe_interface.C` itself):

| order  | cells                                    | upstream libMesh |
|--------|------------------------------------------|------------------|
| FIRST  | TET10/14, HEX20/27                       | yes              |
| FIRST  | PRISM15/18/20/21, PYRAMID13/14/18        | patch            |
| SECOND | TET14, HEX27, PRISM20/21, PYRAMID18      | patch            |

Second order puts degrees of freedom on the faces, so it needs a node at the center of every face:
TET10, HEX20, PRISM15/18 and PYRAMID13/14 stay first order (`max_order()` is 1 for them).

The patch also fixes two sign errors of upstream's FIRST order hexahedron second derivatives.

## Status in pyTEC3D: second order is not used yet

The patch provides and tests second order shapes, but the pyTEC3D solver supports **first order NEDELEC_ONE only**
for now. It assumes one edge DOF per edge midside node, and second order breaks that in several places:

1. **Mesh import** (`solver/src/array_io.cpp`, `elem_type_of`): of the cells with face-center nodes only HEX27
   (`TRIQUADRATIC_HEXAHEDRON`) can be imported. VTK has no fixed cell types matching TET14, PRISM20/21 or
   PYRAMID18 (`VTK_TRIQUADRATIC_PYRAMID` has 19 nodes). Upgrading the cells in the solver (`all_complete_order`)
   would leave the new face nodes without point-array entries (BCs, loads, outputs).
2. **Dirichlet constraints** (`solver/src/constraints.cpp`): an edge variable takes one value per point (component 2
   of `_F`) and constrains only `dof_number(..., 0)`. The second edge DOF and the face DOFs would stay free. Also,
   the prescribed value is a circulation, while the second order DOFs are moments against (1-s, s) plus face
   moments: a boundary value would have to be a tangential field projected onto them.
3. **Couplings / periodicity** (`constraints.cpp`, `build_edge_orientations_map`): an edge pair maps `dof 0 -> dof 0`
   with sign +-1. At second order an opposite orientation also swaps the two edge DOFs, and coupled faces need the
   face orientation mapping (combinations on triangles, signed permutations on quadrilaterals). Master and slave
   faces can be oriented differently, since the orientation follows the lexicographic order of the points.
4. **Tree-cotree gauge** (`python/src/pytec3d/tree_gauge.py`): one DOF per edge, identified by its midside point.
   It removes only the lowest order gradients. At second order the higher order gradients (edges, faces, interior)
   must be gauged too, which is usually done with a *hierarchical* basis. The basis here is dual (interpolatory,
   libMesh's 2D convention) and does not separate them. If a gauge is needed at second order, a hierarchical
   second order basis may be the better choice.
5. **Nodal loads** (`solver/src/element_builder.cpp`, `pytec3d_system.cpp`): `dof 0` only.
6. **Python edge helpers** (`tree_gauge.EDGE_TABLES`, `edges()`, `edge_circulations()`): the quadratic VTK cells
   only, one circulation per edge.

Assembly, element DOF offsets and the nodal solution output (`nodal_soln`) are already independent of the order.

## Files

* `face_coefficients.C`: written by hand. Maps the face degrees of freedom to the face orientation,
  `Elem::oriented_face_vertices()`, so that all elements sharing a face agree on them.
* `first_*.C`, `second_*.C`: **generated** by `generator/generate.py`. Do not edit them by hand.

## Degrees of freedom (second order)

Element dof order is libMesh's: edges (edge node order), then faces (face node order), then interior.

* Edges: `int phi.t l_k ds`, with `t` the edge vector and `l_0 = 1-s`, `l_1 = s`. This is the
  convention of libMesh's 2D SECOND order. An edge with negative orientation reverses its two dofs
  and changes their sign (`sign`, `ii`), as in 2D.
* Triangles: `int phi.(x1-x0)`, `int phi.(x2-x0)`. Quadrilaterals: `int phi.T_s (1-s)`,
  `int phi.T_s s`, `int phi.T_t (1-t)`, `int phi.T_t t`. All are taken over the face, parametrized
  from its oriented vertices.
* Interior:
  * HEX: each component against `(1-x_c), x_c` along its own direction.
  * PRISM: `int phi.e_x`, `int phi.e_y`.
  * PYRAMID: against `e_x, e_y, e_z, xi e_x, eta e_y, zeta e_z`.

Spaces:
* TET: Nedelec first kind.
* HEX: Q(1,2,2) x Q(2,1,2) x Q(2,2,1).
* PRISM: N1_2(triangle) x P2(zeta) horizontally, plus P2(triangle) x P1(zeta) vertically.
* PYRAMID: MFEM's order-2 Fuentes space (rational, 34 functions).

Every second order space contains the linear fields.

## Regenerating

```
python3 -m venv venv && venv/bin/pip install sympy numpy
venv/bin/python generator/generate.py      # about 3 minutes; checks duality and P1 reproduction
```

## Tests

`solver/libmesh/tests/run_nedelec_one_3D_tests.sh <libMesh prefix> <patched libMesh source>` runs:
* derivatives against finite differences;
* tangential continuity, on single-type and hybrid meshes;
* curl-curl h-convergence (rates 1 and 2);
* second derivatives.
