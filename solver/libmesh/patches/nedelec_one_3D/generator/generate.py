"""Generate the NEDELEC_ONE shape function code that patch_nedelec_one_3D.py inserts into libMesh's
src/fe/fe_nedelec_one_shape_3D.C (the fragments ../first_*.C and ../second_*.C).

    FIRST order, PRISM15/18/20/21 and PYRAMID13/14/18: the shape functions of first_order.py
    SECOND order, TET14, HEX27, PRISM20/21 and PYRAMID18: the dual basis of the degrees of freedom
        edges       int phi.t l_k ds, l_0 = 1-s, l_1 = s along the edge (as libMesh in 2D)
        triangles   int phi.(x1-x0), int phi.(x2-x0) over the face
        quads       int phi.T_s (1-s), int phi.T_s s, int phi.T_t (1-t), int phi.T_t t
        interior    see spaces.interior_dofs()
      over the spaces of spaces.py (the pyramid's: MFEM's order 2 Fuentes basis). The face moments
      are for the side node order; the shape code combines them by the face orientation
      (../face_coefficients.C).

Every function is checked (dual to its degrees of freedom; the second order spaces contain the
linear fields) before the code is written. Requires sympy and numpy:

    python3 -m venv venv && venv/bin/pip install sympy numpy && venv/bin/python generate.py
"""
import os
import sys
import time

import numpy as np
import sympy as sp

from nedelec_ref import *
from spaces import SPACES, all_dofs, independent
from first_order import prism1, pyr1
from codegen import element_block, indent

OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))


def log(*a):
    print(*a, flush=True)


def check_contains_linear(el, funcs):
    """The span of funcs contains every linear vector field (needed for order two convergence)."""
    pts = {'TET': lambda a, b, c: (a / 3, b / 3, c / 3),
           'HEX': lambda a, b, c: (a, b, c),
           'PRISM': lambda a, b, c: (a / 2.2, b / 2.2, c),
           'PYRAMID': lambda a, b, c: (a * (1 - c) / 2, b * (1 - c) / 2, c / 2 + 0.5)}[el]
    samples = [pts(a, b, c) for a in (-0.7, -0.1, 0.4, 0.9) for b in (-0.8, 0.2, 0.6)
               for c in (-0.6, 0.0, 0.3, 0.8)]
    fn = [sp.lambdify(X, list(f), 'numpy') for f in funcs]
    M = np.array([[np.array(f(*p), dtype=float) for p in samples] for f in fn]).reshape(len(funcs), -1).T
    worst = 0.
    for m in (1, xi, eta, zeta):
        for e in range(3):
            q = [0, 0, 0]
            q[e] = m
            g = sp.lambdify(X, q, 'numpy')
            b = np.array([np.array(g(*p), dtype=float) * np.ones(3) for p in samples]).ravel()
            c, *_ = np.linalg.lstsq(M, b, rcond=None)
            worst = max(worst, np.abs(M @ c - b).max())
    assert worst < 1e-9, (el, worst)


def second_order(el):
    t0 = time.time()
    space = SPACES[el]()
    n = len(space)
    assert len(independent(space)) == n, el
    dofs = all_dofs(el)
    A = sp.Matrix(n, n, lambda i, j: dofs[i](space[j]))
    assert A.rank() == n, el + ': degrees of freedom not unisolvent'
    C = A.inv()
    funcs = [sp.Matrix([sp.cancel(sp.together(sum(C[j, i] * space[j][c] for j in range(n))))
                        for c in range(3)]) for i in range(n)]
    # dual: the degrees of freedom of the dual basis are the identity
    D = sp.Matrix(n, n, lambda i, j: dofs[i](funcs[j]))
    assert D == sp.eye(n), el + ': not dual'
    check_contains_linear(el, funcs)
    log('  %-8s SECOND: %d shape functions, dual to their degrees of freedom, contain P1^3 (%.0fs)'
        % (el, n, time.time() - t0))
    return funcs


def first_order(el, funcs):
    funcs = [sp.Matrix([sp.cancel(sp.together(c)) for c in f]) for f in funcs]
    check_dual(funcs, edge_dofs(el, 1))
    for c in (vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)):
        interp = sum((d(c) * f for d, f in zip(edge_dofs(el, 1), funcs)), sp.zeros(3, 1))
        assert sp.simplify(interp - c) == sp.zeros(3, 1), el
    log('  %-8s FIRST: %d shape functions, dual to -int phi.t on the edges' % (el, len(funcs)))
    return funcs


def main():
    log('Shape functions:')
    first = {'PRISM': first_order('PRISM', prism1), 'PYRAMID': first_order('PYRAMID', pyr1)}
    second = {el: second_order(el) for el in ('TET', 'HEX', 'PRISM', 'PYRAMID')}

    names = {'shape': 'shape', 'deriv': 'shape_deriv', 'second': 'shape_second_deriv'}
    titles = {'shape': 'shape functions', 'deriv': 'shape function first derivatives',
              'second': 'shape function second derivatives'}
    log('Code:')
    for what, name in names.items():
        # FIRST order: the prism and pyramid cases of switch (elem->type()), after the tetrahedra
        L = []
        for el in ('PRISM', 'PYRAMID'):
            L += element_block(el, 1, first[el], what, 'i') + ['']
        path = os.path.join(OUT, 'first_%s.C' % name)
        with open(path, 'w') as f:
            f.write('\n'.join(indent(L, 10)).rstrip() + '\n\n')
        log('  ' + os.path.relpath(path, OUT))

        # SECOND order: the case SECOND of switch (totalorder), before the unsupported orders
        L = ['      // quadratic Nedelec (first kind) %s' % titles[what],
             '    case SECOND:',
             '      {',
             '        switch (elem->type())',
             '          {']
        B = []
        for el in ('TET', 'HEX', 'PRISM', 'PYRAMID'):
            B += element_block(el, 2, second[el], what, 'ii') + ['']
        L += indent(B, 10)
        L += ['          default:',
              '            libmesh_error_msg("ERROR: Unsupported 3D element type!: " << '
              'Utility::enum_to_string(elem->type()));',
              '          }',
              '      }',
              '', '']
        path = os.path.join(OUT, 'second_%s.C' % name)
        with open(path, 'w') as f:
            f.write('\n'.join(L))
        log('  ' + os.path.relpath(path, OUT))


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    main()
