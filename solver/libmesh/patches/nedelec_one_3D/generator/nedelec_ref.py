"""Reference-element machinery for NEDELEC_ONE in 3D: spaces, degrees of freedom, dual bases.

Coordinates are libMesh's reference coordinates (xi, eta, zeta).
"""
import sympy as sp

xi, eta, zeta = sp.symbols('xi eta zeta', real=True)
X = (xi, eta, zeta)
s_, t_, u_, v_ = sp.symbols('s t u v', real=True)

R = sp.Rational


def vec(*c):
    return sp.Matrix([sp.sympify(a) for a in c])


def grad(f):
    return sp.Matrix([sp.diff(f, x) for x in X])


# ----------------------------------------------------------------------------- geometry
ELEMS = {
    'TET': dict(
        verts=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
        edges=[(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
        sides=[(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)],
        # sides in the order of their face nodes (libMesh dof order)
        face_order=[0, 1, 2, 3]),
    'HEX': dict(
        verts=[(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
               (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)],
        edges=[(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 5), (2, 6), (3, 7),
               (4, 5), (5, 6), (6, 7), (4, 7)],
        sides=[(0, 3, 2, 1), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7)],
        face_order=[0, 1, 2, 3, 4, 5]),
    'PRISM': dict(
        verts=[(0, 0, -1), (1, 0, -1), (0, 1, -1), (0, 0, 1), (1, 0, 1), (0, 1, 1)],
        edges=[(0, 1), (1, 2), (0, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (3, 5)],
        sides=[(0, 2, 1), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5), (3, 4, 5)],
        # face nodes 15, 16, 17 (quads 1, 2, 3), then 18, 19 (triangles 0, 4)
        face_order=[1, 2, 3, 0, 4]),
    'PYRAMID': dict(
        verts=[(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0), (0, 0, 1)],
        edges=[(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)],
        sides=[(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4), (0, 3, 2, 1)],
        # face node 13 (quad 4), then 14..17 (triangles 0..3)
        face_order=[4, 0, 1, 2, 3]),
}


def V(el, k):
    return sp.Matrix(ELEMS[el]['verts'][k])


def at(f, p):
    """Substitute the point p (3 sympy expressions) into the vector or scalar f."""
    return f.subs({xi: p[0], eta: p[1], zeta: p[2]}, simultaneous=True)


# ----------------------------------------------------------------------------- integration
def integrate_volume(el, f):
    """Exact integral of the scalar f over the reference element el."""
    f = sp.sympify(f)
    if el == 'TET':
        return sp.integrate(f, (zeta, 0, 1 - xi - eta), (eta, 0, 1 - xi), (xi, 0, 1))
    if el == 'HEX':
        return sp.integrate(f, (xi, -1, 1), (eta, -1, 1), (zeta, -1, 1))
    if el == 'PRISM':
        return sp.integrate(f, (eta, 0, 1 - xi), (xi, 0, 1), (zeta, -1, 1))
    if el == 'PYRAMID':
        # collapsed coordinates: xi = (1-zeta) a, eta = (1-zeta) b
        a, b, w = sp.symbols('a b w', real=True)
        g = sp.cancel(sp.together(f.subs({xi: (1 - w) * a, eta: (1 - w) * b, zeta: w},
                                         simultaneous=True) * (1 - w) ** 2))
        return sp.integrate(sp.expand(g), (a, -1, 1), (b, -1, 1), (w, 0, 1))
    raise ValueError(el)


def integrate_edge(f, A, B, weight):
    """int_0^1 f(A + s (B-A)) . (B-A) weight(s) ds"""
    p = A + s_ * (B - A)
    g = sp.cancel(sp.together((at(f, p).T * (B - A))[0] * weight))
    return sp.integrate(sp.expand(g), (s_, 0, 1))


def integrate_tri(f, P0, P1, P2, tangent, weight=1):
    """int over the unit triangle of f(P0 + u (P1-P0) + v (P2-P0)) . tangent weight du dv"""
    p = P0 + u_ * (P1 - P0) + v_ * (P2 - P0)
    g = sp.cancel(sp.together((at(f, p).T * tangent)[0] * weight))
    return sp.integrate(sp.expand(g), (v_, 0, 1 - u_), (u_, 0, 1))


def integrate_quad(f, P0, P1, P3, tangent, weight=1):
    """int over the unit square of f(P0 + s (P1-P0) + t (P3-P0)) . tangent weight ds dt"""
    p = P0 + s_ * (P1 - P0) + t_ * (P3 - P0)
    g = sp.cancel(sp.together((at(f, p).T * tangent)[0] * weight))
    return sp.integrate(sp.expand(g), (s_, 0, 1), (t_, 0, 1))


# ----------------------------------------------------------------------------- degrees of freedom
def edge_dofs(el, order):
    """Edge moments against the Lagrange basis of degree order-1 along each edge, in the
    direction of the edge (local edge node 0 to 1). FIRST order follows libMesh's sign: -int."""
    out = []
    for (a, b) in ELEMS[el]['edges']:
        A, B = V(el, a), V(el, b)
        if order == 1:
            out.append(lambda f, A=A, B=B: -integrate_edge(f, A, B, 1))
        elif order == 2:
            out.append(lambda f, A=A, B=B: integrate_edge(f, A, B, 1 - s_))
            out.append(lambda f, A=A, B=B: integrate_edge(f, A, B, s_))
        else:
            raise ValueError(order)
    return out


def face_dofs(el, side):
    """Second order face moments, in the side node order of the face (the 'local' orientation).

    Triangle (v0, v1, v2):     int phi.(v1-v0), int phi.(v2-v0)
    Quadrilateral (v0..v3):    int phi.T_s (1-s), int phi.T_s s, int phi.T_t (1-t), int phi.T_t t
                               with T_s = v1-v0, T_t = v3-v0, (s,t) in the unit square
    """
    vs = [V(el, k) for k in ELEMS[el]['sides'][side]]
    if len(vs) == 3:
        P0, P1, P2 = vs
        return [lambda f: integrate_tri(f, P0, P1, P2, P1 - P0),
                lambda f: integrate_tri(f, P0, P1, P2, P2 - P0)]
    P0, P1, P2, P3 = vs
    assert P2 == P1 + P3 - P0, 'non-parallelogram face'
    Ts, Tt = P1 - P0, P3 - P0
    return [lambda f: integrate_quad(f, P0, P1, P3, Ts, 1 - s_),
            lambda f: integrate_quad(f, P0, P1, P3, Ts, s_),
            lambda f: integrate_quad(f, P0, P1, P3, Tt, 1 - t_),
            lambda f: integrate_quad(f, P0, P1, P3, Tt, t_)]


def dual_basis(space, dofs):
    """The functions of span(space) dual to dofs."""
    n = len(space)
    assert len(dofs) == n, (len(dofs), n)
    A = sp.Matrix(n, n, lambda i, j: dofs[i](space[j]))
    assert A.rank() == n, 'degrees of freedom not unisolvent'
    C = A.inv()
    return [sp.Matrix([sp.cancel(sp.together(sum(C[j, i] * space[j][c] for j in range(n))))
                       for c in range(3)]) for i in range(n)]


def check_dual(funcs, dofs):
    n = len(funcs)
    for i in range(n):
        for j in range(n):
            val = sp.simplify(dofs[i](funcs[j]))
            assert val == (1 if i == j else 0), (i, j, val)
    return True
