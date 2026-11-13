"""Second order NEDELEC_ONE spaces and degrees of freedom on the libMesh reference cells."""
import sympy as sp
from nedelec_ref import *

h = R(1, 2)


def independent(funcs):
    """A maximal linearly independent subset of the (rational) vector functions funcs
    (numerical rank at sample points)."""
    import numpy as np
    pts = [(-0.31 + 0.17 * i, -0.27 + 0.13 * j, 0.07 + 0.11 * k) for i in range(4)
           for j in range(4) for k in range(5)]
    out, rows = [], []
    for f in funcs:
        fn = sp.lambdify(X, list(f), 'numpy')
        row = np.array([fn(*p) for p in pts], dtype=float).ravel()
        M = np.array(rows + [row])
        if np.linalg.matrix_rank(M, tol=1e-9 * np.abs(M).max()) > len(rows):
            rows.append(row)
            out.append(f)
    return out


# ----------------------------------------------------------------------------- spaces
def tet2():
    P1 = [m * e for m in (1, xi, eta, zeta) for e in (vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1))]
    x = vec(xi, eta, zeta)
    S2 = [x.cross(m * e) for m in (xi, eta, zeta) for e in (vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1))]
    sp_ = independent(P1 + S2)
    assert len(sp_) == 20
    return sp_


def hex2():
    out = []
    for c in range(3):
        degs = [2, 2, 2]
        degs[c] = 1
        for a in range(degs[0] + 1):
            for b in range(degs[1] + 1):
                for d in range(degs[2] + 1):
                    e = [0, 0, 0]
                    e[c] = xi ** a * eta ** b * zeta ** d
                    out.append(vec(*e))
    assert len(out) == 54
    return out


def prism2():
    # N1_2 on the triangle (xi, eta)
    nd = [m * e for m in (1, xi, eta) for e in (vec(1, 0, 0), vec(0, 1, 0))]
    nd += [m * vec(-eta, xi, 0) for m in (xi, eta)]
    hor = [f * zeta ** k for f in nd for k in range(3)]
    p2 = [xi ** a * eta ** b for a in range(3) for b in range(3) if a + b <= 2]
    ver = [vec(0, 0, m * zeta ** k) for m in p2 for k in range(2)]
    out = hor + ver
    assert len(independent(out)) == 36
    return out


def pyramid2():
    """MFEM's ND_FuentesPyramidElement basis of order 2 (raw, before its nodal transformation),
    mapped covariantly to the libMesh pyramid."""
    x, y, z = sp.symbols('x y z', real=True)
    XM = (x, y, z)

    def g(f):
        return sp.Matrix([sp.diff(f, c) for c in XM])

    xy = [x, y]
    lam1 = (1 - x - z) * (1 - y - z) / (1 - z)
    lam2 = x * (1 - y - z) / (1 - z)
    lam3 = x * y / (1 - z)
    lam4 = (1 - x - z) * y / (1 - z)
    lam5 = z
    mu0z, mu1z = 1 - z, z

    def mu0(ab): return 1 - xy[ab - 1] / (1 - z)
    def mu1(ab): return xy[ab - 1] / (1 - z)
    def nu0(ab): return 1 - xy[ab - 1] - z
    def nu1(ab): return xy[ab - 1]
    def nu2(ab): return z

    def sds(s0, s1):  # s0 grad s1 - s1 grad s0
        return s0 * g(s1) - s1 * g(s0)

    # homogenized scaled Legendre, degrees 0..1:  P0 = 1, P1 = s1 - s0
    def P(s0, s1): return [1, s1 - s0]
    # homogenized integrated Legendre, degrees 0..2: 0, s1, -s0 s1
    def phiE(s0, s1): return [0, s1, -s0 * s1]

    def E_E(p, s, sd):
        Ps = P(*s)
        return [Ps[i] * sd for i in range(p)]

    W = []
    p = 2
    # mixed edges
    EE = E_E(p, (nu0(1), nu1(1)), sds(nu0(1), nu1(1)))
    W += [mu0(2) * f for f in EE] + [mu1(2) * f for f in EE]
    EE = E_E(p, (nu0(2), nu1(2)), sds(nu0(2), nu1(2)))
    W += [mu0(1) * f for f in EE] + [mu1(1) * f for f in EE]
    # triangle edges
    for lam in (lam1, lam2, lam3, lam4):
        W += E_E(p, (lam, lam5), sds(lam, lam5))
    # quadrilateral face
    mu2 = mu0z ** 2
    def E_Q(s, sd, t):
        EEi = E_E(p, s, sd)
        ph = phiE(*t)
        return {(i, j): ph[j] * EEi[i] for i in range(p) for j in range(2, p + 1)}
    EQ1 = E_Q((mu0(1), mu1(1)), sds(mu0(1), mu1(1)), (mu0(2), mu1(2)))
    EQ2 = E_Q((mu0(2), mu1(2)), sds(mu0(2), mu1(2)), (mu0(1), mu1(1)))
    for EQ in (EQ1, EQ2):
        for j in range(2, p + 1):
            for i in range(p):
                W.append(mu2 * EQ[(i, j)])
    # triangular faces: E_T(p=2) = s2 * (s0 grad s1 - s1 grad s0) for (s0, s1, s2)
    def E_T(s, sd): return s[2] * sd
    for fam in (1, 2):
        for ab, mu in ((1, mu0(2)), (1, mu1(2)), (2, mu0(1)), (2, mu1(1))):
            if fam == 1:
                s = (nu0(ab), nu1(ab), nu2(ab))
                sd = sds(nu0(ab), nu1(ab))
            else:
                s = (nu1(ab), nu2(ab), nu0(ab))
                sd = sds(nu1(ab), nu2(ab))
            W.append(mu * E_T(s, sd))
    # interior
    phQ1 = phiE(mu0(1), mu1(1))[2] * phiE(mu0(2), mu1(2))[2]
    phEz = phiE(mu0z, mu1z)[2]
    W.append(g(phQ1 * phEz))                     # family I
    for i in range(p):
        W.append(mu0z * EQ1[(i, 2)] * phEz)      # family II
    for i in range(p):
        W.append(mu0z * EQ2[(i, 2)] * phEz)      # family III
    phQ2 = phiE(mu0(2), mu1(2))[2] * phiE(mu0(1), mu1(1))[2]
    W.append(2 * mu0z * phQ2 * g(mu0z))          # family IV
    assert len(W) == 34

    # covariant map to libMesh: phi_lm = J^T phi_mfem(x(xi))
    sub = {x: h * (1 + xi - zeta), y: h * (1 + eta - zeta), z: zeta}
    out = []
    for f in W:
        fx, fy, fz = [sp.cancel(sp.together(c.subs(sub, simultaneous=True))) for c in f]
        out.append(sp.Matrix([sp.cancel(h * fx), sp.cancel(h * fy),
                              sp.cancel(-h * fx - h * fy + fz)]))
    return out


SPACES = {'TET': tet2, 'HEX': hex2, 'PRISM': prism2, 'PYRAMID': pyramid2}


# ----------------------------------------------------------------------------- interior dofs
def interior_dofs(el):
    if el == 'TET':
        return []
    if el == 'HEX':
        # each component against the Lagrange basis in its own direction
        out = []
        for c in range(3):
            xc = h * (X[c] + 1)
            for w in (1 - xc, xc):
                out.append(lambda f, c=c, w=w: integrate_volume('HEX', f[c] * w))
        return out
    if el == 'PRISM':
        return [lambda f: integrate_volume('PRISM', f[0]),
                lambda f: integrate_volume('PRISM', f[1])]
    if el == 'PYRAMID':
        # moments against these fields, a set invariant under the symmetries of the pyramid
        # (the same with eta e_x, xi e_y in place of xi e_x, eta e_y is not unisolvent)
        qs = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1), vec(xi, 0, 0), vec(0, eta, 0), vec(0, 0, zeta)]
        return [lambda f, q=q: integrate_volume('PYRAMID', (f.T * q)[0]) for q in qs]
    raise ValueError(el)


def all_dofs(el):
    """The second order degrees of freedom, in libMesh's element order: edges (in the edge node
    order), faces (in the face node order), interior."""
    d = edge_dofs(el, 2)
    for s in ELEMS[el]['face_order']:
        d += face_dofs(el, s)
    return d + interior_dofs(el)
