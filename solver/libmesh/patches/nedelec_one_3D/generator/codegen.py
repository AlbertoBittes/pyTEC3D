"""The C++ switch blocks of fe_nedelec_one_shape_3D.C, in the style of libMesh's shape files."""
import os, pickle, sys
import sympy as sp
from nedelec_ref import *
from emit import cxx, nice, num

# derivative index j of shape_second_deriv -> (variable a, variable b)
SECOND_J = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]
SECOND_J_COMMENT = ['d^2()/dxi^2', 'd^2()/dxideta', 'd^2()/deta^2',
                    'd^2()/dxidzeta', 'd^2()/detadzeta', 'd^2()/dzeta^2']
DERIV_COMMENT = ['d()/dxi', 'd()/deta', 'd()/dzeta']

# --------------------------------------------------------------------------- printing components
def comp_text(e, rational):
    e = sp.cancel(sp.together(sp.sympify(e)))
    if e == 0:
        return '0.0'
    if not rational:
        return cxx(nice(e))
    n, d = sp.fraction(e)
    k = sp.degree(d, zeta) if d.has(zeta) else 0
    c = sp.cancel(d / (1 - zeta) ** k)
    assert not c.free_symbols, (e, d)
    n = sp.factor(sp.expand(n / c))
    if k == 0:
        return cxx(n)

    def over(numer, k):
        dtext = '*'.join(['den'] * k)
        t = cxx(numer)
        # a bare product or single factor needs no parentheses before the division
        if numer.is_Add:
            t = '(' + t + ')'
        return t + '/' + (dtext if k == 1 else '(' + dtext + ')')

    best = over(n, k)
    # Or a polynomial part plus a remainder of lower degree in den = 1-zeta
    d = sp.Symbol('d')
    nd = sp.Poly(sp.expand(n.subs(zeta, 1 - d)), d)
    rem = sum((c_ * d ** r for (r,), c_ in nd.terms() if r < k), sp.Integer(0))
    quo = sp.expand((nd.as_expr() - rem) / d ** k).subs(d, 1 - zeta)
    rem = sp.factor(sp.expand(rem.subs(d, 1 - zeta)))
    quo = sp.factor(sp.expand(quo))
    if quo != 0 and rem != 0:
        r = over(rem, k)
        alt = cxx(quo) + ('-' + r[1:] if r.startswith('-') else '+' + r)
        if len(alt) < len(best):
            best = alt
    return best


def grad(v, rational):
    return 'RealGradient( ' + ', '.join(comp_text(c, rational) for c in v) + ' )'


def zero(v):
    return all(sp.cancel(sp.together(c)) == 0 for c in v)


# --------------------------------------------------------------------------- element layouts
TYPES = {  # libMesh element type cases, per element, at SECOND order
    'TET': ['TET14'], 'HEX': ['HEX27'], 'PRISM': ['PRISM20', 'PRISM21'], 'PYRAMID': ['PYRAMID18']}

ASSERTS = {
    'HEX': ['// Even with a loose inverse_map tolerance we ought to',
            '// be nearly on the element interior in master',
            '// coordinates',
            'libmesh_assert_less_equal ( std::fabs(xi),   1.0+10*TOLERANCE );',
            'libmesh_assert_less_equal ( std::fabs(eta),  1.0+10*TOLERANCE );',
            'libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );'],
    'PRISM': ['// Even with a loose inverse_map tolerance we ought to',
              '// be nearly on the element interior in master',
              '// coordinates',
              'libmesh_assert_greater_equal ( xi,  -10*TOLERANCE );',
              'libmesh_assert_greater_equal ( eta, -10*TOLERANCE );',
              'libmesh_assert_less_equal ( xi+eta, 1.0+10*TOLERANCE );',
              'libmesh_assert_less_equal ( std::fabs(zeta), 1.0+10*TOLERANCE );'],
    'PYRAMID': ['// Even with a loose inverse_map tolerance we ought to',
                '// be nearly on the element interior in master',
                '// coordinates',
                'libmesh_assert_greater_equal ( zeta, -10*TOLERANCE );',
                'libmesh_assert_less_equal ( zeta, 1.0+10*TOLERANCE );',
                'libmesh_assert_less_equal ( std::fabs(xi),  1.0-zeta+10*TOLERANCE );',
                'libmesh_assert_less_equal ( std::fabs(eta), 1.0-zeta+10*TOLERANCE );'],
    'TET': [],
}


def layout(el, order):
    """[(kind, first dof, n dofs, side)] in libMesh element dof order"""
    ne = len(ELEMS[el]['edges'])
    if order == 1:
        return [('edges', 0, ne, None)]
    out = [('edges', 0, 2 * ne, None)]
    k = 2 * ne
    for s in ELEMS[el]['face_order']:
        n = 2 * len(ELEMS[el]['sides'][s]) - 4  # 2 on a triangle, 4 on a quadrilateral
        out.append(('face', k, n, s))
        k += n
    n_int = {'TET': 0, 'HEX': 6, 'PRISM': 2, 'PYRAMID': 6}[el]
    if n_int:
        out.append(('interior', k, n_int, None))
    return out


# --------------------------------------------------------------------------- switch over i
def switch_i(el, order, funcs, ind, var):
    """The switch over the shape function index of funcs (vectors already differentiated)."""
    rational = el == 'PYRAMID'
    L = []
    w = L.append
    w('switch(%s)' % var)
    w('  {')
    # the zero shape functions, grouped at the end as libMesh does
    zeros = []
    for kind, first, n, side in layout(el, order):
        if kind == 'face':
            terms = [(m, funcs[first + m]) for m in range(n) if not zero(funcs[first + m])]
            if not terms:
                zeros += range(first, first + n)
                continue
            for k in range(first, first + n):
                w('  case %d:' % k)
            w('    {')
            w('      const std::array<Real, 4> c =')
            w('        nedelec_one_face_coefficients(*elem, %d, i - %d);' % (side, first))
            for t, (m, f) in enumerate(terms):
                lead = '      return ' if t == 0 else '             '
                tail = ';' if t == len(terms) - 1 else ' +'
                w(lead + 'c[%d] * %s%s' % (m, grad(f, rational), tail))
            w('    }')
            continue
        pre = 'sign * ' if kind == 'edges' else ''
        for k in range(first, first + n):
            if zero(funcs[k]):
                zeros.append(k)
                continue
            w('  case %d:' % k)
            w('    return %s%s;' % (pre, grad(funcs[k], rational)))
    for k in zeros:
        w('  case %d:' % k)
    if zeros:
        w('    return RealGradient();')
    w('')
    w('  default:')
    w('    libmesh_error_msg("Invalid i = " << i);')
    w('  } // switch(%s)' % var)
    return ['  ' * ind + l if l else l for l in L]


def indent(lines, n):
    return [(' ' * n + l) if l else l for l in lines]


# --------------------------------------------------------------------------- pyramid apex
APEX_TOL = '1.e-6'


def apex_limits(funcs):
    out = []
    for f in funcs:
        out.append([sp.limit(sp.cancel(c.subs({xi: 0, eta: 0})), zeta, 1, '-') for c in f])
    return out


# --------------------------------------------------------------------------- element blocks
def element_block(el, order, funcs, what, var):
    """The block `case <types>: { ... }` of one element, for shape / shape_deriv / shape_second_deriv."""
    types = TYPES[el] if order == 2 else {'PRISM': ['PRISM15', 'PRISM18', 'PRISM20', 'PRISM21'],
                                          'PYRAMID': ['PYRAMID13', 'PYRAMID14', 'PYRAMID18']}[el]
    rational = el == 'PYRAMID'
    L = ['case %s:' % t for t in types]
    L.append('  {')
    body = list(ASSERTS[el])
    if rational:
        body.append('')
    if what == 'shape':
        if rational:
            body += ['// The shape functions are rational in den = 1-zeta. At the apex,',
                     '// where den vanishes, they only have limits along each direction:',
                     '// take those along the pyramid axis.',
                     'if (1.0 - zeta < %s)' % APEX_TOL,
                     '  {']
            lim = [sp.Matrix(v) for v in apex_limits(funcs)]
            sw = switch_i(el, order, lim, 0, var)
            body += indent(sw, 4) + ['  }', '', 'const Real den = 1.0 - zeta;', '']
        body += switch_i(el, order, funcs, 0, var)
    elif what == 'deriv':
        if rational:
            body += ['// The derivatives of the rational shape functions have no limit at',
                     '// the apex: only keep them finite there.',
                     'const Real den = std::max(1.0 - zeta, %s);' % APEX_TOL, '']
        body.append('switch (j)')
        body.append('  {')
        for j in range(3):
            d = [sp.diff(f, X[j]) for f in funcs]
            body.append('    // %s' % DERIV_COMMENT[j])
            body.append('  case %d:' % j)
            if all(zero(f) for f in d):
                body.append('    return RealGradient();')
                body.append('')
                continue
            body.append('    {')
            body += switch_i(el, order, d, 3, var)
            body.append('    } // j = %d' % j)
            body.append('')
        body += ['  default:', '    libmesh_error_msg("Invalid j = " << j);', '  }']
    elif what == 'second':
        if rational:
            body += ['// The derivatives of the rational shape functions have no limit at',
                     '// the apex: only keep them finite there.',
                     'const Real den = std::max(1.0 - zeta, %s);' % APEX_TOL, '']
        d2 = [[sp.diff(f, X[a], X[b]) for f in funcs] for (a, b) in SECOND_J]
        zero_j = [j for j in range(6) if all(zero(f) for f in d2[j])]
        body.append('switch (j)')
        body.append('  {')
        if zero_j:
            body.append('    // %s' % ', '.join(SECOND_J_COMMENT[j] for j in zero_j))
            for j in zero_j:
                body.append('  case %d:' % j)
            body.append('    return RealGradient();')
            body.append('')
        for j in range(6):
            if j in zero_j:
                continue
            body.append('    // %s' % SECOND_J_COMMENT[j])
            body.append('  case %d:' % j)
            body.append('    {')
            body += switch_i(el, order, d2[j], 3, var)
            body.append('    } // j = %d' % j)
            body.append('')
        body += ['  default:', '    libmesh_error_msg("Invalid j = " << j);', '  }']
    L += indent(body, 4)
    L.append('  }')
    return L


def kind_of(el, order, k):
    for kind, first, n, side in layout(el, order):
        if first <= k < first + n:
            return kind




# --------------------------------------------------------------------------- main
