"""Print sympy expressions and switch blocks in the style of libMesh's fe_*_shape_*.C files."""
import sympy as sp
from fractions import Fraction


def num(c):
    """A rational number as a C++ double literal."""
    c = sp.Rational(c)
    f = Fraction(int(c.p), int(c.q))
    if f.denominator == 1:
        return '%d.0' % f.numerator
    d = f.denominator
    while d % 2 == 0:
        d //= 2
    while d % 5 == 0:
        d //= 5
    if d == 1:  # terminating decimal
        s = ('%.12f' % float(f)).rstrip('0')
        return s
    return '%d.0/%d.0' % (f.numerator, f.denominator)


def _factors(m):
    """Split a Mul into (rational coefficient, numerator factors, denominator factors)."""
    coeff, rest = m.as_coeff_Mul()
    numer, denom = [], []
    for f in sp.Mul.make_args(rest):
        if f.is_Pow and f.exp.is_Integer and f.exp < 0:
            denom.append(sp.Pow(f.base, -f.exp))
        else:
            numer.append(f)
    return sp.Rational(coeff), numer, denom


def cxx(e, symbols_order=None):
    """C++ text of a polynomial/rational expression e."""
    e = sp.sympify(e)
    return _p(e, top=True)


def _pow(b, n):
    bs = _p(b, factor=True)
    return '*'.join([bs] * n)


def _p(e, top=False, factor=False):
    if e.is_Rational:
        s = num(e)
        return '(' + s + ')' if (factor and (e < 0 or '/' in s)) else s
    if e.is_Symbol:
        return e.name
    if e.is_Pow:
        if e.exp.is_Integer and e.exp > 0:
            return _pow(e.base, int(e.exp))
        if e.exp.is_Integer and e.exp < 0:
            return '1.0/(' + _pow(e.base, int(-e.exp)) + ')'
        raise ValueError(e)
    if e.is_Add:
        terms = e.as_ordered_terms()
        out = ''
        for k, t in enumerate(terms):
            c, _ = t.as_coeff_Mul()
            if c < 0:
                s = _p(-t)
                out += ('-' if k == 0 else '-') + s
            else:
                s = _p(t)
                out += ('' if k == 0 else '+') + s
        return '(' + out + ')' if factor or not top else out
    if e.is_Mul:
        c, numer, denom = _factors(e)
        neg = c < 0
        c = abs(c)
        parts = []
        if c != 1 or not numer:
            parts.append(num(c))
        parts += [_p(f, factor=True) for f in numer]
        s = '*'.join(parts)
        if denom:
            d = '*'.join(_p(f, factor=True) for f in denom)
            s += '/' + (d if len(denom) == 1 and not denom[0].is_Add else '(' + d + ')')
        if neg:
            s = '-' + s
            if factor:
                s = '(' + s + ')'
        return s
    raise ValueError('cannot print %r' % e)


def nice(e):
    """Factor a polynomial, keep a rational function as factored numerator over denominator."""
    e = sp.cancel(sp.together(sp.sympify(e)))
    if e == 0:
        return sp.Integer(0)
    n, d = sp.fraction(e)
    n = sp.factor(n)
    d = sp.factor(d)
    return n / d if d != 1 else n


def grad_text(v, pre=''):
    comps = [cxx(nice(c)) for c in v]
    return pre + 'RealGradient( ' + ', '.join(comps) + ' )'


def is_zero(v):
    return all(sp.simplify(c) == 0 for c in v)
