# Patch selected libMesh FE sources in place (prism/pyramid support for NEDELEC_ONE).
import re, sys, pathlib

root = pathlib.Path(sys.argv[1])

# ---------------------------------------------------------------- helpers
SEP = r'(?:\s|//[^\n]*(?:\n|$))'   # whitespace or // comments


def _isword(c):
    return c.isalnum() or c == '_'


def ws(snippet):
    """Regex matching `snippet` regardless of whitespace/newlines/comments between tokens."""
    toks = re.findall(r'\w+|[^\w\s]', snippet)
    out = []
    for i, t in enumerate(toks):
        if i:
            need = _isword(toks[i - 1][-1]) and _isword(t[0])
            out.append(SEP + ('+' if need else '*'))
        out.append(re.escape(t))
    return ''.join(out)


def die(msg):
    print('ERROR (libMesh patch): ' + msg, file=sys.stderr)
    sys.exit(1)


def apply(rel, guard, edits):
    path = root / rel
    if not path.is_file():
        die(f'{rel} not found')
    text = path.read_text()
    if re.search(guard, text):
        print(f'WARNING: {rel} already seems to support prism/pyramid Nedelec upstream; skipping.')
        return
    for name, pattern, block, where in edits:
        ms = list(re.finditer(pattern, text))
        if len(ms) != 1:
            die(f'{rel}: anchor "{name}" matched {len(ms)} times (expected 1). '
                'Upstream code changed; update the patch.')
        m = ms[0]
        if where == 'after':
            text = text[:m.end()] + block + text[m.end():]
        elif where == 'before':
            text = text[:m.start()] + block + text[m.start():]
        elif where == 'replace':
            text = text[:m.start()] + block + text[m.end():]
        else:
            die('bad where')
        print(f'  patched {rel}: {name}')
    path.write_text(text)


# ------------------------------------------------------------ fe_interface.C
apply('src/fe/fe_interface.C',
      r'case\s+NEDELEC_ONE\s*:\s*switch\s*\(\s*el_t\s*\)\s*\{[^}]*PRISM15',
      [(
          'max_order NEDELEC_ONE',
          ws('case NEDELEC_ONE: switch (el_t) { case TRI6: case TRI7: case QUAD8: case QUAD9: '
             'return 5; case TET10: case TET14: case HEX20: case HEX27:'),
          '''
            case PRISM15:
            case PRISM18:
            case PYRAMID13:
            case PYRAMID14:''',
          'after')])

# ------------------------------------------------------------ fe_nedelec_one.C
NDOFS_BLOCK = '''
      case PRISM15:
      case PRISM18:
        libmesh_assert_less(o, 2);
        return 9;
      case PYRAMID13:
      case PYRAMID14:
        libmesh_assert_less(o, 2);
        return 8;'''

AT_NODE_BLOCK = '''
      case PRISM15:
      case PRISM18:
      {
        libmesh_assert_less(o, 2);

        switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
          return 0;
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
          return o;
        case 15:
        case 16:
        case 17:
          libmesh_assert_equal_to(t, PRISM18);
          return 0;
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n);
        }
      }
      case PYRAMID13:
      case PYRAMID14:
      {
        libmesh_assert_less(o, 2);

        switch (n)
        {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
          return 0;
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
          return o;
        case 13:
          libmesh_assert_equal_to(t, PYRAMID14);
          return 0;
        default:
          libmesh_error_msg("ERROR: Invalid node ID " << n);
        }
      }'''

PER_ELEM_BLOCK = '''
      case PRISM15:
      case PRISM18:
      case PYRAMID13:
      case PYRAMID14:
        libmesh_assert_less(o, 2);
        return 0;'''

apply('src/fe/fe_nedelec_one.C',
      r'PRISM15',
      [
          ('nodal_soln element-type check',
           ws('elem_type != HEX20 && elem_type != HEX27') + r'(?=\s*\))',
           ''' &&
          elem_type != PRISM15 && elem_type != PRISM18 &&
          elem_type != PYRAMID13 && elem_type != PYRAMID14''',
           'after'),
          ('n_dofs',
           ws('return 3 * o * (o + 1) * (o + 1);'),
           NDOFS_BLOCK, 'after'),
          ('n_dofs_at_node',
           ws('case 26: libmesh_assert_equal_to(t, HEX27); return 0; '
              'default: libmesh_error_msg("ERROR: Invalid node ID " << n); } }'),
           AT_NODE_BLOCK, 'after'),
          ('n_dofs_per_elem',
           ws('return 3 * o * (o - 1) * (o - 1);'),
           PER_ELEM_BLOCK, 'after'),
      ])

# ------------------------------------------------------ fe_nedelec_one_shape_3D.C
LOCALS = ''

SHAPE_BLOCK = '''

      case PRISM15:
      case PRISM18:
      {
''' + LOCALS + '''
        // Reference prism: triangle in (xi,eta), zeta in [-1,1].
        libmesh_assert_greater_equal(xi, -10. * TOLERANCE);
        libmesh_assert_greater_equal(eta, -10. * TOLERANCE);
        libmesh_assert_less_equal(xi + eta, 1. + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1. + 10. * TOLERANCE);

        switch (i)
        {
        case 0:
          return sign * RealGradient(-0.5 * (1. - eta) * (1. - zeta), -0.5 * xi * (1. - zeta), 0.0);
        case 1:
          return sign * RealGradient(0.5 * eta * (1. - zeta), -0.5 * xi * (1. - zeta), 0.0);
        case 2:
          return sign * RealGradient(-0.5 * eta * (1. - zeta), -0.5 * (1. - xi) * (1. - zeta), 0.0);
        case 3:
          return sign * RealGradient(0.0, 0.0, -0.5 * (1. - xi - eta));
        case 4:
          return sign * RealGradient(0.0, 0.0, -0.5 * xi);
        case 5:
          return sign * RealGradient(0.0, 0.0, -0.5 * eta);
        case 6:
          return sign * RealGradient(-0.5 * (1. - eta) * (1. + zeta), -0.5 * xi * (1. + zeta), 0.0);
        case 7:
          return sign * RealGradient(0.5 * eta * (1. + zeta), -0.5 * xi * (1. + zeta), 0.0);
        case 8:
          return sign * RealGradient(-0.5 * eta * (1. + zeta), -0.5 * (1. - xi) * (1. + zeta), 0.0);

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }

      case PYRAMID13:
      case PYRAMID14:
      {
''' + LOCALS + '''
        // Reference pyramid: |xi|, |eta| <= 1-zeta, zeta in [0,1].
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1. + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(xi), 1. - zeta + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1. - zeta + 10. * TOLERANCE);

        const Real oz = 1. - zeta;

        // Along the apex the rational basis has direction-dependent
        // limits.  Use the same centerline limit as MFEM, transformed
        // to libMesh coordinates and sign convention.
        if (oz <= 1.e-6)
        {
          switch (i)
          {
          case 0:
          case 1:
          case 2:
          case 3:
            return RealGradient();
          case 4:
            return sign * RealGradient(-0.25, -0.25, -0.25);
          case 5:
            return sign * RealGradient(0.25, -0.25, -0.25);
          case 6:
            return sign * RealGradient(0.25, 0.25, -0.25);
          case 7:
            return sign * RealGradient(-0.25, 0.25, -0.25);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

        const Real x = 0.5 * (1. + xi - zeta);
        const Real y = 0.5 * (1. + eta - zeta);
        const Real ox = 0.5 * (1. - xi - zeta);
        const Real oy = 0.5 * (1. - eta - zeta);
        const Real ozi = 1. / oz;
        const Real q = 1. - 2. * zeta;

        switch (i)
        {
        case 0:
          return sign * RealGradient(-0.5 * oy, 0.0, -x * oy * ozi + 0.5 * oy);
        case 1:
          return sign * RealGradient(0.0, -0.5 * x, -x * y * ozi + 0.5 * x);
        case 2:
          // MFEM's third base edge is oriented 3->2, while libMesh's
          // edge 2 is 2->3.  This is therefore the only base edge
          // which is not negated relative to the MFEM basis.
          return sign * RealGradient(0.5 * y, 0.0, x * y * ozi - 0.5 * y);
        case 3:
          return sign * RealGradient(0.0, -0.5 * ox, -ox * y * ozi + 0.5 * ox);
        case 4:
        {
          const Real vx = oy * zeta * ozi;
          const Real vy = ox * zeta * ozi;
          const Real vz = 1. - x - y + x * y * q * ozi * ozi;
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 5:
        {
          const Real vx = -oy * zeta * ozi;
          const Real vy = x * zeta * ozi;
          const Real vz = x * (1. - y * q * ozi * ozi);
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 6:
        {
          const Real vx = -y * zeta * ozi;
          const Real vy = -x * zeta * ozi;
          const Real vz = x * y * q * ozi * ozi;
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 7:
        {
          const Real vx = y * zeta * ozi;
          const Real vy = -ox * zeta * ozi;
          const Real vz = y * (1. - x * q * ozi * ozi);
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }'''

SHAPE_DERIV_BLOCK = '''

      case PRISM15:
      case PRISM18:
      {
''' + LOCALS + '''
        // Reference prism: triangle in (xi,eta), zeta in [-1,1].
        libmesh_assert_greater_equal(xi, -10. * TOLERANCE);
        libmesh_assert_greater_equal(eta, -10. * TOLERANCE);
        libmesh_assert_less_equal(xi + eta, 1. + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1. + 10. * TOLERANCE);

        switch (j)
        {
          // d()/dxi
        case 0:
        {
          switch (i)
          {
          case 0:
          case 1:
            return sign * RealGradient(0.0, -0.5 * (1. - zeta), 0.0);
          case 2:
            return sign * RealGradient(0.0, 0.5 * (1. - zeta), 0.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 0.5);
          case 4:
            return sign * RealGradient(0.0, 0.0, -0.5);
          case 5:
            return RealGradient();
          case 6:
          case 7:
            return sign * RealGradient(0.0, -0.5 * (1. + zeta), 0.0);
          case 8:
            return sign * RealGradient(0.0, 0.5 * (1. + zeta), 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

          // d()/deta
        case 1:
        {
          switch (i)
          {
          case 0:
          case 1:
            return sign * RealGradient(0.5 * (1. - zeta), 0.0, 0.0);
          case 2:
            return sign * RealGradient(-0.5 * (1. - zeta), 0.0, 0.0);
          case 3:
            return sign * RealGradient(0.0, 0.0, 0.5);
          case 4:
            return RealGradient();
          case 5:
            return sign * RealGradient(0.0, 0.0, -0.5);
          case 6:
          case 7:
            return sign * RealGradient(0.5 * (1. + zeta), 0.0, 0.0);
          case 8:
            return sign * RealGradient(-0.5 * (1. + zeta), 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

          // d()/dzeta
        case 2:
        {
          switch (i)
          {
          case 0:
            return sign * RealGradient(0.5 * (1. - eta), 0.5 * xi, 0.0);
          case 1:
            return sign * RealGradient(-0.5 * eta, 0.5 * xi, 0.0);
          case 2:
            return sign * RealGradient(0.5 * eta, 0.5 * (1. - xi), 0.0);
          case 3:
          case 4:
          case 5:
            return RealGradient();
          case 6:
            return sign * RealGradient(-0.5 * (1. - eta), -0.5 * xi, 0.0);
          case 7:
            return sign * RealGradient(0.5 * eta, -0.5 * xi, 0.0);
          case 8:
            return sign * RealGradient(-0.5 * eta, -0.5 * (1. - xi), 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }

      case PYRAMID13:
      case PYRAMID14:
      {
''' + LOCALS + '''
        // Reference pyramid: |xi|, |eta| <= 1-zeta, zeta in [0,1].
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1. + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(xi), 1. - zeta + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(eta), 1. - zeta + 10. * TOLERANCE);

        const Real eps_apex = 1.e-6;
        const Real oz = std::max(1.0 - zeta, eps_apex); // oz is guaranteed >= 1.e-6
        const Real z = 1.0 - oz;                        // z is guaranteed <= 0.999999
        const Real ozi = 1.0 / oz;                      // Safe division! Never 1/0

        const Real x = 0.5 * (1. + xi - z);
        const Real y = 0.5 * (1. + eta - z);
        const Real ox = 0.5 * (1. - xi - z);
        const Real oy = 0.5 * (1. - eta - z);
        const Real q = 1. - 2. * z;

        const Real dx = (j == 0) ? 0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dy = (j == 1) ? 0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dox = (j == 0) ? -0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real doy = (j == 1) ? -0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dz = (j == 2) ? 1.0 : 0.0;
        const Real dq = (j == 2) ? -2.0 : 0.0;
        const Real dozi = (j == 2) ? ozi * ozi : 0.0;

        switch (i)
        {
        case 0:
        {
          const Real dvx = doy;
          const Real dvy = 0.0;
          const Real dvz = (dx * oy + x * doy) * ozi + x * oy * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 1:
        {
          const Real dvx = 0.0;
          const Real dvy = dx;
          const Real dvz = (dx * y + x * dy) * ozi + x * y * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 2:
        {
          const Real dvx = dy;
          const Real dvy = 0.0;
          const Real dvz = (dx * y + x * dy) * ozi + x * y * dozi;
          return sign * RealGradient(0.5 * dvx, 0.5 * dvy, dvz - 0.5 * (dvx + dvy));
        }
        case 3:
        {
          const Real dvx = 0.0;
          const Real dvy = dox;
          const Real dvz = (dox * y + ox * dy) * ozi + ox * y * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 4:
        {
          const Real dvx = (doy * z + oy * dz) * ozi + oy * z * dozi;
          const Real dvy = (dox * z + ox * dz) * ozi + ox * z * dozi;
          const Real dvz = -dx - dy +
                           (dx * y + x * dy) * q * ozi * ozi +
                           x * y * dq * ozi * ozi +
                           2. * x * y * q * ozi * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 5:
        {
          const Real dvx = -((doy * z + oy * dz) * ozi + oy * z * dozi);
          const Real dvy = (dx * z + x * dz) * ozi + x * z * dozi;
          const Real dvz = dx * (1. - y * q * ozi * ozi) -
                           x * (dy * q * ozi * ozi + y * dq * ozi * ozi +
                                2. * y * q * ozi * dozi);
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 6:
        {
          const Real dvx = -((dy * z + y * dz) * ozi + y * z * dozi);
          const Real dvy = -((dx * z + x * dz) * ozi + x * z * dozi);
          const Real dvz = (dx * y + x * dy) * q * ozi * ozi +
                           x * y * dq * ozi * ozi +
                           2. * x * y * q * ozi * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 7:
        {
          const Real dvx = (dy * z + y * dz) * ozi + y * z * dozi;
          const Real dvy = -((dox * z + ox * dz) * ozi + ox * z * dozi);
          const Real dvz = dy * (1. - x * q * ozi * ozi) -
                           y * (dx * q * ozi * ozi + x * dq * ozi * ozi +
                                2. * x * q * ozi * dozi);
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }'''

SECOND_DERIV_BLOCK = '''      case PRISM15:
      case PRISM18:
      {
        switch (j)
        {
          // d^2()/dxi^2, d^2()/dxideta, d^2()/deta^2, d^2()/dzeta^2
        case 0:
        case 1:
        case 2:
        case 5:
          return RealGradient();

          // d^2()/dxidzeta
        case 3:
        {
          switch (i)
          {
          case 0:
          case 1:
            return sign * RealGradient(0.0, 0.5, 0.0);
          case 2:
            return sign * RealGradient(0.0, -0.5, 0.0);
          case 3:
          case 4:
          case 5:
            return RealGradient();
          case 6:
          case 7:
            return sign * RealGradient(0.0, -0.5, 0.0);
          case 8:
            return sign * RealGradient(0.0, 0.5, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

          // d^2()/detadzeta
        case 4:
        {
          switch (i)
          {
          case 0:
          case 1:
            return sign * RealGradient(-0.5, 0.0, 0.0);
          case 2:
            return sign * RealGradient(0.5, 0.0, 0.0);
          case 3:
          case 4:
          case 5:
            return RealGradient();
          case 6:
          case 7:
            return sign * RealGradient(0.5, 0.0, 0.0);
          case 8:
            return sign * RealGradient(-0.5, 0.0, 0.0);

          default:
            libmesh_error_msg("Invalid i = " << i);
          }
        }

        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }

      case PYRAMID13:
      case PYRAMID14:
        libmesh_error_msg("ERROR: Second derivatives for NEDELEC_ONE pyramid shapes are not implemented.");

'''

apply('src/fe/fe_nedelec_one_shape_3D.C',
      r'PYRAMID13',
      [
          ('shape: prism + pyramid',
           ws('case 5: return sign * RealGradient(0.0, zeta, -eta); '
              'default: libmesh_error_msg("Invalid i = " << i); } }'),
           SHAPE_BLOCK, 'after'),
          ('shape_deriv: prism + pyramid',
           ws('case 3: return sign * RealGradient(-1.0, -1.0, 0.0); '
              'case 4: return sign * RealGradient(1.0, 0.0, 0.0); '
              'default: libmesh_error_msg("Invalid i = " << i); } } '
              'default: libmesh_error_msg("Invalid j = " << j); } }'),
           SHAPE_DERIV_BLOCK, 'after'),
          ('shape_second_deriv: prism + pyramid error',
           r'(?://[^\n]*\n\s*)?' + ws('case TET10: case TET14: return RealGradient();'),
           SECOND_DERIV_BLOCK, 'before'),
      ])

print('libMesh Nedelec prism/pyramid patch applied.')