#!/bin/bash
set -exuo pipefail

SRC_DIR=${SRC_DIR:-$PWD}
if [ ! -d "${SRC_DIR}/src" ] || [ ! -d "${SRC_DIR}/python/apec_m" ]; then
  echo "ERROR: SRC_DIR='${SRC_DIR}' is not the project root."
  echo "Run from the repository root or set SRC_DIR explicitly."
  exit 1
fi
SRC_DIR="$(cd "${SRC_DIR}" && pwd)"
RECIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -z "${RECIPE_DIR}" ] || [ ! -f "${RECIPE_DIR}/fe_interface_new.C" ] || [ ! -f "${RECIPE_DIR}/fe_nedelec_one_new.C" ] || [ ! -f "${RECIPE_DIR}/fe_nedelec_one_shape_3D_new.C" ]; then
  echo "ERROR: Could not resolve recipe directory with override files."
  echo "Resolved RECIPE_DIR='${RECIPE_DIR}'"
  exit 1
fi
PREFIX=${PREFIX:-${CONDA_PREFIX:-}}
if [ -z "${PREFIX}" ]; then
  echo "ERROR: PREFIX is not set. Activate a conda environment or export PREFIX explicitly."
  exit 1
fi
export PREFIX

LIBMESH_PREFIX=${LIBMESH_PREFIX:-${PREFIX}}
LIBMESH_SRC_DIR=${SRC_DIR}/libmesh
PYTHON=${PYTHON:-python}
if [ -z "${CPU_COUNT:-}" ]; then
  if command -v nproc >/dev/null 2>&1; then
    CPU_COUNT="$(nproc)"
  elif command -v getconf >/dev/null 2>&1; then
    CPU_COUNT="$(getconf _NPROCESSORS_ONLN)"
  else
    CPU_COUNT=1
  fi
fi
if ! [[ "${CPU_COUNT}" =~ ^[0-9]+$ ]] || [ "${CPU_COUNT}" -le 0 ]; then
  CPU_COUNT=1
fi
APEC_M_LIBMESH_METHOD=${APEC_M_LIBMESH_METHOD:-opt}
APEC_M_LIBMESH_METHOD="$(printf '%s' "${APEC_M_LIBMESH_METHOD}" | tr '[:upper:]' '[:lower:]')"
APEC_M_INSTALL_MODE=${APEC_M_INSTALL_MODE:-copy}

case "${APEC_M_LIBMESH_METHOD}" in
  opt|dbg|devel) ;;
  *)
    echo "ERROR: Invalid APEC_M_LIBMESH_METHOD='${APEC_M_LIBMESH_METHOD}'. Allowed: opt, dbg, devel."
    exit 1
    ;;
esac
export APEC_M_LIBMESH_METHOD

case "${APEC_M_INSTALL_MODE}" in
  copy|editable) ;;
  *)
    echo "ERROR: Invalid APEC_M_INSTALL_MODE='${APEC_M_INSTALL_MODE}'. Allowed: copy, editable."
    exit 1
    ;;
esac

if [ "${APEC_M_INSTALL_MODE}" = "editable" ] && [ -n "${CONDA_BUILD_STATE:-}" ]; then
  echo "ERROR: APEC_M_INSTALL_MODE=editable is for local debugging only and is not supported during conda-build."
  exit 1
fi

export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export FC=${FC:-mpifort}
export F77=${F77:-mpifort}
export F90=${F90:-mpifort}
export PETSC_DIR=${PETSC_DIR:-${PREFIX}}
export METIS_DIR=${METIS_DIR:-${PREFIX}}
export PARMETIS_DIR=${PARMETIS_DIR:-${PREFIX}}
export MPI_DIR=${MPI_DIR:-${PREFIX}}
export LIBMESH_DIR=${LIBMESH_PREFIX}

VTK_CMAKE_DIR=${VTK_CMAKE_DIR:-}
if [ -z "${VTK_CMAKE_DIR}" ]; then
  for candidate in "${PREFIX}"/lib/cmake/vtk-* "${PREFIX}"/lib/cmake/vtk; do
    if [ -d "${candidate}" ]; then
      VTK_CMAKE_DIR="${candidate}"
      break
    fi
  done
fi
if [ -z "${VTK_CMAKE_DIR}" ]; then
  echo "ERROR: Could not find VTK CMake config directory under ${PREFIX}/lib/cmake"
  exit 1
fi

VTK_VERSION="$("${PYTHON}" - <<'PY'
import sys

try:
    import vtk
except Exception as exc:
    print(f"ERROR: failed to import vtk: {exc}", file=sys.stderr)
    sys.exit(2)

version = getattr(vtk, "__version__", "0.0.0")
parts = version.split(".")
major = int(parts[0]) if len(parts) > 0 and parts[0].isdigit() else 0
minor = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
if (major, minor) < (9, 5):
    print(f"ERROR: vtk>=9.5 is required, found vtk=={version}", file=sys.stderr)
    sys.exit(3)

print(version)
PY
)"
echo "Detected vtk version: ${VTK_VERSION}"

CMAKE_BUILD_TYPE=Release
APEC_M_METHOD_CFLAGS_DEFAULT="-O3 -g0"
if [ "${APEC_M_LIBMESH_METHOD}" = "dbg" ]; then
  CMAKE_BUILD_TYPE=Debug
  APEC_M_METHOD_CFLAGS_DEFAULT="-O0 -g3"
elif [ "${APEC_M_LIBMESH_METHOD}" = "devel" ]; then
  CMAKE_BUILD_TYPE=RelWithDebInfo
  APEC_M_METHOD_CFLAGS_DEFAULT="-O2 -g"
fi

APEC_M_METHOD_CFLAGS=${APEC_M_METHOD_CFLAGS:-${APEC_M_METHOD_CFLAGS_DEFAULT}}
APEC_M_METHOD_CXXFLAGS=${APEC_M_METHOD_CXXFLAGS:-${APEC_M_METHOD_CFLAGS}}
APEC_M_EXTRA_CFLAGS=${APEC_M_EXTRA_CFLAGS:-}
APEC_M_EXTRA_CXXFLAGS=${APEC_M_EXTRA_CXXFLAGS:-}
APEC_M_EXTRA_CPPFLAGS=${APEC_M_EXTRA_CPPFLAGS:-}
APEC_M_EXTRA_LDFLAGS=${APEC_M_EXTRA_LDFLAGS:-}

export CFLAGS="${CFLAGS:-} -I${PREFIX}/include ${APEC_M_METHOD_CFLAGS} ${APEC_M_EXTRA_CFLAGS}"
export CPPFLAGS="${CPPFLAGS:-} -I${PREFIX}/include ${APEC_M_EXTRA_CPPFLAGS}"
export CXXFLAGS="${CXXFLAGS:-} -I${PREFIX}/include ${APEC_M_METHOD_CXXFLAGS} ${APEC_M_EXTRA_CXXFLAGS}"
export LDFLAGS="${LDFLAGS:-} -L${PREFIX}/lib ${APEC_M_EXTRA_LDFLAGS}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LIBMESH_LIB="${LIBMESH_PREFIX}/lib/libmesh_${APEC_M_LIBMESH_METHOD}.so"
export TIMPI_LIB="${LIBMESH_PREFIX}/lib/libtimpi_${APEC_M_LIBMESH_METHOD}.so"

echo "Cleaning previous local libMesh/APEC-M build artifacts"
rm -rf \
  "${LIBMESH_PREFIX}/include/libmesh" \
  "${LIBMESH_PREFIX}/share/libmesh" \
  "${LIBMESH_PREFIX}/lib/cmake/libmesh" \
  "${SRC_DIR}/libmesh" \
  "${SRC_DIR}/build/apec_m"
rm -f \
  "${LIBMESH_PREFIX}/bin/libmesh-config" \
  "${LIBMESH_PREFIX}/lib"/libmesh_*.a \
  "${LIBMESH_PREFIX}/lib"/libmesh_*.la \
  "${LIBMESH_PREFIX}/lib"/libmesh_*.so \
  "${LIBMESH_PREFIX}/lib"/libmesh_*.so.* \
  "${LIBMESH_PREFIX}/lib"/libtimpi_*.a \
  "${LIBMESH_PREFIX}/lib"/libtimpi_*.la \
  "${LIBMESH_PREFIX}/lib"/libtimpi_*.so \
  "${LIBMESH_PREFIX}/lib"/libtimpi_*.so.* \
  "${LIBMESH_PREFIX}/lib/pkgconfig"/libmesh*.pc \
  "${LIBMESH_PREFIX}/lib/pkgconfig"/timpi*.pc \
  "${SRC_DIR}/python/apec_m"/apec_m_solver*.so \
  "${SRC_DIR}/python/apec_m"/apec_m_solver*.pyd \
  "${SRC_DIR}/python/apec_m"/apec_m_solver*.dylib

if [ "${APEC_M_LIBMESH_METHOD}" = "dbg" ]; then
  # Conda toolchain activation may inject -DNDEBUG/-D_FORTIFY_SOURCE; sanitize for dbg builds.
  for flag_var in CFLAGS CPPFLAGS CXXFLAGS; do
    current_val="$(eval "printf '%s' \"\${${flag_var}}\"")"
    current_val="$(printf '%s' "${current_val}" | sed -E 's/(^|[[:space:]])-DNDEBUG(=[^[:space:]]*)?([[:space:]]|$)/ /g')"
    current_val="$(printf '%s' "${current_val}" | sed -E 's/(^|[[:space:]])-D_FORTIFY_SOURCE(=[^[:space:]]*)?([[:space:]]|$)/ /g')"
    current_val="$(printf '%s' "${current_val}" | sed -E 's/(^|[[:space:]])-U_FORTIFY_SOURCE([[:space:]]|$)/ /g')"
    current_val="$(echo "${current_val}" | xargs)"
    eval "export ${flag_var}=\"${current_val} -UNDEBUG -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0\""
  done
fi

# # # 1. Clone and build libMesh with the requested options
git clone --depth 1 --branch master https://github.com/libMesh/libmesh.git "${LIBMESH_SRC_DIR}"
pushd "${LIBMESH_SRC_DIR}"
git submodule update --init --recursive

# # Replace selected libMesh FE sources with recipe-maintained overrides.
# cp -f "${RECIPE_DIR}/fe_interface_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_interface.C"
# cp -f "${RECIPE_DIR}/fe_nedelec_one_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one.C"
# cp -f "${RECIPE_DIR}/fe_nedelec_one_shape_3D_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one_shape_3D.C"

#!/usr/bin/env bash
set -e

PYRAMID13_SRC="${LIBMESH_SRC_DIR}/src/geom/cell_pyramid13.C"

python3 -c "
import sys, re

file_path = sys.argv[1]

with open(file_path, 'r') as f:
    content = f.read()

# Define the updated connectivity function
new_func = '''void Pyramid13::connectivity(const unsigned int sc,
                           const IOPackage iop,
                           std::vector< dof_id_type > & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        // TODO
        libmesh_not_implemented();
      }

    case VTK:
      {
        switch (sc)
          {
          case 0:
            // libMesh and VTK node orderings match perfectly for Pyramid13.
            conn.resize(13);
            conn[0]  = this->node_id(0);
            conn[1]  = this->node_id(1);
            conn[2]  = this->node_id(2);
            conn[3]  = this->node_id(3);
            conn[4]  = this->node_id(4);
            conn[5]  = this->node_id(5);
            conn[6]  = this->node_id(6);
            conn[7]  = this->node_id(7);
            conn[8]  = this->node_id(8);
            conn[9]  = this->node_id(9);
            conn[10] = this->node_id(10);
            conn[11] = this->node_id(11);
            conn[12] = this->node_id(12);
            return;

          default:
            libmesh_error_msg(\"Unknown sc = \" << sc);
          }
      }

    default:
      libmesh_error_msg(\"Unsupported IO package \" << iop);
    }
}'''

# Regex to find and replace the void Pyramid13::connectivity method definition
pattern = r'void\s+Pyramid13::connectivity\s*\([\s\S]*?^\}'

if not re.search(pattern, content, re.MULTILINE):
    print('Error: Could not locate Pyramid13::connectivity function in target file.', file=sys.stderr)
    sys.exit(1)

new_content = re.sub(pattern, new_func, content, flags=re.MULTILINE)

with open(file_path, 'w') as f:
    f.write(new_content)

print(f'Successfully patched Pyramid13::connectivity in {file_path}')
" "${PYRAMID13_SRC}"

#!/usr/bin/env bash
set -e

NEDELEC_SHAPE_3D_SRC="${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one_shape_3D.C"

python3 - "$NEDELEC_SHAPE_3D_SRC" << 'EOF'
import sys
import re

file_path = sys.argv[1]

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Replace non-breaking spaces (\xa0) with standard spaces
content = content.replace('\xa0', ' ')

prism_pyramid_shape = '''      case PRISM15:
      case PRISM18:
      {
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
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1. + 10. * TOLERANCE);

        const Real eps_apex = 1.e-6;
        const Real clamped_zeta = std::min(zeta, 1.0 - eps_apex);
        const Real oz = 1.0 - clamped_zeta;

        const Real clamped_xi  = std::max(-oz, std::min(xi, oz));
        const Real clamped_eta = std::max(-oz, std::min(eta, oz));

        const Real x   = 0.5 * (1. + clamped_xi - clamped_zeta);
        const Real y   = 0.5 * (1. + clamped_eta - clamped_zeta);
        const Real ox  = 0.5 * (1. - clamped_xi - clamped_zeta);
        const Real oy  = 0.5 * (1. - clamped_eta - clamped_zeta);
        const Real ozi = 1. / oz;
        const Real q   = 1. - 2. * clamped_zeta;

        switch (i)
        {
        case 0:
          return sign * RealGradient(-0.5 * oy, 0.0, -x * oy * ozi + 0.5 * oy);
        case 1:
          return sign * RealGradient(0.0, -0.5 * x, -x * y * ozi + 0.5 * x);
        case 2:
          return sign * RealGradient(0.5 * y, 0.0, x * y * ozi - 0.5 * y);
        case 3:
          return sign * RealGradient(0.0, -0.5 * ox, -ox * y * ozi + 0.5 * ox);
        case 4:
        {
          const Real vx = oy * clamped_zeta * ozi;
          const Real vy = ox * clamped_zeta * ozi;
          const Real vz = 1. - x - y + x * y * q * ozi * ozi;
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 5:
        {
          const Real vx = -oy * clamped_zeta * ozi;
          const Real vy = x * clamped_zeta * ozi;
          const Real vz = x * (1. - y * q * ozi * ozi);
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 6:
        {
          const Real vx = -y * clamped_zeta * ozi;
          const Real vy = -x * clamped_zeta * ozi;
          const Real vz = x * y * q * ozi * ozi;
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }
        case 7:
        {
          const Real vx = y * clamped_zeta * ozi;
          const Real vy = -ox * clamped_zeta * ozi;
          const Real vz = y * (1. - x * q * ozi * ozi);
          return sign * RealGradient(-0.5 * vx, -0.5 * vy, -vz + 0.5 * (vx + vy));
        }

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }'''

prism_pyramid_deriv = '''      case PRISM15:
      case PRISM18:
      {
        libmesh_assert_greater_equal(xi, -10. * TOLERANCE);
        libmesh_assert_greater_equal(eta, -10. * TOLERANCE);
        libmesh_assert_less_equal(xi + eta, 1. + 10. * TOLERANCE);
        libmesh_assert_less_equal(std::fabs(zeta), 1. + 10. * TOLERANCE);

        switch (j)
        {
        case 0:
        {
          switch (i)
          {
          case 0: case 1: return sign * RealGradient(0.0, -0.5 * (1. - zeta), 0.0);
          case 2: return sign * RealGradient(0.0, 0.5 * (1. - zeta), 0.0);
          case 3: return sign * RealGradient(0.0, 0.0, 0.5);
          case 4: return sign * RealGradient(0.0, 0.0, -0.5);
          case 5: return RealGradient();
          case 6: case 7: return sign * RealGradient(0.0, -0.5 * (1. + zeta), 0.0);
          case 8: return sign * RealGradient(0.0, 0.5 * (1. + zeta), 0.0);
          default: libmesh_error_msg("Invalid i = " << i);
          }
        }
        case 1:
        {
          switch (i)
          {
          case 0: case 1: return sign * RealGradient(0.5 * (1. - zeta), 0.0, 0.0);
          case 2: return sign * RealGradient(-0.5 * (1. - zeta), 0.0, 0.0);
          case 3: return sign * RealGradient(0.0, 0.0, 0.5);
          case 4: return RealGradient();
          case 5: return sign * RealGradient(0.0, 0.0, -0.5);
          case 6: case 7: return sign * RealGradient(0.5 * (1. + zeta), 0.0, 0.0);
          case 8: return sign * RealGradient(-0.5 * (1. + zeta), 0.0, 0.0);
          default: libmesh_error_msg("Invalid i = " << i);
          }
        }
        case 2:
        {
          switch (i)
          {
          case 0: return sign * RealGradient(0.5 * (1. - eta), 0.5 * xi, 0.0);
          case 1: return sign * RealGradient(-0.5 * eta, 0.5 * xi, 0.0);
          case 2: return sign * RealGradient(0.5 * eta, 0.5 * (1. - xi), 0.0);
          case 3: case 4: case 5: return RealGradient();
          case 6: return sign * RealGradient(-0.5 * (1. - eta), -0.5 * xi, 0.0);
          case 7: return sign * RealGradient(0.5 * eta, -0.5 * xi, 0.0);
          case 8: return sign * RealGradient(-0.5 * eta, -0.5 * (1. - xi), 0.0);
          default: libmesh_error_msg("Invalid i = " << i);
          }
        }
        default:
          libmesh_error_msg("Invalid j = " << j);
        }
      }

      case PYRAMID13:
      case PYRAMID14:
      {
        libmesh_assert_greater_equal(zeta, -10. * TOLERANCE);
        libmesh_assert_less_equal(zeta, 1. + 10. * TOLERANCE);

        const Real eps_apex = 1.e-6;
        const Real clamped_zeta = std::min(zeta, 1.0 - eps_apex);
        const Real oz = 1.0 - clamped_zeta;

        const Real clamped_xi  = std::max(-oz, std::min(xi, oz));
        const Real clamped_eta = std::max(-oz, std::min(eta, oz));

        const Real x   = 0.5 * (1. + clamped_xi - clamped_zeta);
        const Real y   = 0.5 * (1. + clamped_eta - clamped_zeta);
        const Real ox  = 0.5 * (1. - clamped_xi - clamped_zeta);
        const Real oy  = 0.5 * (1. - clamped_eta - clamped_zeta);
        const Real ozi = 1. / oz;
        const Real q   = 1. - 2. * clamped_zeta;

        const Real dx   = (j == 0) ? 0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dy   = (j == 1) ? 0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dox  = (j == 0) ? -0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real doy  = (j == 1) ? -0.5 : ((j == 2) ? -0.5 : 0.0);
        const Real dz   = (j == 2) ? 1.0 : 0.0;
        const Real dq   = (j == 2) ? -2.0 : 0.0;
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
          const Real dvx = (doy * clamped_zeta + oy * dz) * ozi + oy * clamped_zeta * dozi;
          const Real dvy = (dox * clamped_zeta + ox * dz) * ozi + ox * clamped_zeta * dozi;
          const Real dvz = -dx - dy +
                           (dx * y + x * dy) * q * ozi * ozi +
                           x * y * dq * ozi * ozi +
                           2. * x * y * q * ozi * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 5:
        {
          const Real dvx = -((doy * clamped_zeta + oy * dz) * ozi + oy * clamped_zeta * dozi);
          const Real dvy = (dx * clamped_zeta + x * dz) * ozi + x * clamped_zeta * dozi;
          const Real dvz = dx * (1. - y * q * ozi * ozi) -
                           x * (dy * q * ozi * ozi + y * dq * ozi * ozi +
                                2. * y * q * ozi * dozi);
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 6:
        {
          const Real dvx = -((dy * clamped_zeta + y * dz) * ozi + y * clamped_zeta * dozi);
          const Real dvy = -((dx * clamped_zeta + x * dz) * ozi + x * clamped_zeta * dozi);
          const Real dvz = (dx * y + x * dy) * q * ozi * ozi +
                           x * y * dq * ozi * ozi +
                           2. * x * y * q * ozi * dozi;
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }
        case 7:
        {
          const Real dvx = (dy * clamped_zeta + y * dz) * ozi + y * clamped_zeta * dozi;
          const Real dvy = -((dox * clamped_zeta + ox * dz) * ozi + ox * clamped_zeta * dozi);
          const Real dvz = dy * (1. - x * q * ozi * ozi) -
                           y * (dx * q * ozi * ozi + x * dq * ozi * ozi +
                                2. * x * q * ozi * dozi);
          return sign * RealGradient(-0.5 * dvx, -0.5 * dvy, -dvz + 0.5 * (dvx + dvy));
        }

        default:
          libmesh_error_msg("Invalid i = " << i);
        }
      }'''

def patch_section(section_text, new_code):
    # 1. If case PRISM15 already exists, replace from PRISM15 up to default:
    prism_match = re.search(r'case\s+PRISM15:', section_text)
    if prism_match:
        pattern = re.compile(r'case\s+PRISM15:[\s\S]*?(?=\n\s*default:)', re.MULTILINE)
        return pattern.sub(new_code + '\n', section_text, count=1)
    
    # 2. Otherwise, insert right before the default: case in the switch statement
    default_match = re.search(r'(\n\s*default:)', section_text)
    if default_match:
        return section_text[:default_match.start()] + '\n' + new_code + '\n' + section_text[default_match.start():]
    
    return None

deriv_idx = content.find("shape_deriv(")
if deriv_idx == -1:
    print(f"Error: Could not locate shape_deriv() in {file_path}", file=sys.stderr)
    sys.exit(1)

part_shape = content[:deriv_idx]
part_deriv = content[deriv_idx:]

new_part_shape = patch_section(part_shape, prism_pyramid_shape)
new_part_deriv = patch_section(part_deriv, prism_pyramid_deriv)

if new_part_shape is None or new_part_deriv is None:
    print("Error: Could not find insertion point ('default:' block) in target file.", file=sys.stderr)
    sys.exit(1)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(new_part_shape + new_part_deriv)

print(f"Successfully patched PRISM and PYRAMID shape functions in {file_path}")
EOF

#!/usr/bin/env bash
set -e

NEDELEC_ONE_SRC="${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one.C"

python3 - "$NEDELEC_ONE_SRC" << 'EOF'
import sys
import re

file_path = sys.argv[1]

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Clean up non-breaking spaces (\xa0)
content = content.replace('\xa0', ' ')

# 1. Patch valid element check inside nedelec_one_nodal_soln
if "PRISM15" not in content:
    nodal_pattern = r'if\s*\(\s*elem_type\s*!=\s*TRI6[\s\S]*?elem_type\s*!=\s*HEX27\s*\)'
    nodal_replacement = '''if (elem_type != TRI6 && elem_type != TRI7 &&
          elem_type != QUAD8 && elem_type != QUAD9 &&
          elem_type != TET10 && elem_type != TET14 &&
          elem_type != HEX20 && elem_type != HEX27 &&
          elem_type != PRISM15 && elem_type != PRISM18 &&
          elem_type != PYRAMID13 && elem_type != PYRAMID14)'''
    
    content, count = re.subn(nodal_pattern, nodal_replacement, content, count=1)
    if count == 0:
        # Fallback regex matching standard element checks
        content, count = re.subn(
            r'(elem_type\s*!=\s*HEX27\s*\))',
            r'elem_type != HEX27 &&\n          elem_type != PRISM15 && elem_type != PRISM18 &&\n          elem_type != PYRAMID13 && elem_type != PYRAMID14)',
            content,
            count=1
        )
        if count == 0:
            print("Error: Could not locate nodal check block in nedelec_one_nodal_soln.", file=sys.stderr)
            sys.exit(1)

# Helper function to insert case block cleanly inside a specific function switch statement
def insert_case_block(func_name, code_to_insert, content_str):
    func_match = re.search(rf'{func_name}\s*\([^)]*\)\s*\{{\s*[\s\S]*?switch\s*\([^\)]+\)\s*\{{', content_str)
    if not func_match:
        print(f"Error: Could not find function signature or switch statement for {func_name}.", file=sys.stderr)
        sys.exit(1)
    
    pos = func_match.end()
    target_match = re.search(r'(\n\s*case\s+INVALID_ELEM:|\n\s*default:)', content_str[pos:])
    if not target_match:
        print(f"Error: Could not find insertion point (INVALID_ELEM or default) in {func_name}.", file=sys.stderr)
        sys.exit(1)
    
    target_pos = pos + target_match.start()
    
    # Check idempotency
    if "PRISM15" in content_str[pos:target_pos]:
        return content_str
    
    return content_str[:target_pos] + "\n" + code_to_insert + content_str[target_pos:]

# 2. Patch nedelec_one_n_dofs
n_dofs_code = '''      case PRISM15:
      case PRISM18:
        libmesh_assert_less(o, 2);
        return 9;
      case PYRAMID13:
      case PYRAMID14:
        libmesh_assert_less(o, 2);
        return 8;'''
content = insert_case_block('nedelec_one_n_dofs', n_dofs_code, content)

# 3. Patch nedelec_one_n_dofs_at_node
n_dofs_node_code = '''      case PRISM15:
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
content = insert_case_block('nedelec_one_n_dofs_at_node', n_dofs_node_code, content)

# 4. Patch nedelec_one_n_dofs_per_elem
n_dofs_per_elem_code = '''      case PRISM15:
      case PRISM18:
      case PYRAMID13:
      case PYRAMID14:
        libmesh_assert_less(o, 2);
        return 0;'''
content = insert_case_block('nedelec_one_n_dofs_per_elem', n_dofs_per_elem_code, content)

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)

print(f"Successfully patched PRISM and PYRAMID FE info in {file_path}")
EOF

#!/usr/bin/env bash
set -e

FE_INTERFACE_SRC="${LIBMESH_SRC_DIR}/src/fe/fe_interface.C"

python3 - "$FE_INTERFACE_SRC" << 'EOF'
import sys
import re

file_path = sys.argv[1]

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Clean up non-breaking spaces (\xa0)
content = content.replace('\xa0', ' ')

# Match the switch statement inside case NEDELEC_ONE
pattern = r'(case\s+NEDELEC_ONE:\s*\n?\s*switch\s*\(\s*el_t\s*\)\s*\{)[\s\S]*?(default:\s*\n?\s*return\s+0;\s*\n?\s*\}\s*\n?\s*break;)'

new_cases = '''
            case TRI6:
            case TRI7:
            case QUAD8:
            case QUAD9:
                return 5;
            case TET10:
            case TET14:
            case HEX20:
            case HEX27:
            case PRISM15:
            case PRISM18:
            case PYRAMID13:
            case PYRAMID14:
                return 1;'''

match = re.search(pattern, content)
if not match:
    print('Error: Could not locate NEDELEC_ONE block in FEInterface::max_order.', file=sys.stderr)
    sys.exit(1)

# Reconstruct block keeping original header/footer formatting intact
new_block = match.group(1) + new_cases + '\n            ' + match.group(2)
new_content = content[:match.start()] + new_block + content[match.end():]

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(new_content)

print(f'Successfully patched NEDELEC_ONE order cases in {file_path}')
EOF

mkdir build
pushd build
../configure \
    --prefix="${LIBMESH_PREFIX}" \
    --with-methods="${APEC_M_LIBMESH_METHOD}" \
    --with-petsc="${PETSC_DIR}" \
    --with-metis="${METIS_DIR}" \
    --with-parmetis="${PARMETIS_DIR}" \
    --with-mpi="${MPI_DIR}" \
    --disable-vtk \
    --disable-exodus \
    --disable-libHilbert \
    --disable-tetgen \
    --disable-triangle \
    --disable-netgen \
    --disable-poly2tri \
    --disable-nlopt \
    --disable-gmv \
    --disable-netcdf \
    --enable-distmesh \
    --disable-perflog \
    --disable-eigen \
    --disable-fparser \
    --disable-metaphysicl \
    --enable-periodic \
    --disable-nanoflann \
    --disable-qhull \
    --disable-second \
    --disable-capnproto \
    --disable-amr
make -j"${CPU_COUNT}"
make install
popd
popd

# 2. Build the apec_m_solver pybind module with CMake
mkdir -p "${SRC_DIR}/build/apec_m"
pushd "${SRC_DIR}/build/apec_m"
cmake "${SRC_DIR}/src" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DLIBMESH_DIR="${LIBMESH_PREFIX}" \
    -DLIBMESH_METHOD="${APEC_M_LIBMESH_METHOD}" \
    -DLIBMESH_LIB:FILEPATH="${LIBMESH_LIB}" \
    -DTIMPI_LIB:FILEPATH="${TIMPI_LIB}" \
    -DVTK_DIR="${VTK_CMAKE_DIR}" \
    -DAPEC_M_BUILD_PYBIND_MODULE=ON
ninja -j"${CPU_COUNT}"

shopt -s nullglob
solver_modules=(
  "${SRC_DIR}/python/apec_m/apec_m_solver"*.so
  "${SRC_DIR}/python/apec_m/apec_m_solver"*.pyd
  "${SRC_DIR}/python/apec_m/apec_m_solver"*.dylib
)
if [ "${#solver_modules[@]}" -eq 0 ]; then
  echo "ERROR: Expected a compiled apec_m_solver extension under ${SRC_DIR}/python/apec_m but none was produced."
  exit 1
fi
shopt -u nullglob
popd

# 3. Install the Python module into the conda prefix
pushd "${SRC_DIR}"
if [ "${APEC_M_INSTALL_MODE}" = "editable" ]; then
  "${PYTHON}" -m pip install --no-build-isolation --no-deps --config-settings editable_mode=compat -e .
else
  "${PYTHON}" -m pip install --no-build-isolation --no-deps --prefix "${PREFIX}" .
fi
popd

# 4. Install Python runtime deps not reliably available via conda.
# In conda-build, pandas/openpyxl come from meta.yaml host deps.
if [ -n "${CONDA_BUILD_STATE:-}" ]; then
  "${PYTHON}" -m pip install --prefix "${PREFIX}" "ansys-mapdl-reader==0.55.2"
else
  # Local build convenience path: ensure apec_m Python deps are present.
  "${PYTHON}" -m pip install --prefix "${PREFIX}" --upgrade \
    pandas \
    openpyxl \
    "ansys-mapdl-reader==0.55.2"
fi

# 5. Vendor openpyxl into repository-local python/apec_m/local_packages 
"${PYTHON}" -m pip install --no-deps --upgrade openpyxl et_xmlfile --target "${SRC_DIR}/python/apec_m/local_packages"
