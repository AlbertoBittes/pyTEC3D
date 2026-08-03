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
cp -f "${RECIPE_DIR}/fe_interface_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_interface.C"
cp -f "${RECIPE_DIR}/fe_nedelec_one_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one.C"
cp -f "${RECIPE_DIR}/fe_nedelec_one_shape_3D_new.C" "${LIBMESH_SRC_DIR}/src/fe/fe_nedelec_one_shape_3D.C"

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

LAGRANGE_3D_SRC="${LIBMESH_SRC_DIR}/src/fe/fe_lagrange_shape_3D.C"

python3 -c "
import sys, re

file_path = sys.argv[1]

with open(file_path, 'r') as f:
    content = f.read()

# Target specifically the second-order pyramid block using its unique 13-node assertion anchor
pattern = r'(libmesh_assert_less\s*\(\s*i\s*,\s*13\s*\)\s*;\s*const\s+Real\s+xi\s*=\s*p\(0\)\s*;\s*const\s+Real\s+eta\s*=\s*p\(1\)\s*;\s*)const\s+Real\s+zeta\s*=\s*p\(2\)\s*;'

replacement = r'\1const Real zeta = (p(2) > 1.0 - 1.e-6) ? (1.0 - 1.e-6) : p(2);'

if not re.search(pattern, content):
    print('Error: Could not locate second-order pyramid (PYRAMID13) block in target file.', file=sys.stderr)
    sys.exit(1)

new_content, count = re.subn(pattern, replacement, content)

with open(file_path, 'w') as f:
    f.write(new_content)

print(f'Successfully patched {count} instance(s) in second-order pyramid block of {file_path}')
" "${LAGRANGE_3D_SRC}"

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
