#!/bin/bash
set -exuo pipefail

SRC_DIR=${SRC_DIR:-$PWD}
if [ ! -d "${SRC_DIR}/src" ] || [ ! -d "${SRC_DIR}/python/apec_m" ]; then
  echo "ERROR: SRC_DIR='${SRC_DIR}' is not the project root."
  echo "Run from the repository root or set SRC_DIR explicitly."
  exit 1
fi
SRC_DIR="$(cd "${SRC_DIR}" && pwd)"
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
