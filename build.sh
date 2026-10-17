#!/bin/bash
# pyTEC3D build: libMesh -> solver module. Run `pip install .` or
# `pip install -e .` afterwards.
#
# Env vars:
#   pyTEC3D_LIBMESH_METHOD=opt|dbg|devel   libMesh method (default: opt)
#   FORCE_LIBMESH=1                        wipe and rebuild libMesh (default: reuse if present)
#   PREFIX, PYTHON, CPU_COUNT, DEPS_DIR, PATCH_DIR
#   pyTEC3D_METHOD_CFLAGS / _CXXFLAGS, pyTEC3D_EXTRA_{CFLAGS,CXXFLAGS,CPPFLAGS,LDFLAGS}
#
# Layout (derived from the repo root):
#   .deps/src/libmesh            libMesh checkout (patched once)
#   .deps/libmesh-<method>/      libMesh install prefix (baked into the module's RPATH)
#   build/pytec3d-<method>/      CMake build tree (wiped on every run)
#   src/pytec3d_solver/*.so      compiled module, installed there by `cmake --install`
#
# Incremental rebuild after editing the C++ sources:
#   cmake --build build/pytec3d-<method> && cmake --install build/pytec3d-<method>
set -euo pipefail

# --- Options (env only) ----------------------------------------------------------
pyTEC3D_LIBMESH_METHOD="$(printf '%s' "${pyTEC3D_LIBMESH_METHOD:-opt}" | tr '[:upper:]' '[:lower:]')"
case "${pyTEC3D_LIBMESH_METHOD}" in opt|dbg|devel) ;;
  *) echo "ERROR: invalid pyTEC3D_LIBMESH_METHOD='${pyTEC3D_LIBMESH_METHOD}'. Allowed: opt, dbg, devel."; exit 1 ;;
esac
export pyTEC3D_LIBMESH_METHOD
FORCE_LIBMESH="${FORCE_LIBMESH:-0}"

# --- Paths ---------------------------------------------------------------------
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -d "${SRC_DIR}/src/pytec3d" ]; then
  echo "ERROR: '${SRC_DIR}' is not the project root (missing src/pytec3d)."
  exit 1
fi

PREFIX="${PREFIX:-${CONDA_PREFIX:-}}"
if [ -z "${PREFIX}" ]; then
  echo "ERROR: Could not resolve PREFIX. Activate the conda env or export PREFIX."
  exit 1
fi
export PREFIX
PYTHON="${PYTHON:-python}"

DEPS_DIR="${DEPS_DIR:-${SRC_DIR}/.deps}"
PATCH_DIR="${PATCH_DIR:-${SRC_DIR}/libmesh_patches}"
LIBMESH_SRC_DIR="${DEPS_DIR}/src/libmesh"
LIBMESH_PREFIX="${DEPS_DIR}/libmesh-${pyTEC3D_LIBMESH_METHOD}"
LIBMESH_LIB="${LIBMESH_PREFIX}/lib/libmesh_${pyTEC3D_LIBMESH_METHOD}.so"
SOLVER_BUILD_DIR="${SRC_DIR}/build/pytec3d-${pyTEC3D_LIBMESH_METHOD}"
PYTEC3D_PYTHON_PACKAGE_ROOT="${SRC_DIR}/src"
SOLVER_PKG_DIR="${PYTEC3D_PYTHON_PACKAGE_ROOT}/pytec3d_solver"

if [ -z "${CPU_COUNT:-}" ]; then
  CPU_COUNT="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi
if ! [[ "${CPU_COUNT}" =~ ^[0-9]+$ ]] || [ "${CPU_COUNT}" -le 0 ]; then CPU_COUNT=1; fi

# --- Toolchain / flags -----------------------------------------------------------
export CC="${CC:-mpicc}" CXX="${CXX:-mpicxx}"
export FC="${FC:-mpifort}" F77="${F77:-mpifort}" F90="${F90:-mpifort}"
export PETSC_DIR="${PETSC_DIR:-${PREFIX}}"
export METIS_DIR="${METIS_DIR:-${PREFIX}}"
export PARMETIS_DIR="${PARMETIS_DIR:-${PREFIX}}"
export MPI_DIR="${MPI_DIR:-${PREFIX}}"

CMAKE_BUILD_TYPE=Release
METHOD_FLAGS_DEFAULT="-O3 -g0"
if [ "${pyTEC3D_LIBMESH_METHOD}" = "dbg" ]; then
  CMAKE_BUILD_TYPE=Debug;          METHOD_FLAGS_DEFAULT="-O0 -g3"
elif [ "${pyTEC3D_LIBMESH_METHOD}" = "devel" ]; then
  CMAKE_BUILD_TYPE=RelWithDebInfo; METHOD_FLAGS_DEFAULT="-O2 -g"
fi
pyTEC3D_METHOD_CFLAGS="${pyTEC3D_METHOD_CFLAGS:-${METHOD_FLAGS_DEFAULT}}"
pyTEC3D_METHOD_CXXFLAGS="${pyTEC3D_METHOD_CXXFLAGS:-${pyTEC3D_METHOD_CFLAGS}}"

export CFLAGS="${CFLAGS:-} -I${PREFIX}/include ${pyTEC3D_METHOD_CFLAGS} ${pyTEC3D_EXTRA_CFLAGS:-}"
export CPPFLAGS="${CPPFLAGS:-} -I${PREFIX}/include ${pyTEC3D_EXTRA_CPPFLAGS:-}"
export CXXFLAGS="${CXXFLAGS:-} -I${PREFIX}/include ${pyTEC3D_METHOD_CXXFLAGS} ${pyTEC3D_EXTRA_CXXFLAGS:-}"
export LDFLAGS="${LDFLAGS:-} -L${PREFIX}/lib ${pyTEC3D_EXTRA_LDFLAGS:-}"
# Build-time convenience only (configure runs small test binaries). The final
# module must not depend on it: it is verified below with LD_LIBRARY_PATH unset.
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

# Toolchains may inject -DNDEBUG/-D_FORTIFY_SOURCE; strip them for dbg builds.
if [ "${pyTEC3D_LIBMESH_METHOD}" = "dbg" ]; then
  for flag_var in CFLAGS CPPFLAGS CXXFLAGS; do
    val="${!flag_var}"
    val="$(printf '%s' "${val}" | sed -E 's/(^|[[:space:]])-DNDEBUG(=[^[:space:]]*)?([[:space:]]|$)/ /g')"
    val="$(printf '%s' "${val}" | sed -E 's/(^|[[:space:]])-D_FORTIFY_SOURCE(=[^[:space:]]*)?([[:space:]]|$)/ /g')"
    val="$(printf '%s' "${val}" | sed -E 's/(^|[[:space:]])-U_FORTIFY_SOURCE([[:space:]]|$)/ /g')"
    val="$(echo "${val}" | xargs)"
    export "${flag_var}=${val} -UNDEBUG -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0"
  done
fi

# =================================================================================
# Stage 1: libMesh
# =================================================================================
build_libmesh() {
  if [ -f "${LIBMESH_LIB}" ] && [ "${FORCE_LIBMESH}" != "1" ]; then
    echo "==> libMesh (${pyTEC3D_LIBMESH_METHOD}) already installed in ${LIBMESH_PREFIX}; skipping (FORCE_LIBMESH=1 to rebuild)."
    return
  fi
  echo "==> Building libMesh (${pyTEC3D_LIBMESH_METHOD})"
  [ -d "${PATCH_DIR}" ] || { echo "ERROR: patch directory not found: '${PATCH_DIR}'"; exit 1; }

  rm -rf "${LIBMESH_PREFIX}"
  mkdir -p "${LIBMESH_PREFIX}" "$(dirname "${LIBMESH_SRC_DIR}")"

  # Clone + patch only when the checkout is missing or a rebuild is forced.
  if [ ! -d "${LIBMESH_SRC_DIR}/.git" ] || [ "${FORCE_LIBMESH}" = "1" ]; then
    rm -rf "${LIBMESH_SRC_DIR}"
    git clone --depth 1 --branch master https://github.com/libMesh/libmesh.git "${LIBMESH_SRC_DIR}"
    git -C "${LIBMESH_SRC_DIR}" submodule update --init --recursive

    apply_patch() {  # <script> <target>
      [ -f "${PATCH_DIR}/$1" ] || { echo "ERROR: missing patch script '${PATCH_DIR}/$1'"; exit 1; }
      "${PYTHON}" "${PATCH_DIR}/$1" "$2"
    }
    apply_patch patch_nedelec_prism_pyramid.py    "${LIBMESH_SRC_DIR}"
    apply_patch cell_pyramid13_patch.py           "${LIBMESH_SRC_DIR}/src/geom/cell_pyramid13.C"
    apply_patch lagrange3D_shape_pyramid_patch.py "${LIBMESH_SRC_DIR}/src/fe/fe_lagrange_shape_3D.C"
  fi

  # Out-of-tree build, one build dir per method.
  local build_dir="${LIBMESH_SRC_DIR}/build-${pyTEC3D_LIBMESH_METHOD}"
  rm -rf "${build_dir}"; mkdir -p "${build_dir}"
  ( cd "${build_dir}"
    ../configure \
      --prefix="${LIBMESH_PREFIX}" \
      --with-methods="${pyTEC3D_LIBMESH_METHOD}" \
      --with-petsc="${PETSC_DIR}" \
      --with-metis="${METIS_DIR}" \
      --with-parmetis="${PARMETIS_DIR}" \
      --with-mpi="${MPI_DIR}" \
      --disable-vtk --disable-exodus --disable-libHilbert --disable-tetgen \
      --disable-triangle --disable-netgen --disable-poly2tri --disable-nlopt \
      --disable-gmv --disable-netcdf --enable-distmesh --disable-perflog \
      --disable-eigen --disable-fparser --disable-metaphysicl --enable-periodic \
      --disable-nanoflann --disable-qhull --disable-second --disable-capnproto \
      --disable-amr
    make -j"${CPU_COUNT}"
    make install
  )
  [ -f "${LIBMESH_LIB}" ] || { echo "ERROR: expected ${LIBMESH_LIB} after install."; exit 1; }
}

# =================================================================================
# Stage 2: pytec3d_solver (pybind11 module)
# =================================================================================
build_solver() {
  echo "==> Building pytec3d_solver"
  [ -f "${LIBMESH_LIB}" ] || { echo "ERROR: ${LIBMESH_LIB} not found."; exit 1; }
  mkdir -p "${SOLVER_PKG_DIR}"

  # Always start clean: no stale module, no stale CMake cache.
  rm -f "${SOLVER_PKG_DIR}"/pytec3d_solver*.{so,pyd,dylib}
  rm -rf "${SOLVER_BUILD_DIR}"

  # LIBMESH_DIR is the only libMesh input; CMake locates libmesh/libtimpi itself
  # and bakes their directories into the module's RPATH.
  cmake -S "${SRC_DIR}/src" -B "${SOLVER_BUILD_DIR}" -G Ninja \
    -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DPython3_EXECUTABLE="$(command -v "${PYTHON}")" \
    -DPYTEC3D_PYTHON_PACKAGE_ROOT="${PYTEC3D_PYTHON_PACKAGE_ROOT}" \
    -DPYTEC3D_RUNTIME_PREFIX="${PREFIX}" \
    -DLIBMESH_DIR="${LIBMESH_PREFIX}" \
    -DLIBMESH_METHOD="${pyTEC3D_LIBMESH_METHOD}"
  cmake --build "${SOLVER_BUILD_DIR}" -j"${CPU_COUNT}"
  cmake --install "${SOLVER_BUILD_DIR}"

  shopt -s nullglob
  local modules=("${SOLVER_PKG_DIR}"/pytec3d_solver*.{so,pyd,dylib})
  shopt -u nullglob
  [ "${#modules[@]}" -gt 0 ] || { echo "ERROR: build did not produce pytec3d_solver in ${SOLVER_PKG_DIR}."; exit 1; }

  # The module must resolve its libraries without LD_LIBRARY_PATH.
  if command -v ldd >/dev/null 2>&1; then
    local missing
    missing="$(env -u LD_LIBRARY_PATH ldd "${modules[0]}" | grep 'not found' || true)"
    if [ -n "${missing}" ]; then
      echo "ERROR: unresolved shared libraries (RPATH problem):"; echo "${missing}"; exit 1
    fi
  fi
  echo "Built ${modules[0]}"
}

echo "pyTEC3D build: method=${pyTEC3D_LIBMESH_METHOD} prefix=${PREFIX}"
build_libmesh
build_solver
echo "==> Done. Now run: pip install -e .   (or: pip install .)"