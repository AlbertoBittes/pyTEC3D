#!/bin/bash
# Build and install the libMesh the pyTEC3D solver links: the libMesh development branch (or LIBMESH_REF), patched
# (patches/), configured with only what the solver uses. Shared by build.sh, the conda recipe (conda/recipe.yaml) and
# the Windows CI (MSYS2). The commit built is written to <install prefix>/share/pytec3d/libmesh-commit.
#
# Usage: build_libmesh.sh <source dir> <install prefix> [opt|dbg|devel]
#   <source dir>      cloned (LIBMESH_REF, with submodules) and patched when it has no checkout yet; a checkout is
#                     patched once (marker file .pytec3d_patched)
#   <install prefix>  receives include/, lib/ (libmesh_<method>, timpi_<method>), bin/
#
# Environment:
#   LIBMESH_REF                                   branch, tag or commit of libMesh to build (default: master)
#   PETSC_DIR, MPI_DIR, METIS_DIR, PARMETIS_DIR   the dependencies (default: CONDA_PREFIX / PREFIX)
#   CC, CXX, FC                                   the compilers (MPI wrappers, or compilers with MPI_DIR)
#   PYTHON                                        runs the patch scripts (default: python3)
#   CPU_COUNT                                     parallel jobs (default: all cores)
#   LIBMESH_CONFIGURE_EXTRA                       more configure options (Windows: --disable-shared --enable-static)
set -euo pipefail

LIBMESH_REPOSITORY="https://github.com/libMesh/libmesh.git"
# The solver follows the libMesh development: master by default. The patches stop the build if they no longer apply.
LIBMESH_REF="${LIBMESH_REF:-master}"

SRC_DIR="${1:?source directory}"
INSTALL_PREFIX="${2:?install prefix}"
METHOD="${3:-opt}"
case "${METHOD}" in opt|dbg|devel) ;; *) echo "ERROR: invalid method '${METHOD}' (opt, dbg, devel)."; exit 1 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPS_PREFIX="${CONDA_PREFIX:-${PREFIX:-}}"
PYTHON="${PYTHON:-python3}"
CPU_COUNT="${CPU_COUNT:-$(nproc 2>/dev/null || echo 2)}"

# --- Source: LIBMESH_REF with its submodules, patched once -------------------------------------------------------------
if [ ! -e "${SRC_DIR}/.git" ]; then
  echo "==> Cloning libMesh ${LIBMESH_REF} into ${SRC_DIR}"
  mkdir -p "${SRC_DIR}"
  git -C "${SRC_DIR}" init -q
  git -C "${SRC_DIR}" remote add origin "${LIBMESH_REPOSITORY}"
  git -C "${SRC_DIR}" fetch -q --depth 1 origin "${LIBMESH_REF}"
  git -C "${SRC_DIR}" checkout -q FETCH_HEAD
  git -C "${SRC_DIR}" submodule update -q --init --recursive --jobs 8
  # Git on Windows (MSYS2) may check out the symbolic links of the tree as text files.
  if [ -n "${MSYSTEM:-}" ] && [ -x "${SRC_DIR}/contrib/bin/fix_windows_symlinks.sh" ]; then
    ( cd "${SRC_DIR}" && ./contrib/bin/fix_windows_symlinks.sh )
  fi
fi
if [ ! -f "${SRC_DIR}/.pytec3d_patched" ]; then
  echo "==> Patching libMesh"
  "${PYTHON}" "${HERE}/patches/patch_nedelec_prism_pyramid.py"    "${SRC_DIR}"
  "${PYTHON}" "${HERE}/patches/cell_pyramid13_patch.py"           "${SRC_DIR}/src/geom/cell_pyramid13.C"
  "${PYTHON}" "${HERE}/patches/lagrange3D_shape_pyramid_patch.py" "${SRC_DIR}/src/fe/fe_lagrange_shape_3D.C"
  touch "${SRC_DIR}/.pytec3d_patched"
fi

# --- Configure, build, install (out of tree, one build directory per method) -----------------------------------------
# Only what the solver uses; no Boost (an external Boost of the system would put its include directory, e.g.
# /usr/include, before the sysroot of the conda compilers).
BUILD_DIR="${SRC_DIR}/build-${METHOD}"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}" "${INSTALL_PREFIX}"
cd "${BUILD_DIR}"
# shellcheck disable=SC2086  # LIBMESH_CONFIGURE_EXTRA is a list of options
../configure \
  --prefix="${INSTALL_PREFIX}" \
  --with-methods="${METHOD}" \
  --with-petsc="${PETSC_DIR:-${DEPS_PREFIX}}" \
  --with-metis="${METIS_DIR:-${DEPS_PREFIX}}" \
  --with-parmetis="${PARMETIS_DIR:-${DEPS_PREFIX}}" \
  --with-mpi="${MPI_DIR:-${DEPS_PREFIX}}" \
  --disable-vtk --disable-exodus --disable-libHilbert --disable-tetgen \
  --disable-triangle --disable-netgen --disable-poly2tri --disable-nlopt \
  --disable-gmv --disable-netcdf --enable-distmesh --disable-perflog \
  --disable-eigen --disable-fparser --disable-metaphysicl --enable-periodic \
  --disable-nanoflann --disable-qhull --disable-second --disable-capnproto \
  --disable-amr --disable-examples --disable-boost \
  ${LIBMESH_CONFIGURE_EXTRA:-}
make -j"${CPU_COUNT}"
make install
COMMIT="$(git -C "${SRC_DIR}" rev-parse HEAD)"
mkdir -p "${INSTALL_PREFIX}/share/pytec3d"
echo "${COMMIT}" > "${INSTALL_PREFIX}/share/pytec3d/libmesh-commit"
echo "==> libMesh ${COMMIT} (${METHOD}) installed in ${INSTALL_PREFIX}"
