#!/bin/bash
# Build and install the MFEM the pyTEC3D solver links: MFEM_REF (default: the pinned release), parallel (MPI, hypre,
# METIS), shared, configured with only what the solver uses. The commit built is written to
# <install prefix>/share/pytec3d/mfem-commit.
#
# Usage: build_mfem.sh <source dir> <install prefix> [opt|dbg]
#
# Environment:
#   MFEM_REF                        tag, branch or commit of MFEM to build (default: v4.9)
#   HYPRE_DIR, METIS_DIR, MPI_DIR   the dependencies (default: CONDA_PREFIX / PREFIX)
#   CC, CXX                         the MPI compiler wrappers (default: <MPI_DIR>/bin/mpicc, mpicxx)
#   CPU_COUNT                       parallel jobs (default: all cores)
#   MFEM_CMAKE_EXTRA                more CMake options
set -euo pipefail

MFEM_REPOSITORY="https://github.com/mfem/mfem.git"
MFEM_REF="${MFEM_REF:-v4.9}"

SRC_DIR="${1:?source directory}"
INSTALL_PREFIX="${2:?install prefix}"
METHOD="${3:-opt}"
case "${METHOD}" in opt) BUILD_TYPE=Release ;; dbg) BUILD_TYPE=Debug ;;
  *) echo "ERROR: invalid method '${METHOD}' (opt, dbg)."; exit 1 ;; esac

DEPS_PREFIX="${CONDA_PREFIX:-${PREFIX:-}}"
HYPRE_DIR="${HYPRE_DIR:-${DEPS_PREFIX}}"
METIS_DIR="${METIS_DIR:-${DEPS_PREFIX}}"
MPI_DIR="${MPI_DIR:-${DEPS_PREFIX}}"
# The MPI wrappers of the environment: Intel oneAPI may be first on PATH and would mix its MPI with the env's.
export CC="${CC:-${MPI_DIR}/bin/mpicc}" CXX="${CXX:-${MPI_DIR}/bin/mpicxx}"
CPU_COUNT="${CPU_COUNT:-$(nproc 2>/dev/null || echo 2)}"

if [ ! -e "${SRC_DIR}/.git" ]; then
  echo "==> Cloning MFEM ${MFEM_REF} into ${SRC_DIR}"
  mkdir -p "${SRC_DIR}"
  git -C "${SRC_DIR}" init -q
  git -C "${SRC_DIR}" remote add origin "${MFEM_REPOSITORY}"
  git -C "${SRC_DIR}" fetch -q --depth 1 origin "${MFEM_REF}"
  git -C "${SRC_DIR}" checkout -q FETCH_HEAD
fi

BUILD_DIR="${SRC_DIR}/build-${METHOD}"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}" "${INSTALL_PREFIX}"
# shellcheck disable=SC2086  # MFEM_CMAKE_EXTRA is a list of options
cmake -S "${SRC_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
  -DCMAKE_C_COMPILER="${CC}" -DCMAKE_CXX_COMPILER="${CXX}" \
  -DBUILD_SHARED_LIBS=ON \
  -DMFEM_USE_MPI=ON \
  -DHYPRE_DIR="${HYPRE_DIR}" \
  -DMFEM_USE_METIS_5=ON -DMETIS_DIR="${METIS_DIR}" \
  -DMFEM_ENABLE_EXAMPLES=OFF -DMFEM_ENABLE_MINIAPPS=OFF -DMFEM_ENABLE_TESTING=OFF \
  ${MFEM_CMAKE_EXTRA:-}
cmake --build "${BUILD_DIR}" -j "${CPU_COUNT}"
cmake --install "${BUILD_DIR}"

mkdir -p "${INSTALL_PREFIX}/share/pytec3d"
git -C "${SRC_DIR}" rev-parse HEAD > "${INSTALL_PREFIX}/share/pytec3d/mfem-commit"
echo "==> MFEM $(cat "${INSTALL_PREFIX}/share/pytec3d/mfem-commit") installed in ${INSTALL_PREFIX}"
