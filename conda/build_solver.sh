#!/bin/bash
# conda (rattler-build) build of pytec3d-solver on Linux: pip install ./solver (scikit-build-core, CMake) against the
# MFEM of pytec3d-mfem and the hypre / PETSc / Open MPI of the host environment.
set -euxo pipefail

export OMPI_CC="${CC}" OMPI_CXX="${CXX}"
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export CMAKE_ARGS="${CMAKE_ARGS:-} -DMPI_CXX_COMPILER=${PREFIX}/bin/mpicxx -DCMAKE_PREFIX_PATH=${PREFIX}"

"${PYTHON}" -m pip install ./solver -vv --no-deps --no-build-isolation \
  -Cbuild-dir=build-solver \
  -Ccmake.define.MFEM_DIR="${PREFIX}" \
  -Ccmake.define.PYTEC3D_RUNTIME_PREFIX="${PREFIX}" \
  -Ccmake.define.Python3_EXECUTABLE="${PYTHON}"
