#!/bin/bash
# conda (rattler-build) build of pytec3d-mfem: solver/mfem/build_mfem.sh with the Open MPI wrappers of the host
# environment around the conda compilers, against the hypre and METIS of the host environment, installed into the
# package prefix.
set -euxo pipefail

export OMPI_CC="${CC}" OMPI_CXX="${CXX}"
export CC="${PREFIX}/bin/mpicc" CXX="${PREFIX}/bin/mpicxx"
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export HYPRE_DIR="${PREFIX}" METIS_DIR="${PREFIX}" MPI_DIR="${PREFIX}"
export MFEM_CMAKE_EXTRA="-DCMAKE_PREFIX_PATH=${PREFIX}"

# MFEM_REF: the MFEM tag / commit (recipe.yaml passes it from the environment of the build; default the pinned release).
bash solver/mfem/build_mfem.sh "${SRC_DIR}/.mfem-src" "${PREFIX}" opt
