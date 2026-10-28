#!/bin/bash
# conda (rattler-build) build of pytec3d-libmesh: solver/libmesh/build_libmesh.sh with the Open MPI wrappers of the host
# environment around the conda compilers, installed into the package prefix.
set -euxo pipefail

export OMPI_CC="${CC}" OMPI_CXX="${CXX}" OMPI_FC="${FC}"
export CC="${PREFIX}/bin/mpicc" CXX="${PREFIX}/bin/mpicxx" FC="${PREFIX}/bin/mpifort" F77="${PREFIX}/bin/mpifort"
export F90="${PREFIX}/bin/mpifort"
# configure runs small MPI programs (singletons): allow it in containers.
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_plm=isolated OMPI_MCA_rmaps_base_oversubscribe=yes
export PETSC_DIR="${PREFIX}" MPI_DIR="${PREFIX}" METIS_DIR="${PREFIX}" PARMETIS_DIR="${PREFIX}"
export PYTHON="${BUILD_PREFIX}/bin/python"

# LIBMESH_REF: the libMesh branch / commit (recipe.yaml passes it from the environment of the build; default master).
bash solver/libmesh/build_libmesh.sh "${SRC_DIR}/.libmesh-src" "${PREFIX}" opt

# The package is the library: no build tree, no libtool archives.
rm -f "${PREFIX}"/lib/libmesh_opt.la "${PREFIX}"/lib/libtimpi_opt.la
