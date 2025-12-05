#!/bin/bash
set -exuo pipefail

SRC_DIR=${SRC_DIR:-$PWD}

LIBMESH_PREFIX=${LIBMESH_PREFIX:-${PREFIX}}
LIBMESH_SRC_DIR=${SRC_DIR}/libmesh
PYTHON=${PYTHON:-python}

export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export FC=${FC:-mpifort}
export F77=${F77:-mpifort}
export F90=${F90:-mpifort}
export PETSC_DIR=${PETSC_DIR:-${PREFIX}}
export METIS_DIR=${METIS_DIR:-${PREFIX}}
export PARMETIS_DIR=${PARMETIS_DIR:-${PREFIX}}
export MPI_DIR=${MPI_DIR:-${PREFIX}}
export VTK_DIR=${VTK_DIR:-${PREFIX}/lib/cmake/vtk}
export BOOST_DIR=${BOOST_DIR:-${PREFIX}}
export LIBMESH_DIR=${LIBMESH_PREFIX}
export CPPFLAGS="${CPPFLAGS:-} -I${PREFIX}/include"
export CXXFLAGS="${CXXFLAGS:-} -I${PREFIX}/include"
export LDFLAGS="${LDFLAGS:-} -L${PREFIX}/lib"
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

# 1. Clone and build libMesh with the requested options
rm -rf "${LIBMESH_SRC_DIR}"
git clone --depth 1 https://github.com/libMesh/libmesh.git "${LIBMESH_SRC_DIR}"
pushd "${LIBMESH_SRC_DIR}"
git submodule update --init --recursive
mkdir build
pushd build
../configure \
    --prefix="${LIBMESH_PREFIX}" \
    --with-petsc="${PETSC_DIR}" \
    --with-metis="${METIS_DIR}" \
    --with-parmetis="${PARMETIS_DIR}" \
    --with-mpi="${MPI_DIR}" \
    --with-vtk="${VTK_DIR}" \
    --with-boost="${BOOST_DIR}" \
    --disable-strict-lgpl \
    --enable-distmesh \
    --enable-ifem \
    --enable-perflog \
    --disable-eigen \
    --disable-amr
make -j"${CPU_COUNT}"
make install
popd
popd

# 2. Build and install apec_m with CMake
mkdir -p "${SRC_DIR}/build_apec_m"
pushd "${SRC_DIR}/build_apec_m"
cmake "${SRC_DIR}/src" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DLIBMESH_DIR="${LIBMESH_PREFIX}" \
    -DPETSC_DIR="${PETSC_DIR}" \
    -DVTK_DIR="${VTK_DIR}" \
    -DBOOST_ROOT="${BOOST_DIR}" \
    -DBOOST_INCLUDEDIR="${BOOST_DIR}/include" \
    -DBOOST_LIBRARYDIR="${BOOST_DIR}/lib"
ninja -j"${CPU_COUNT}"
install -Dm755 "${SRC_DIR}/bin/apec_m" "${PREFIX}/bin/apec_m"
popd

# 3. Install the Python module (and its pip-managed dependencies)
pushd "${SRC_DIR}"
"${PYTHON}" -m pip install --no-build-isolation --prefix "${PREFIX}" .
popd
