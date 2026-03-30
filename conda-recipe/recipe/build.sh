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
export VTK_PREFIX=${VTK_PREFIX:-${PREFIX}}
export VTK_CMAKE_DIR=${VTK_CMAKE_DIR:-${PREFIX}/lib/cmake/vtk-9.3}
export VTK_INCLUDE_DIR=${VTK_INCLUDE_DIR:-${PREFIX}/include/vtk-9.3}
export BOOST_DIR=${BOOST_DIR:-${PREFIX}}
export LIBMESH_DIR=${LIBMESH_PREFIX}
export CPPFLAGS="${CPPFLAGS:-} -I${PREFIX}/include"
export CXXFLAGS="${CXXFLAGS:-} -I${PREFIX}/include"
export LDFLAGS="${LDFLAGS:-} -L${PREFIX}/lib"
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

export CPPFLAGS="${CPPFLAGS} -I${VTK_INCLUDE_DIR}"
export CXXFLAGS="${CXXFLAGS} -I${VTK_INCLUDE_DIR}"

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
    --with-vtk="${VTK_PREFIX}" \
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

# 2. Build the apec_m_solver pybind module with CMake
mkdir -p "${SRC_DIR}/build_apec_m"
pushd "${SRC_DIR}/build_apec_m"
cmake "${SRC_DIR}/src" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DLIBMESH_DIR="${LIBMESH_PREFIX}" \
    -DPETSC_DIR="${PETSC_DIR}" \
    -DVTK_DIR="${VTK_CMAKE_DIR}" \
    -DBOOST_ROOT="${BOOST_DIR}" \
    -DBOOST_INCLUDEDIR="${BOOST_DIR}/include" \
    -DBOOST_LIBRARYDIR="${BOOST_DIR}/lib" \
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
"${PYTHON}" -m pip install --no-build-isolation --no-deps --prefix "${PREFIX}" .
popd

# 4. Install ansys-mapdl-reader from PyPI (py312-compatible release).
"${PYTHON}" -m pip install --prefix "${PREFIX}" "ansys-mapdl-reader==0.55.2"

# 5. Vendor openpyxl into repository-local python/apec_m/local_packages for ParaView usage.
"${PYTHON}" -m pip install --upgrade openpyxl et_xmlfile --target "${SRC_DIR}/python/apec_m/local_packages"
