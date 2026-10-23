#!/bin/bash
# Windows build of pyTEC3D under MSYS2 UCRT64 (feasibility): patched libMesh -> pytec3d_solver.pyd -> smoke test.
#
# Usage (MSYS2 UCRT64 shell, packages of .github/workflows/windows-msys2.yml installed):
#   ci/msys2/build.sh libmesh|solver|test|all
#
# PETSc is the MSYS2 package, flavour dmo (double, real, MS-MPI, optimised; METIS and ParMETIS), through its build tree
# (petsc-build: PETSC_DIR=$MINGW_PREFIX/src/petsc, PETSC_ARCH=dmo). It has no MUMPS: enough to build and import the
# module, not for the default direct solver. libMesh is a static library, linked into the module.
#
# Env vars:
#   PYTEC3D_DEPS   install folder of libMesh (default: $HOME/deps); an installed libMesh is skipped
set -euo pipefail

STAGE="${1:-all}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEPS="${PYTEC3D_DEPS:-${HOME}/deps}"
PY="${MINGW_PREFIX}/bin/python"
JOBS="$(nproc)"

# The libMesh of the Linux build (.deps/src/libmesh); the PETSc of MSYS2.
LIBMESH_COMMIT=9e5b0b870a585427a1bd400fec715d7789c4d976
export PETSC_DIR="${MINGW_PREFIX}/src/petsc" PETSC_ARCH=dmo
LIBMESH_PREFIX="${DEPS}/libmesh-opt"
SOLVER_BUILD_DIR="${ROOT}/build/msys2-opt"

mkdir -p "${DEPS}/src"
echo "pyTEC3D MSYS2 build: stage=${STAGE} root=${ROOT} deps=${DEPS} jobs=${JOBS} (${MSYSTEM:-?}) PETSc ${PETSC_DIR} ${PETSC_ARCH}"

build_libmesh() {
  if [ -f "${LIBMESH_PREFIX}/lib/libmesh_opt.a" ]; then echo "==> libMesh: installed, skipped"; return; fi
  echo "==> libMesh ${LIBMESH_COMMIT}"
  local src="${DEPS}/src/libmesh"
  rm -rf "${src}"
  git clone --filter=blob:none https://github.com/libMesh/libmesh.git "${src}"
  git -C "${src}" checkout --quiet "${LIBMESH_COMMIT}"
  git -C "${src}" submodule update --init --recursive
  "${PY}" "${ROOT}/libmesh_patches/patch_nedelec_prism_pyramid.py" "${src}"
  "${PY}" "${ROOT}/libmesh_patches/cell_pyramid13_patch.py" "${src}/src/geom/cell_pyramid13.C"
  "${PY}" "${ROOT}/libmesh_patches/lagrange3D_shape_pyramid_patch.py" "${src}/src/fe/fe_lagrange_shape_3D.C"
  mkdir -p "${src}/build-opt"
  cd "${src}/build-opt"
  # The options of build.sh, plus static libraries (no DLL exports to deal with) and METIS / ParMETIS from PETSc.
  ../configure --prefix="${LIBMESH_PREFIX}" --with-methods=opt \
    CC=mpicc CXX=mpicxx FC=mpifort F77=mpifort \
    --disable-shared --enable-static \
    --with-petsc="${PETSC_DIR}" --with-metis=PETSc --with-parmetis=PETSc \
    --disable-vtk --disable-exodus --disable-libHilbert --disable-tetgen \
    --disable-triangle --disable-netgen --disable-poly2tri --disable-nlopt \
    --disable-gmv --disable-netcdf --enable-distmesh --disable-perflog \
    --disable-eigen --disable-fparser --disable-metaphysicl --enable-periodic \
    --disable-nanoflann --disable-qhull --disable-second --disable-capnproto \
    --disable-amr
  make -j"${JOBS}"
  make install
  [ -f "${LIBMESH_PREFIX}/lib/libmesh_opt.a" ] || { echo "ERROR: no libmesh_opt.a in ${LIBMESH_PREFIX}/lib"; exit 1; }
}

build_solver() {
  echo "==> pytec3d_solver"
  # The module finds PETSc as the pkg-config package petsc: an alias of the flavour's petsc-dmo.pc.
  local pc_dir="${SOLVER_BUILD_DIR}-pkgconfig"
  mkdir -p "${pc_dir}"
  cp "${MINGW_PREFIX}/lib/pkgconfig/petsc-${PETSC_ARCH}.pc" "${pc_dir}/petsc.pc"
  export PKG_CONFIG_PATH="${pc_dir}:${PKG_CONFIG_PATH:-}"
  # Static libMesh: its dependencies (PETSc included) go at the end of the link line.
  local static_libs
  static_libs="$("${LIBMESH_PREFIX}/bin/libmesh-config" --libs)"
  echo "    static link libraries: ${static_libs}"
  rm -rf "${SOLVER_BUILD_DIR}"
  cmake -S "${ROOT}/src" -B "${SOLVER_BUILD_DIR}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx \
    -DCMAKE_PREFIX_PATH="${MINGW_PREFIX}" \
    -DPython3_EXECUTABLE="${PY}" \
    -DLIBMESH_DIR="${LIBMESH_PREFIX}" -DLIBMESH_METHOD=opt \
    -DCMAKE_CXX_STANDARD_LIBRARIES="${static_libs}"
  cmake --build "${SOLVER_BUILD_DIR}" -j"${JOBS}"
  cmake --install "${SOLVER_BUILD_DIR}"
  ls -l "${ROOT}/src/pytec3d_solver/"
  echo "    DLLs the module needs:"
  ldd "${ROOT}"/src/pytec3d_solver/pytec3d_solver*.pyd | grep -vi "/c/windows" || true
}

smoke_test() {
  echo "==> Smoke test"
  cd /
  PYTHONPATH="${ROOT}/src" "${PY}" -c "
import pytec3d_solver.pytec3d_solver as m
print('module:', m.__file__)
import pytec3d
print('pytec3d', pytec3d.__version__, 'from', pytec3d.__file__)
"
}

case "${STAGE}" in
  libmesh) build_libmesh ;;
  solver) build_solver ;;
  test) smoke_test ;;
  all) build_libmesh; build_solver; smoke_test ;;
  *) echo "ERROR: unknown stage '${STAGE}'"; exit 1 ;;
esac
