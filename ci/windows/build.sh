#!/bin/bash
# Windows build of the pyTEC3D solver, in an MSYS2 UCRT64 shell (GitHub Actions: .github/workflows/conda.yml).
#
# conda-forge has no Windows PETSc, and the PETSc of MSYS2 has no MUMPS. So PETSc is built here, as a static library,
# against the MUMPS, ScaLAPACK, METIS, ParMETIS and OpenBLAS packages of MSYS2 (MS-MPI builds), and libMesh as a static
# library against that PETSc. The module of each Python version links them, with the GCC and C++ runtimes statically;
# the DLLs of MSYS2 it needs (MUMPS, ScaLAPACK, OpenBLAS, the Fortran runtime, the Fortran interface of MS-MPI...) are
# copied next to it. msmpi.dll comes from the conda package msmpi, the rest from Windows. The modules are staged for
# conda/recipe.yaml (install_solver_windows.bat):
#
#   $STAGE_DIR/py<XY>/pytec3d_solver/{__init__.py, pytec3d_solver.cp<XY>-win_amd64.pyd, *.dll}
#
# Usage: build.sh [petsc|libmesh|solver|all]   (default all; the CI runs them as separate steps, to cache each)
#
# Environment:
#   DEPS_DIR    prefix of the PETSc and libMesh builds (cached between runs): $DEPS_DIR/{petsc,libmesh}
#   STAGE_DIR   receives the staged modules
#   PYTHONS     the python.exe of each conda environment to build for (Windows paths, separated by ';')
#   CPU_COUNT   parallel jobs (default: all cores)
#   PETSC_VERSION  PETSc release to build (default 3.26.0)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# Windows paths (C:/..., C:\...) as POSIX paths of MSYS2 (/c/...), for the shell and the autotools builds.
DEPS_DIR="$(cygpath -u "${DEPS_DIR:?prefix of the PETSc / libMesh builds}")"
STAGE_DIR="$(cygpath -u "${STAGE_DIR:?stage directory}")"
PYTHONS="${PYTHONS:?python.exe of each target Python, separated by ;}"
CPU_COUNT="${CPU_COUNT:-$(nproc)}"
PETSC_VERSION="${PETSC_VERSION:-3.26.0}"
PETSC_PREFIX="${DEPS_DIR}/petsc"
LIBMESH_PREFIX="${DEPS_DIR}/libmesh"
MINGW="${MINGW_PREFIX:-/ucrt64}"

[ "${MSYSTEM:-}" = "UCRT64" ] || { echo "ERROR: run in an MSYS2 UCRT64 shell (MSYSTEM=${MSYSTEM:-unset})."; exit 1; }
mkdir -p "${DEPS_DIR}/src" "${STAGE_DIR}"

# --- MPI: the compiler wrappers and the import libraries of MS-MPI of MSYS2 (mingw-w64-ucrt-x86_64-msmpi) ------------
for wrapper in mpicc mpicxx mpifort; do
  command -v "${wrapper}" >/dev/null || { echo "ERROR: ${wrapper} not found (pacman -S mingw-w64-ucrt-x86_64-msmpi)."; exit 1; }
done
MSMPI_LIBRARY="$(ls "${MINGW}"/lib/libmsmpi.dll.a "${MINGW}"/lib/libmsmpi.a 2>/dev/null | head -n1 || true)"
[ -n "${MSMPI_LIBRARY}" ] || { echo "ERROR: no MS-MPI library in ${MINGW}/lib (pacman -S mingw-w64-ucrt-x86_64-msmpi)."; exit 1; }
if ! command -v mpiexec >/dev/null; then
  msmpi_bin="$(cygpath -u "${MSMPI_BIN:-C:/Program Files/Microsoft MPI/Bin}")"
  export PATH="${PATH}:${msmpi_bin}"
fi

STAGE="${1:-all}"
case "${STAGE}" in petsc|libmesh|solver|all) ;; *) echo "ERROR: unknown stage '${STAGE}'."; exit 1 ;; esac
want() { [ "${STAGE}" = "all" ] || [ "${STAGE}" = "$1" ]; }

# --- PETSc, static, with the MUMPS of MSYS2 ---------------------------------------------------------------------------
if want petsc && [ ! -f "${PETSC_PREFIX}/lib/libpetsc.a" ]; then
  echo "==> Building PETSc ${PETSC_VERSION}"
  cd "${DEPS_DIR}/src"
  rm -rf "petsc-${PETSC_VERSION}"
  curl -fsSL "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz" | tar xz
  cd "petsc-${PETSC_VERSION}"
  # The import libraries (.dll.a) of the MSYS2 packages: their DLLs are shipped with the module.
  /usr/bin/python3 ./configure \
    --prefix="${PETSC_PREFIX}" \
    --with-shared-libraries=0 --with-debugging=0 \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort \
    COPTFLAGS=-O2 CXXOPTFLAGS=-O2 FOPTFLAGS=-O2 \
    --with-mpiexec=mpiexec \
    --with-blaslapack-lib="${MINGW}/lib/libopenblas.dll.a" \
    --with-mumps-include="${MINGW}/include" --with-mumps-lib="${MINGW}/lib/libmumps-dmo.dll.a" \
    --with-scalapack-include="${MINGW}/include" --with-scalapack-lib="${MINGW}/lib/libscalapack.dll.a" \
    --with-metis-include="${MINGW}/include" --with-metis-lib="${MINGW}/lib/libmetis.dll.a" \
    --with-parmetis-include="${MINGW}/include" --with-parmetis-lib="${MINGW}/lib/libparmetis.dll.a" \
    --with-x=0 --with-ssl=0 \
    || { echo "==> PETSc configure.log (end):"; tail -n 120 configure.log; exit 1; }
  make PETSC_DIR="$(pwd)" all -j"${CPU_COUNT}"
  make PETSC_DIR="$(pwd)" install
fi

# --- libMesh, static, against that PETSc (and its METIS / ParMETIS) ---------------------------------------------------
# --allow-multiple-definition: in a static libMesh, an explicit specialization (FE<2, SUBDIVISION>::side_map) and the
# explicit instantiation of its template are both definitions for the MinGW linker (on ELF the specialization wins);
# pyTEC3D does not use subdivision elements.
if want libmesh && [ ! -f "${LIBMESH_PREFIX}/lib/libmesh_opt.a" ]; then
  echo "==> Building libMesh"
  LDFLAGS="${LDFLAGS:-} -Wl,--allow-multiple-definition" \
  CC=mpicc CXX=mpicxx FC=mpifort F77=mpifort \
  PETSC_DIR="${PETSC_PREFIX}" MPI_DIR="${MINGW}" METIS_DIR=PETSc PARMETIS_DIR=PETSc \
  PYTHON=/usr/bin/python3 CPU_COUNT="${CPU_COUNT}" \
  LIBMESH_CONFIGURE_EXTRA="--disable-shared --enable-static" \
    "${ROOT}/solver/libmesh/build_libmesh.sh" "${DEPS_DIR}/src/libmesh" "${LIBMESH_PREFIX}" opt
fi

want solver || exit 0

# --- The solver module, for each Python --------------------------------------------------------------------------------
export PKG_CONFIG_PATH="${PETSC_PREFIX}/lib/pkgconfig:${MINGW}/lib/pkgconfig"
IFS=';' read -r -a pythons <<< "${PYTHONS}"
for python in "${pythons[@]}"; do
  [ -n "${python}" ] || continue
  python="$(cygpath -u "${python}")"
  tag="$("${python}" -c 'import sys; print(f"{sys.version_info[0]}{sys.version_info[1]}")' | tr -d '\r')"
  pybind11_dir="$("${python}" -m pybind11 --cmakedir | tr -d '\r')"
  build="${ROOT}/build/windows-py${tag}"
  echo "==> Building the solver for Python ${tag} (${python})"
  rm -rf "${build}"
  cmake -G Ninja -S "${ROOT}/solver" -B "${build}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DMPI_CXX_SKIP_MPICXX=ON -DMPI_CXX_INCLUDE_PATH="${MINGW}/include" \
    -DMPI_CXX_LIB_NAMES=msmpi -DMPI_msmpi_LIBRARY="${MSMPI_LIBRARY}" \
    -DCMAKE_PREFIX_PATH="${MINGW}" \
    -DLIBMESH_DIR="$(cygpath -m "${LIBMESH_PREFIX}")" -DLIBMESH_METHOD=opt \
    -DPYTEC3D_EXTRA_LIBS=z \
    -DPython3_EXECUTABLE="$(cygpath -m "${python}")" \
    -Dpybind11_DIR="$(cygpath -m "${pybind11_dir}")"
  cmake --build "${build}" -j"${CPU_COUNT}"
  rm -rf "${STAGE_DIR}/py${tag}"
  cmake --install "${build}" --prefix "${STAGE_DIR}/py${tag}" --component python

  # The DLLs of MSYS2 the module needs (transitively) go next to it; then everything must resolve to that directory,
  # Windows, MS-MPI (msmpi.dll) or Python.
  module="$(ls "${STAGE_DIR}/py${tag}"/pytec3d_solver/pytec3d_solver*.pyd)"
  package_dir="$(dirname "${module}")"
  ldd "${module}" | awk '{print $3}' | grep -i "^${MINGW}/bin/.*\.dll$" | sort -u | while read -r dll; do
    cp -v "${dll}" "${package_dir}/"
  done
  for binary in "${module}" "${package_dir}"/*.dll; do
    [ -f "${binary}" ] || continue
    foreign="$(PATH="${package_dir}:/c/Windows/System32" ldd "${binary}" \
      | grep -v -i -E '/c/windows/|msmpi\.dll|python3[0-9]*\.dll' | grep -v -F "${package_dir}" | grep -v '^\s*$' || true)"
    if [ -n "${foreign}" ]; then
      echo "ERROR: $(basename "${binary}") needs DLLs that are not shipped:"; echo "${foreign}"; exit 1
    fi
  done
  ls -l "${package_dir}"
  echo "Staged ${module}"
done
