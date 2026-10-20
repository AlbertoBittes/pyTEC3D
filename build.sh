#!/bin/bash
# pyTEC3D build: libMesh -> solver module -> pip install of the package -> ansys-mapdl-reader.
#
# Usage: ./build.sh [--pip=editable|regular|none] [--extras=dev,notebooks]
#   --pip=editable   pip install -e .   (default: the environment imports src/ directly)
#   --pip=regular    pip install .      (copy into site-packages; rerun after every change of the sources)
#   --pip=none       no pip install (only libMesh and the module; ansys-mapdl-reader is still installed)
#   --extras=LIST    optional dependencies of pyproject.toml, e.g. --extras=dev,notebooks
#
# Env vars:
#   pyTEC3D_LIBMESH_METHOD=opt|dbg|devel   libMesh method (default: opt)
#   pyTEC3D_PIP_INSTALL, pyTEC3D_PIP_EXTRAS same as --pip / --extras (the flags take precedence)
#   pyTEC3D_CDB_READER=0                   do not install ansys-mapdl-reader (default 1: installed with
#                                          --no-deps, see install_cdb_reader)
#   FORCE_LIBMESH=1                        wipe and rebuild libMesh (default: reuse if present; also rebuilt
#                                          automatically when it links an MPI other than the one of PETSc)
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

# --- Options --------------------------------------------------------------------
usage() { sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
pyTEC3D_PIP_INSTALL="${pyTEC3D_PIP_INSTALL:-editable}"
pyTEC3D_PIP_EXTRAS="${pyTEC3D_PIP_EXTRAS:-}"
for arg in "$@"; do
  case "${arg}" in
    --pip=*)    pyTEC3D_PIP_INSTALL="${arg#--pip=}" ;;
    --extras=*) pyTEC3D_PIP_EXTRAS="${arg#--extras=}" ;;
    -h|--help)  usage; exit 0 ;;
    *) echo "ERROR: unknown argument '${arg}'."; usage; exit 1 ;;
  esac
done
case "${pyTEC3D_PIP_INSTALL}" in editable|regular|none) ;;
  *) echo "ERROR: invalid --pip='${pyTEC3D_PIP_INSTALL}'. Allowed: editable, regular, none."; exit 1 ;;
esac
if ! [[ "${pyTEC3D_PIP_EXTRAS}" =~ ^[A-Za-z0-9_,-]*$ ]]; then
  echo "ERROR: invalid --extras='${pyTEC3D_PIP_EXTRAS}' (comma-separated names, e.g. dev,notebooks)."; exit 1
fi

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

# --- MPI: one implementation only ------------------------------------------------
# libMesh, the solver module, PETSc and MUMPS must use the same MPI library. Other MPI installations on PATH or
# LD_LIBRARY_PATH (e.g. Intel oneAPI loaded by setvars.sh) would otherwise provide the mpicc/mpicxx wrappers, and the
# module would load two MPI libraries at run time (crash in the first parallel call, typically inside MUMPS).
# The MPI of the build is the one PETSc is linked to; its wrappers are taken from the same prefix.
mpi_libs_of() {  # <shared object> -> resolved MPI libraries (libmpi*.so*), one per line
  env -u LD_LIBRARY_PATH ldd "$1" 2>/dev/null | awk '/libmpi[a-z_]*\.so/ {print $3}' | sort -u
}
PETSC_LIB="$(ls "${PETSC_DIR:-${PREFIX}}"/lib/libpetsc.so 2>/dev/null || true)"
[ -n "${PETSC_LIB}" ] || { echo "ERROR: libpetsc.so not found in ${PETSC_DIR:-${PREFIX}}/lib."; exit 1; }
PETSC_MPI_DIR="$(mpi_libs_of "${PETSC_LIB}" | head -n1 | xargs -r dirname)"
[ -n "${PETSC_MPI_DIR}" ] || { echo "ERROR: cannot find the MPI library linked by ${PETSC_LIB}."; exit 1; }
MPI_PREFIX="$(dirname "${PETSC_MPI_DIR}")"
for w in mpicc mpicxx mpifort; do
  [ -x "${MPI_PREFIX}/bin/${w}" ] || { echo "ERROR: ${MPI_PREFIX}/bin/${w} not found (MPI of PETSc: ${PETSC_MPI_DIR})."; exit 1; }
done
# Drop from PATH / LD_LIBRARY_PATH the directories of other MPI installations (they provide mpicc or libmpi).
strip_foreign_mpi() {  # <colon list> <file to look for> -> list without the directories of other MPIs
  local out="" d
  IFS=':' read -r -a dirs <<< "$1"
  for d in "${dirs[@]}"; do
    [ -z "${d}" ] && continue
    if [ "${d}" != "${MPI_PREFIX}/bin" ] && [ "${d}" != "${PETSC_MPI_DIR}" ] && compgen -G "${d}/$2" >/dev/null; then
      echo "    (ignoring ${d}: another MPI installation)" >&2
      continue
    fi
    out="${out:+${out}:}${d}"
  done
  printf '%s' "${out}"
}
echo "==> MPI of PETSc: ${PETSC_MPI_DIR} (wrappers from ${MPI_PREFIX}/bin)"
export PATH="${MPI_PREFIX}/bin:$(strip_foreign_mpi "${PATH}" 'mpicc')"
export LD_LIBRARY_PATH="$(strip_foreign_mpi "${LD_LIBRARY_PATH:-}" 'libmpi.so*')"
unset I_MPI_ROOT I_MPI_CC I_MPI_CXX I_MPI_F90 I_MPI_FC MPICC MPICXX 2>/dev/null || true

# Fail when a library links an MPI other than the one of PETSc (or two MPI libraries at once).
check_mpi() {  # <shared object>
  local libs foreign
  libs="$(mpi_libs_of "$1")"
  foreign="$(printf '%s\n' "${libs}" | grep -v "^${PETSC_MPI_DIR}/" | grep -v '^$' || true)"
  if [ -n "${foreign}" ]; then
    echo "ERROR: $1 links MPI libraries outside ${PETSC_MPI_DIR} (the MPI of PETSc):"; echo "${foreign}"
    return 1
  fi
}

# --- Toolchain / flags -----------------------------------------------------------
export CC="${CC:-${MPI_PREFIX}/bin/mpicc}" CXX="${CXX:-${MPI_PREFIX}/bin/mpicxx}"
export FC="${FC:-${MPI_PREFIX}/bin/mpifort}" F77="${F77:-${MPI_PREFIX}/bin/mpifort}" F90="${F90:-${MPI_PREFIX}/bin/mpifort}"
export PETSC_DIR="${PETSC_DIR:-${PREFIX}}"
export METIS_DIR="${METIS_DIR:-${PREFIX}}"
export PARMETIS_DIR="${PARMETIS_DIR:-${PREFIX}}"
export MPI_DIR="${MPI_DIR:-${MPI_PREFIX}}"

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
    if check_mpi "${LIBMESH_LIB}" >/dev/null; then
      echo "==> libMesh (${pyTEC3D_LIBMESH_METHOD}) already installed in ${LIBMESH_PREFIX}; skipping (FORCE_LIBMESH=1 to rebuild)."
      return
    fi
    echo "==> The installed libMesh links another MPI than PETSc: rebuilding it (same patched checkout)."
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
  check_mpi "${LIBMESH_LIB}" || exit 1
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
  check_mpi "${modules[0]}" || exit 1
  echo "Built ${modules[0]} (MPI: ${PETSC_MPI_DIR})"
}

# =================================================================================
# Stage 3: the Python packages (pytec3d + the module just built)
# =================================================================================
pip_install() {
  local target="${SRC_DIR}${pyTEC3D_PIP_EXTRAS:+[${pyTEC3D_PIP_EXTRAS}]}"
  case "${pyTEC3D_PIP_INSTALL}" in
    none)     echo "==> pip install: skipped (--pip=none)"; return ;;
    editable) echo "==> pip install -e ${target}"; "${PYTHON}" -m pip install -e "${target}" ;;
    regular)  echo "==> pip install ${target}";    "${PYTHON}" -m pip install "${target}" ;;
  esac
  # Imported from outside the repository, so that the installed copy is checked, not ./src.
  ( cd / && "${PYTHON}" -c "import pytec3d, pytec3d_solver.pytec3d_solver; print('    pytec3d imported from', pytec3d.__file__)" ) \
    || { echo "ERROR: pytec3d does not import after pip install."; exit 1; }
}

# =================================================================================
# Stage 4: ansys-mapdl-reader (conversion of ANSYS .cdb meshes, pytec3d-convert-cdb)
# =================================================================================
# Its declared requirements (vtk / numpy / pyvista pins) conflict with the conda packages the solver is built against;
# the package itself works with them. It is installed without dependencies, plus its two small pure-Python modules,
# after the pip install of pytec3d, so that no later resolution of dependencies touches it.
install_cdb_reader() {
  if [ "${pyTEC3D_CDB_READER:-1}" = "0" ]; then
    echo "==> ansys-mapdl-reader: skipped (pyTEC3D_CDB_READER=0)"
    return
  fi
  echo "==> Installing ansys-mapdl-reader (--no-deps)"
  "${PYTHON}" -m pip install --no-deps --quiet ansys-mapdl-reader appdirs tqdm
  "${PYTHON}" -c "from ansys.mapdl import reader" \
    || { echo "ERROR: ansys-mapdl-reader does not import in this environment."; exit 1; }
}

echo "pyTEC3D build: method=${pyTEC3D_LIBMESH_METHOD} prefix=${PREFIX} pip=${pyTEC3D_PIP_INSTALL}${pyTEC3D_PIP_EXTRAS:+ extras=${pyTEC3D_PIP_EXTRAS}}"
build_libmesh
build_solver
pip_install
install_cdb_reader
if [ "${pyTEC3D_PIP_INSTALL}" = "none" ]; then
  echo "==> Done. Now run: pip install -e .   (or: pip install .)"
else
  echo "==> Done."
fi