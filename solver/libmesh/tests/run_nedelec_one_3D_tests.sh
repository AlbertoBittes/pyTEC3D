#!/bin/bash
# Build and run the NEDELEC_ONE 3D tests against an installed, patched libMesh.
#
# Usage: run_nedelec_one_3D_tests.sh <libMesh install prefix> <patched libMesh source dir> [build dir]
#   The source dir provides src/fe/fe_nedelec_one_shape_3D.C for the second derivative test, which
#   compiles it with second derivatives enabled (the solver's libMesh is built without them).
#   CXX: the MPI C++ wrapper libMesh was built with (default: libmesh-config --cxx, else mpicxx).
set -euo pipefail

PREFIX="$(cd "${1:?libMesh install prefix}" && pwd)"
SRC="$(cd "${2:?patched libMesh source dir}" && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${3:-${HERE}/build}"
LC="${PREFIX}/bin/libmesh-config"
CXX="${CXX:-$("${LC}" --cxx 2>/dev/null || true)}"
CXX="${CXX:-mpicxx}"
# shellcheck disable=SC2046  # flag lists
FLAGS=($("${LC}" --cppflags --cxxflags --include) -std=c++17)
# shellcheck disable=SC2207
LIBS=($("${LC}" --libs))

mkdir -p "${BUILD}"
"${CXX}" "${FLAGS[@]}" -o "${BUILD}/nedelec_one_3D_test" "${HERE}/nedelec_one_3D_test.C" "${LIBS[@]}"
"${CXX}" "${FLAGS[@]}" -DLIBMESH_ENABLE_SECOND_DERIVATIVES \
  -DNEDELEC_ONE_SHAPE_3D="\"${SRC}/src/fe/fe_nedelec_one_shape_3D.C\"" \
  -o "${BUILD}/nedelec_one_3D_second_deriv_test" "${HERE}/nedelec_one_3D_second_deriv_test.C" "${LIBS[@]}"

status=0
"${BUILD}/nedelec_one_3D_test" derivatives conformity convergence -ksp_type preonly -pc_type lu || status=1
"${BUILD}/nedelec_one_3D_second_deriv_test" || status=1
exit "${status}"
