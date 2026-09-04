#!/bin/bash
set -exuo pipefail

SRC_DIR=${SRC_DIR:-$PWD}
SRC_DIR="$(cd "${SRC_DIR}" && pwd)"
CORE_BUILD_SCRIPT="${SRC_DIR}/conda-recipe/recipes/apec-m-core/build.sh"
PYTHON=${PYTHON:-python}

if [ ! -d "${SRC_DIR}/src" ] || [ ! -d "${SRC_DIR}/python/apec_m_core" ]; then
  echo "ERROR: SRC_DIR='${SRC_DIR}' is not the project root."
  exit 1
fi

if [ ! -f "${CORE_BUILD_SCRIPT}" ]; then
  echo "ERROR: Missing core build script '${CORE_BUILD_SCRIPT}'."
  exit 1
fi

# Build and install libMesh + solver extension into the active PREFIX.
bash "${CORE_BUILD_SCRIPT}"

"${PYTHON}" - <<'PY'
import apec_m_core.apec_m_solver
print("Installed apec-m solver package successfully")
PY
