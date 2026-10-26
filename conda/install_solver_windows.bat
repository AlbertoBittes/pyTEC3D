:: conda (rattler-build) build of pytec3d-solver on Windows: packages the solver compiled in MSYS2 (ci/windows/build.sh)
:: and staged in %PYTEC3D_WINDOWS_STAGE%\py<XY>\pytec3d_solver, with the metadata of the distribution pytec3d-solver.
@echo on
setlocal EnableDelayedExpansion
if "%PYTEC3D_WINDOWS_STAGE%"=="" (
  echo ERROR: PYTEC3D_WINDOWS_STAGE is not set: build the solver with ci\windows\build.sh first.
  exit 1
)
set "STAGE=%PYTEC3D_WINDOWS_STAGE:/=\%\py%PY_VER:.=%"
if not exist "%STAGE%\pytec3d_solver\__init__.py" (
  echo ERROR: no solver staged for Python %PY_VER% in %STAGE%.
  exit 1
)
robocopy "%STAGE%\pytec3d_solver" "%SP_DIR%\pytec3d_solver" /E /NFL /NDL /NJH /NJS
if %ERRORLEVEL% GEQ 8 exit 1

set "DIST=%SP_DIR%\pytec3d_solver-%PKG_VERSION%.dist-info"
mkdir "%DIST%"
(
  echo Metadata-Version: 2.1
  echo Name: pytec3d-solver
  echo Version: %PKG_VERSION%
  echo Summary: pyTEC3D C++ solver ^(libMesh, PETSc / MUMPS, MPI^) and its Python module pytec3d_solver
  echo Requires-Dist: numpy
) > "%DIST%\METADATA"
echo conda> "%DIST%\INSTALLER"
type nul > "%DIST%\RECORD"
exit 0
