cmake_minimum_required(VERSION 3.16)
project(pyTEC3D VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Dependencies ----------------------------------------------------------------
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(MPI REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC QUIET IMPORTED_TARGET PETSc)
if(NOT PETSC_FOUND)
  pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)
endif()

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

find_package(VTK REQUIRED
  COMPONENTS
    CommonCore
    CommonDataModel
    WrappingPythonCore
)

# --- libMesh -----------------------------------------------------------------------
# LIBMESH_DIR is the install prefix (contains include/ and lib/).
set(LIBMESH_DIR "" CACHE PATH "Path to the libMesh installation prefix.")
if(NOT LIBMESH_DIR AND DEFINED ENV{LIBMESH_DIR})
  set(LIBMESH_DIR "$ENV{LIBMESH_DIR}")
endif()
if(NOT LIBMESH_DIR)
  message(FATAL_ERROR "LIBMESH_DIR must be provided (e.g., -DLIBMESH_DIR=/path/to/libmesh).")
endif()

set(LIBMESH_METHOD "$ENV{pyTEC3D_LIBMESH_METHOD}" CACHE STRING "libMesh method to link against: opt, dbg, devel")
if(NOT LIBMESH_METHOD)
  set(LIBMESH_METHOD "opt")
endif()
string(TOLOWER "${LIBMESH_METHOD}" LIBMESH_METHOD)
if(NOT LIBMESH_METHOD MATCHES "^(opt|dbg|devel)$")
  message(FATAL_ERROR "Invalid LIBMESH_METHOD='${LIBMESH_METHOD}'. Allowed: opt, dbg, devel.")
endif()
message(STATUS "Using libMesh method: ${LIBMESH_METHOD}")

# NO_DEFAULT_PATH: never silently pick up a different libMesh from the env.
find_library(LIBMESH_LIB
  NAMES mesh_${LIBMESH_METHOD}
  PATHS "${LIBMESH_DIR}/lib"
  NO_DEFAULT_PATH
  REQUIRED
)
find_library(TIMPI_LIB
  NAMES timpi_${LIBMESH_METHOD}
  PATHS "${LIBMESH_DIR}/lib"
  NO_DEFAULT_PATH
  REQUIRED
)

# --- Sources / output location --------------------------------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)
get_filename_component(pyTEC3D_REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

# JSON definitions are compiled into the module so no runtime file lookup is needed.
set(PYTEC3D_RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources")
foreach(_res element_types fields domains)
  string(TOUPPER "${_res}" _res_var)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PYTEC3D_RESOURCE_DIR}/${_res}.json")
  file(READ "${PYTEC3D_RESOURCE_DIR}/${_res}.json" PYTEC3D_${_res_var}_JSON)
endforeach()
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/pytec3d_solver/core/embedded_resources.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/embedded_resources.cpp"
  @ONLY)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/embedded_resources.cpp")

set(PYTEC3D_PYTHON_PACKAGE_ROOT "${pyTEC3D_REPO_ROOT}/src" CACHE PATH
  "Root directory that contains pytec3d_* Python packages for extension module output")
set(PYTEC3D_SOLVER_PACKAGE_DIR "${PYTEC3D_PYTHON_PACKAGE_ROOT}/pytec3d_solver")

# --- Runtime library search path (RPATH) ----------------------------------------
# The module is copied as-is by pip (no relinking), so the RPATH must be baked in
# at build time. Absolute paths: libMesh/TIMPI install dirs + the env's lib/
# (PETSc, VTK, MPI, ...). Works for `pip install -e .` and `pip install .` on the
# machine that built libMesh. For relocatable wheels, bundle the libs and use
# auditwheel instead.
get_filename_component(LIBMESH_LIB_DIR "${LIBMESH_LIB}" DIRECTORY)
get_filename_component(TIMPI_LIB_DIR   "${TIMPI_LIB}"   DIRECTORY)

set(PYTEC3D_RUNTIME_PREFIX "$ENV{CONDA_PREFIX}" CACHE PATH
  "Env prefix whose lib/ holds PETSc/VTK/MPI at runtime")

set(PYTEC3D_RPATH "${LIBMESH_LIB_DIR}" "${TIMPI_LIB_DIR}")
if(PYTEC3D_RUNTIME_PREFIX)
  list(APPEND PYTEC3D_RPATH "${PYTEC3D_RUNTIME_PREFIX}/lib")
endif()
list(REMOVE_DUPLICATES PYTEC3D_RPATH)
message(STATUS "Module RPATH: ${PYTEC3D_RPATH}")

# --- Target ----------------------------------------------------------------------
set(pyTEC3D_INCLUDE_DIRS
  "${LIBMESH_DIR}/include"
  "${pyTEC3D_REPO_ROOT}/include"
  ${PETSC_INCLUDE_DIRS}
)

set(pyTEC3D_LINK_LIBS
  spdlog::spdlog_header_only
  nlohmann_json::nlohmann_json
  ${LIBMESH_LIB}
  ${TIMPI_LIB}
  PkgConfig::PETSC
  MPI::MPI_CXX
  Threads::Threads
  OpenMP::OpenMP_CXX
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::WrappingPythonCore
)

add_library(pytec3d_solver MODULE ${SOURCES})

set_target_properties(pytec3d_solver PROPERTIES
  PREFIX ""
  OUTPUT_NAME "pytec3d_solver"
  LIBRARY_OUTPUT_DIRECTORY "${PYTEC3D_SOLVER_PACKAGE_DIR}"
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH "${PYTEC3D_RPATH}"
)

# DT_RPATH (not DT_RUNPATH) so the search path also applies to transitive deps.
if(UNIX AND NOT APPLE)
  target_link_options(pytec3d_solver PRIVATE "LINKER:--disable-new-dtags")
endif()

target_include_directories(pytec3d_solver PRIVATE ${pyTEC3D_INCLUDE_DIRS})

target_compile_definitions(pytec3d_solver PRIVATE
  pyTEC3D_VERSION=${PROJECT_VERSION}
  pyTEC3D_BUILD_TYPE="$<CONFIG>"
)

if(PETSC_CFLAGS_OTHER)
  target_compile_options(pytec3d_solver PRIVATE ${PETSC_CFLAGS_OTHER})
endif()

target_link_libraries(pytec3d_solver
  PRIVATE
    pybind11::module
    ${pyTEC3D_LINK_LIBS}
)

vtk_module_autoinit(
  TARGETS pytec3d_solver
  MODULES
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::WrappingPythonCore
)