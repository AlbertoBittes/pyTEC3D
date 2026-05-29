cmake_minimum_required(VERSION 3.16)
project(ApecM VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(MPI REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC QUIET IMPORTED_TARGET PETSc)
if(NOT PETSC_FOUND)
  pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET petsc)
endif()
option(APEC_M_BUILD_PYBIND_MODULE "Build the Python solver extension module (pybind11)." ON)
option(APEC_M_ENABLE_SANITIZERS "Enable AddressSanitizer+UBSan for debug checks." OFF)
option(APEC_M_ENABLE_IMPORT_TEST "Enable Python import smoke test for the solver module." ON)
option(APEC_M_ENABLE_VALGRIND_SMOKE "Add a valgrind smoke target for module import." OFF)

if(APEC_M_BUILD_PYBIND_MODULE)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG REQUIRED)
endif()

find_package(VTK REQUIRED
  COMPONENTS
    CommonCore
    CommonDataModel
    CommonExecutionModel
    FiltersCore
    IOXML
    WrappingPythonCore
)

set(LIBMESH_DIR "${LIBMESH_DIR}" CACHE PATH "Path to the libMesh installation.")
if(NOT LIBMESH_DIR AND DEFINED ENV{LIBMESH_DIR})
  set(LIBMESH_DIR "$ENV{LIBMESH_DIR}")
endif()
if(NOT LIBMESH_DIR)
  message(FATAL_ERROR "LIBMESH_DIR must be provided (e.g., -DLIBMESH_DIR=/path/to/libmesh).")
endif()

set(LIBMESH_METHOD "$ENV{APEC_M_LIBMESH_METHOD}" CACHE STRING "libMesh method to link against: opt, dbg, devel")
if(NOT LIBMESH_METHOD)
  set(LIBMESH_METHOD "opt")
endif()
string(TOLOWER "${LIBMESH_METHOD}" LIBMESH_METHOD)
if(NOT LIBMESH_METHOD STREQUAL "opt" AND
   NOT LIBMESH_METHOD STREQUAL "dbg" AND
   NOT LIBMESH_METHOD STREQUAL "devel")
  message(FATAL_ERROR "Invalid LIBMESH_METHOD='${LIBMESH_METHOD}'. Allowed: opt, dbg, devel.")
endif()
message(STATUS "Using libMesh method: ${LIBMESH_METHOD}")

find_library(LIBMESH_LIB
  NAMES mesh_${LIBMESH_METHOD}
  HINTS "${LIBMESH_DIR}/lib" "${LIBMESH_DIR}/lib64"
  REQUIRED
)
find_library(TIMPI_LIB
  NAMES timpi_${LIBMESH_METHOD}
  HINTS "${LIBMESH_DIR}/lib" "${LIBMESH_DIR}/lib64"
  REQUIRED
)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)
get_filename_component(APEC_M_REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

set(APEC_M_INCLUDE_DIRS
  "${LIBMESH_DIR}/include"
  "${APEC_M_REPO_ROOT}/include"
  ${PETSC_INCLUDE_DIRS}
)

if(EXISTS "${APEC_M_REPO_ROOT}/include/dependencies")
  list(APPEND APEC_M_INCLUDE_DIRS "${APEC_M_REPO_ROOT}/include/dependencies")
endif()

set(APEC_M_LINK_LIBS
  spdlog::spdlog_header_only
  ${LIBMESH_LIB}
  ${TIMPI_LIB}
  PkgConfig::PETSC
  MPI::MPI_CXX
  Threads::Threads
  OpenMP::OpenMP_CXX
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::CommonExecutionModel
  VTK::FiltersCore
  VTK::IOXML
  VTK::WrappingPythonCore
)

if(APEC_M_BUILD_PYBIND_MODULE)
  add_library(apec_m_solver MODULE ${SOURCES})
  set_target_properties(apec_m_solver PROPERTIES
    PREFIX ""
    OUTPUT_NAME "apec_m_solver"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../python/apec_m"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "$ORIGIN/../../../;$ORIGIN/../lib"
  )

  target_include_directories(apec_m_solver
    PRIVATE ${APEC_M_INCLUDE_DIRS}
  )

  target_compile_definitions(apec_m_solver PRIVATE APEC_M_ENABLE_PYBIND)

  if(PETSC_CFLAGS_OTHER)
    target_compile_options(apec_m_solver PRIVATE ${PETSC_CFLAGS_OTHER})
  endif()

  target_link_libraries(apec_m_solver
    PRIVATE
      pybind11::module
      ${APEC_M_LINK_LIBS}
  )

  if(APEC_M_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(apec_m_solver PRIVATE -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(apec_m_solver PRIVATE -fsanitize=address,undefined)
  endif()

  if(APEC_M_ENABLE_IMPORT_TEST)
    enable_testing()
    add_test(
      NAME apec_m_solver_import_smoke
      COMMAND ${Python3_EXECUTABLE} -c "import sys; sys.path.insert(0, r'${CMAKE_SOURCE_DIR}/../python'); import apec_m.apec_m_solver as mod; print(mod.__doc__)"
    )
  endif()

  if(APEC_M_ENABLE_VALGRIND_SMOKE)
    find_program(VALGRIND_EXECUTABLE valgrind)
    if(VALGRIND_EXECUTABLE)
      add_custom_target(valgrind_solver_smoke
        COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_SOURCE_DIR}/../python
          ${VALGRIND_EXECUTABLE} --error-exitcode=1 --leak-check=full
          ${Python3_EXECUTABLE} -c "import apec_m.apec_m_solver"
        DEPENDS apec_m_solver
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running valgrind smoke check for apec_m_solver import")
    else()
      message(WARNING "APEC_M_ENABLE_VALGRIND_SMOKE is ON but valgrind was not found.")
    endif()
  endif()

  vtk_module_autoinit(
    TARGETS apec_m_solver
    MODULES
      VTK::CommonCore
      VTK::CommonDataModel
      VTK::CommonExecutionModel
      VTK::FiltersCore
      VTK::IOXML
      VTK::WrappingPythonCore
  )
endif()
