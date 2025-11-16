cmake_minimum_required(VERSION 3.16)
project(ApecM VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)
find_package(VTK REQUIRED
  COMPONENTS
    CommonCore
    CommonDataModel
    CommonExecutionModel
    FiltersCore
    IOXML
    IOParallelXML
    ParallelMPI
)

set(LIBMESH_DIR "${LIBMESH_DIR}" CACHE PATH "Path to the libMesh installation.")
if(NOT LIBMESH_DIR AND DEFINED ENV{LIBMESH_DIR})
  set(LIBMESH_DIR "$ENV{LIBMESH_DIR}")
endif()
if(NOT LIBMESH_DIR)
  message(FATAL_ERROR "LIBMESH_DIR must be provided (e.g., -DLIBMESH_DIR=/path/to/libmesh).")
endif()

find_library(LIBMESH_LIB
  NAMES mesh_opt mesh_dbg
  HINTS "${LIBMESH_DIR}/lib" "${LIBMESH_DIR}/lib64"
  REQUIRED
)
find_library(TIMPI_LIB
  NAMES timpi_opt timpi_dbg
  HINTS "${LIBMESH_DIR}/lib" "${LIBMESH_DIR}/lib64"
  REQUIRED
)

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(apec_m ${SOURCES})
set_target_properties(apec_m PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin"
  BUILD_WITH_INSTALL_RPATH TRUE
  INSTALL_RPATH "$ORIGIN/../lib"
)

target_include_directories(apec_m
  PRIVATE
    "${LIBMESH_DIR}/include"
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/dependencies
    ${PETSC_INCLUDES}
)

target_compile_definitions(apec_m PRIVATE ${PETSC_DEFINITIONS})

target_link_libraries(apec_m
  PRIVATE
    ${LIBMESH_LIB}
    ${TIMPI_LIB}
    ${PETSC_LIBRARIES}
    MPI::MPI_CXX
    Threads::Threads
    OpenMP::OpenMP_CXX
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::CommonExecutionModel
    VTK::FiltersCore
    VTK::IOXML
    VTK::IOParallelXML
    VTK::ParallelMPI
)

vtk_module_autoinit(
  TARGETS apec_m
  MODULES
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::CommonExecutionModel
    VTK::FiltersCore
    VTK::IOXML
    VTK::IOParallelXML
    VTK::ParallelMPI
)
