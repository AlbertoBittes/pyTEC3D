cmake_minimum_required(VERSION 3.10)
project(ApecM VERSION 1.0)

# Set the C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -g -fopenmp")

# Set the MPI C++ compiler
set(CMAKE_CXX_COMPILER mpicxx)

# Include directories
# include_directories(${LIBMESH_DIR}/include)
# include_directories(${PETSC_DIR}/include)
# include_directories(/usr/include/vtk-9.1)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/dependencies)

# Source files
file(GLOB SOURCES "*.cpp")

# Library directories
# link_directories(${LIBMESH_DIR}/lib)
link_directories(/usr/lib/x86_64-linux-gnu)
# link_directories(${CMAKE_SOURCE_DIR}/lib)

# Libraries to link against
set(LIBS
    mesh_dbg
    timpi_dbg
    pthread
    gomp
    vtkCommonCore-9.1
    vtkCommonDataModel-9.1
    vtkFiltersCore-9.1
    vtkIOXML-9.1
    vtkCommonExecutionModel-9.1
)

# Add executable
add_executable(plasma_arc ${SOURCES})

# Set the output directory for the executable
set_target_properties(plasma_arc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)

# Link libraries
target_link_libraries(plasma_arc ${LIBS})

# Set RPATH
set_target_properties(plasma_arc PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "$ORIGIN/../lib"
)

# Add a custom target for cleaning up build files
add_custom_target(clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/../bin/plasma_arc
)
