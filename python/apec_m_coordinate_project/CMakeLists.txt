cmake_minimum_required(VERSION 3.18)
project(apec_m_coordinate VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set(APEC_M_REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")

add_library(apec_m_coordinate_system MODULE
  "${APEC_M_REPO_ROOT}/src/coordinate/coordinate_system.cpp"
  "${APEC_M_REPO_ROOT}/src/coordinate/coordinate_system_bindings.cpp"
  "${APEC_M_REPO_ROOT}/src/coordinate/coordinate_system_module.cpp"
)

set_target_properties(apec_m_coordinate_system PROPERTIES
  PREFIX ""
  OUTPUT_NAME "apec_m_coordinate_system"
)

target_include_directories(apec_m_coordinate_system PRIVATE
  "${APEC_M_REPO_ROOT}/include"
)

target_link_libraries(apec_m_coordinate_system PRIVATE pybind11::module)

install(TARGETS apec_m_coordinate_system LIBRARY DESTINATION apec_m_coordinate)
install(FILES
  "${APEC_M_REPO_ROOT}/python/apec_m_coordinate/__init__.py"
  "${APEC_M_REPO_ROOT}/python/apec_m_coordinate/coordinate_system.py"
  DESTINATION apec_m_coordinate
)
