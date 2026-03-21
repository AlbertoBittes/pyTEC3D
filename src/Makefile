# Default to local libMesh install; environment overrides still work
LIBMESH_DIR ?= /home/bitteal/libmesh_inst
LIBMESH_CONFIG ?= $(LIBMESH_DIR)/bin/libmesh-config
METHOD ?= dbg
PYTHON ?= python3

LIBMESH_CXX := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --cxx)
LIBMESH_INCLUDE := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --include)
LIBMESH_CXXFLAGS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --cxxflags)
LIBMESH_LIBS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --libs)
PYBIND11_INCLUDES ?= $(shell $(PYTHON) -m pybind11 --includes 2>/dev/null)
ifeq ($(strip $(PYBIND11_INCLUDES)),)
PYBIND11_INCLUDES := $(shell pybind11-config --includes 2>/dev/null)
endif
PYBIND11_INCLUDE_DIR ?=
PYTHON_INCLUDES ?= $(shell $(PYTHON)-config --includes 2>/dev/null)
PYTHON_EXT_SUFFIX ?= $(shell $(PYTHON)-config --extension-suffix 2>/dev/null)
PYTHON_LDFLAGS ?= $(shell $(PYTHON)-config --ldflags 2>/dev/null)
VTK_LIBS ?= -lvtkFiltersGeneral-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1

# Compiler settings - default to libMesh compiler so flags/libs stay consistent
CXX=$(LIBMESH_CXX)

# Include path for libMesh/PETSc pulled from libmesh-config
INCLUDES=-I../include $(LIBMESH_INCLUDE) $(PYBIND11_INCLUDES) $(if $(PYBIND11_INCLUDE_DIR),-I$(PYBIND11_INCLUDE_DIR),) $(PYTHON_INCLUDES)

# Compiler flags (add -MMD -MP for dependency generation)
CXXFLAGS=-std=c++17 -Wall $(INCLUDES) $(LIBMESH_CXXFLAGS) -fPIC -pthread -MMD -MP

# Linker flags, add your project-specific libraries not related to MPI
LDFLAGS=$(LIBMESH_LIBS) $(VTK_LIBS) $(PYTHON_LDFLAGS)

# Python module output
API_DIR=../api
MODULE_NAME=apec_m_solver
TARGET=$(API_DIR)/$(MODULE_NAME)$(PYTHON_EXT_SUFFIX)

SOURCES=$(wildcard *.cpp)
ELEMENT_BUILDERS_SOURCES := $(wildcard element_builders/*.cpp)
ELEMENT_TYPES_SOURCES := $(wildcard element_types/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o) $(ELEMENT_BUILDERS_SOURCES:.cpp=.o) $(ELEMENT_TYPES_SOURCES:.cpp=.o)
DEPS=$(SOURCES:.cpp=.d) $(ELEMENT_BUILDERS_SOURCES:.cpp=.d) $(ELEMENT_TYPES_SOURCES:.cpp=.d)

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

element_builders/%.o: element_builders/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

element_types/%.o: element_types/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	mkdir -p $(API_DIR)
	$(CXX) -shared $(OBJECTS) -o $(TARGET) $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f *.o *.d element_builders/*.o element_builders/*.d element_types/*.o element_types/*.d \
	      $(API_DIR)/$(MODULE_NAME)*.so $(API_DIR)/$(MODULE_NAME)*.pyd $(API_DIR)/$(MODULE_NAME)*.dylib
