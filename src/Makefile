# Default to local libMesh install; environment overrides still work
LIBMESH_DIR ?= /home/bitteal/libmesh_inst
LIBMESH_CONFIG ?= $(LIBMESH_DIR)/bin/libmesh-config
METHOD ?= dbg

LIBMESH_CXX := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --cxx)
LIBMESH_INCLUDE := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --include)
LIBMESH_CXXFLAGS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --cxxflags)
LIBMESH_LIBS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --libs)
VTK_LIBS ?= -lvtkFiltersGeneral-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1

# Compiler settings - default to libMesh compiler so flags/libs stay consistent
CXX=$(LIBMESH_CXX)

# Include path for libMesh/PETSc pulled from libmesh-config
INCLUDES=-I../include $(LIBMESH_INCLUDE)

# Compiler flags (add -MMD -MP for dependency generation)
CXXFLAGS=-std=c++17 -Wall $(INCLUDES) $(LIBMESH_CXXFLAGS) -pthread -MMD -MP

# Linker flags, add your project-specific libraries not related to MPI
LDFLAGS=$(LIBMESH_LIBS) $(VTK_LIBS)

# The target executable name
BIN_DIR=../bin
TARGET=$(BIN_DIR)/apec_m

SOURCES=$(wildcard *.cpp)
ELEMENT_BUILDERS_SOURCES := $(wildcard element_builders/*.cpp)
ELEMENT_TYPES_SOURCES := $(wildcard element_types/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o) $(ELEMENT_BUILDERS_SOURCES:.cpp=.o) $(ELEMENT_TYPES_SOURCES:.cpp=.o)
DEPS=$(SOURCES:.cpp=.d) $(ELEMENT_BUILDERS_SOURCES:.cpp=.d) $(ELEMENT_TYPES_SOURCES:.cpp=.d)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

element_builders/%.o: element_builders/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

element_types/%.o: element_types/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

-include $(DEPS)

clean:
	rm -f *.o *.d element_builders/*.o element_builders/*.d element_types/*.o element_types/*.d ../apec_m
