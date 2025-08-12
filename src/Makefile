# Compiler settings - change this to your MPI C++ compiler
CXX=mpicxx

# Include path for libMesh
INCLUDES=-I../include -I/$(LIBMESH_DIR)/include -I/$(PETSC_DIR)/include -I/usr/include/vtk-9.1

# Compiler flags (add -MMD -MP for dependency generation)
CXXFLAGS=-std=c++17 -Wall $(INCLUDES) -pthread -O0 -g -fopenmp -MMD -MP

# Linker flags, add your project-specific libraries not related to MPI
LDFLAGS=-L/$(LIBMESH_DIR)/lib -L/usr/lib/x86_64-linux-gnu -lmesh_dbg -ltimpi_dbg -lpetsc -pthread -fopenmp \
-lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersGeneral-9.1 -lvtkIOXML-9.1 \
-lvtkCommonExecutionModel-9.1 -lvtkIOParallelXML-9.1 -lvtkIOCore-9.1 -lvtkIOGeometry-9.1 \
-lvtkIOLegacy-9.1 -lvtkParallelMPI-9.1 

# The target executable name
BIN_DIR=../bin
TARGET=$(BIN_DIR)/apec_m

SOURCES=$(wildcard *.cpp)
ELEMENT_BUILDERS_SOURCES := $(wildcard element_builders/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o) $(ELEMENT_BUILDERS_SOURCES:.cpp=.o)
DEPS=$(SOURCES:.cpp=.d) $(ELEMENT_BUILDERS_SOURCES:.cpp=.d)

%.o: %.cpp
    $(CXX) $(CXXFLAGS) -c $< -o $@

element_builders/%.o: element_builders/%.cpp
    $(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
    $(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

-include $(DEPS)

clean:
    rm -f *.o *.d element_builders/*.o element_builders/*.d ../apec_m