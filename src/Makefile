# Compiler settings - change this to your MPI C++ compiler
CXX=mpicxx

# Include path for libMesh
INCLUDES=-I/$(LIBMESH_DIR)/include -I/$(PETSC_DIR)/include -I/usr/include/vtk-9.1

# Compiler flags
CXXFLAGS=-std=c++17 -Wall $(INCLUDES) -pthread -g -fopenmp

# Linker flags, add your project-specific libraries not related to MPI
LDFLAGS=-L/$(LIBMESH_DIR)/lib -L/usr/lib/x86_64-linux-gnu -lmesh_dbg -ltimpi_dbg -lpetsc -pthread -fopenmp \
-lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkFiltersCore-9.1 -lvtkIOXML-9.1 \
-lvtkCommonExecutionModel-9.1 -lvtkIOParallelXML-9.1 -lvtkIOCore-9.1 -lvtkIOGeometry-9.1 \
-lvtkIOLegacy-9.1 -lvtkParallelCore-9.1 -lvtkRenderingCore-9.1 -lvtkRenderingFreeType-9.1 \
-lvtkRenderingOpenGL2-9.1 -lvtkFiltersExtraction-9.1 -lvtkParallelMPI-9.1 -lvtkFiltersParallelGeometry-9.1 \
-lvtkFiltersGeneral-9.1

# The target executable name
BIN_DIR=../bin
TARGET=$(BIN_DIR)/plasma_arc

# Source files
SOURCES=$(wildcard *.cpp)

# Object files
OBJECTS=$(SOURCES:.cpp=.o)

# Rule to make the object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link the program
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f *.o ../plasma_arc