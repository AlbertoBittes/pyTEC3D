# Compiler settings - change this to your MPI C++ compiler
CXX=/home/bitteal/petsc/arch-linux-c-debug/bin/mpicxx

# Include path for libMesh
INCLUDES=-I/home/bitteal/libmesh_inst/include -I/home/bitteal/petsc/include

# Compiler flags
CXXFLAGS=-std=c++17 -Wall -O2 $(INCLUDES) -pthread

# Linker flags, add your project-specific libraries not related to MPI
LDFLAGS=-L/home/bitteal/libmesh_inst/lib -lmesh_dbg -ltimpi_dbg -pthread

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
