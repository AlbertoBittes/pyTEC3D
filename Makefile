.PHONY: all clean lib src tests module module-dbg module-opt

# Default to using the standard Makefile in src
SRC_MAKEFILE ?= Makefile

all: lib src tests

lib:
# Command to build the library
# Example: $(MAKE) -C lib

src:
# Command to build the source (executable)
	$(MAKE) -C src -f $(SRC_MAKEFILE)

module:
# Command to build only the Python pybind11 module
	$(MAKE) -C src -f $(SRC_MAKEFILE) module

module-dbg:
# Build pybind11 module with debug flags from libMesh config
	$(MAKE) -C src -f $(SRC_MAKEFILE) METHOD=dbg module

module-opt:
# Build pybind11 module with optimized flags from libMesh config
	$(MAKE) -C src -f $(SRC_MAKEFILE) METHOD=opt module

tests:
# Command to build and run tests
# Example: $(MAKE) -C tests

clean:
# Clean commands for each subdirectory
	$(MAKE) -C src -f $(SRC_MAKEFILE) clean
# Example: $(MAKE) -C lib clean
# Example: $(MAKE) -C tests clean
