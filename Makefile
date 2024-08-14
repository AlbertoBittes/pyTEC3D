.PHONY: all clean lib src tests

# Default to using the standard Makefile in src
SRC_MAKEFILE ?= Makefile

all: lib src tests

lib:
# Command to build the library
# Example: $(MAKE) -C lib

src:
# Command to build the source (executable)
	$(MAKE) -C src -f $(SRC_MAKEFILE)

tests:
# Command to build and run tests
# Example: $(MAKE) -C tests

clean:
# Clean commands for each subdirectory
	$(MAKE) -C src -f $(SRC_MAKEFILE) clean
# Example: $(MAKE) -C lib clean
# Example: $(MAKE) -C tests clean