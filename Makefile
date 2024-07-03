.PHONY: all clean lib src tests

all: lib src tests

lib:
# Command to build the library
# Example: $(MAKE) -C lib

src:
# Command to build the source (executable)
	$(MAKE) -C src

tests:
# Command to build and run tests
# Example: $(MAKE) -C tests

clean:
# Clean commands for each subdirectory
	$(MAKE) -C src clean
# Example: $(MAKE) -C lib clean
# Example: $(MAKE) -C tests clean