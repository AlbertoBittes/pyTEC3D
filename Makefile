.PHONY: all clean module module-dbg module-opt

# Default to using the standard Makefile in src
SRC_MAKEFILE ?= Makefile

all: module

module:
	$(MAKE) -C src -f $(SRC_MAKEFILE) module

module-dbg:
	$(MAKE) -C src -f $(SRC_MAKEFILE) METHOD=dbg module

module-opt:
	$(MAKE) -C src -f $(SRC_MAKEFILE) METHOD=opt module

clean:
	$(MAKE) -C src -f $(SRC_MAKEFILE) clean
