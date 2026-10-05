CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS ?= -Iinclude
CORE = src/config.c src/rng.c src/lattice.c src/union_find.c src/molecule.c
.PHONY: all test clean
all: build/tetra-peg-network
build:
	mkdir -p build
build/tetra-peg-network: src/main.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_foundation: tests/test_foundation.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/demo_molecule: src/demo_molecule.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/demo_two_molecules: src/demo_two_molecules.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_molecule: tests/test_molecule.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
test: build/test_foundation build/test_molecule
	./build/test_foundation
	./build/test_molecule
clean:
	rm -rf build
