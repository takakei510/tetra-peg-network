CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS ?= -Iinclude
CORE = src/config.c src/rng.c src/lattice.c src/union_find.c src/molecule.c src/endpoint_search.c src/bond.c
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
build/demo_candidates: src/demo_candidates.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_endpoint_search: tests/test_endpoint_search.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/demo_bonds: src/demo_bonds.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_bond: tests/test_bond.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
test: build/test_bond build/test_foundation build/test_molecule build/test_endpoint_search
	./build/test_bond
	./build/test_foundation
	./build/test_molecule
	./build/test_endpoint_search
clean:
	rm -rf build
