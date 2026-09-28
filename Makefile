CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS ?= -Iinclude
CORE = src/config.c src/rng.c src/lattice.c src/union_find.c
.PHONY: all test clean
all: build/tetra-peg-network
build:
	mkdir -p build
build/tetra-peg-network: src/main.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_foundation: tests/test_foundation.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
test: build/test_foundation
	./build/test_foundation
clean:
	rm -rf build
