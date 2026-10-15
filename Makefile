CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS ?= -Iinclude
CORE = src/config.c src/rng.c src/lattice.c src/union_find.c src/molecule.c src/endpoint_search.c src/bond.c src/cluster.c src/network_config.c src/placement.c
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
build/demo_clusters: src/demo_clusters.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_cluster: tests/test_cluster.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/check_network_config: src/check_network_config.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_network_config: tests/test_network_config.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/demo_network: src/demo_network.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
build/test_placement: tests/test_placement.c $(CORE) | build
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ -o $@
test: build/demo_network build/test_placement build/test_network_config build/test_cluster build/test_bond build/test_foundation build/test_molecule build/test_endpoint_search
	./build/test_placement
	./build/test_network_config
	./build/test_cluster
	./build/test_bond
	./build/test_foundation
	./build/test_molecule
	./build/test_endpoint_search
	python3 tests/test_network_run.py
	python3 tests/test_graph_analysis.py
clean:
	rm -rf build

.PHONY: test-viewer
# ビューアを変更するときのみNode.jsで実行する。C処理の依存は増やさない。
test-viewer:
	node tests/test_viewer.mjs
	node tests/test_viewer_ui.mjs
	node tests/test_viewer_ui.mjs --canvas

.PHONY: test-benchmark
test-benchmark:
	python3 tests/test_benchmark.py
