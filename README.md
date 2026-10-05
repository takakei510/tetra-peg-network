# tetra-peg-network

Research code for Tetra-PEG type AB-SAW gelation and molecular network analysis.

## Current status

The foundation and a minimal single-molecule demo are implemented: four sequential SAW arms, shared site ownership, full rollback on trapping, trajectory CSV, and Python 3D visualization. Growth currently chooses uniformly among unoccupied nearest neighbors with open boundaries. These are provisional conditions to confirm with the supervisor. Multi-molecule placement, AB bonding, and graph measurements are not implemented yet.

For CSV generation and plotting instructions, see [molecule visualization](docs/molecule-visualization.md). The Python script is `scripts/visualization/plot_molecule.py`; Matplotlib is needed only for visualization.

## Build and run

Requires a C11 compiler and Make.

```sh
make
make test
./build/tetra-peg-network configs/default.cfg
```

The configuration requires `dim` (2 or 3), `L` (at least 2), and an explicit unsigned 64-bit `seed`. Duplicate or unknown keys are rejected. The sample config uses a 3D cubic lattice with **open boundaries**. The command checks initialization and prints the chosen parameters; it does not run a physical simulation.

## Modules

- `config`: strict minimal configuration reader.
- `rng`: explicit SplitMix64 streams, rejection sampling for unbiased bounded integers, deterministic trial seeds. Record the master seed and trial number for repeatability.
- `lattice`: dense 2D/3D indexing, nearest neighbors, and ownership array (`-1` means empty; nonnegative values are future molecule IDs). Open boundaries only at this stage.
- `union_find`: connected components over molecule IDs; all allocated IDs initially represent active singleton nodes.

Lattice storage costs O(L^dim) time and memory. Union-find uses O(N) memory and amortized O(alpha(N)) per union/find. For large sweeps, assess memory before selecting `L` and molecule count.

The old `percolation-model` repository was consulted for module boundaries. Its site-occupation, single-walk statistics, and active-site union-find assumptions were deliberately not copied. This code has no dependency on that repository.
