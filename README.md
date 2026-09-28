# tetra-peg-network

Research code for Tetra-PEG type AB-SAW gelation and molecular network analysis.

## Current status

This first stage provides only the simulation foundation. Molecular placement, SAW growth, AB bonding, graph measurements, and CSV output are **not implemented yet**. The physical definitions and boundary conditions for those stages must be fixed against the supervisor's model before implementation.

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
