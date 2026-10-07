# tetra-peg-network

Research code for Tetra-PEG type AB-SAW gelation and molecular network analysis.

## Current status

The foundation, a single-molecule demo, and a fixed two-molecule placement demo are implemented: four sequential SAW arms, shared site ownership, full rollback on trapping, trajectory CSV, and Python 3D visualization. Growth currently chooses uniformly among unoccupied nearest neighbors with open boundaries. These are provisional conditions to confirm with the supervisor. Manhattan-radius-2 endpoint candidate search is now implemented (stage 3-A). Endpoint-order random AB bonding is implemented (stage 3-B), with at most one bond per endpoint and optional molecular multibonds. General configurable multi-molecule placement and graph measurements are not implemented yet.

For CSV generation and plotting instructions, see [molecule visualization](docs/molecule-visualization.md). The Python script is `scripts/visualization/plot_molecule.py`; Matplotlib is needed only for visualization.

For the two-molecule demo (one fixed A and one fixed B), see [two molecules](docs/two-molecules.md). Use `make test build/demo_two_molecules`, then its CSV with `scripts/visualization/plot_molecules.py` to preserve their positions in shared absolute coordinates.

## Build and run

To keep each demo's configuration, CSV, plot, logs, and source snapshot together without overwriting previous results, use:

```sh
python scripts/run_demo.py --mode two --config configs/demos/two_L8_seed12345.cfg
```

Results go to a new `data/runs/` folder per execution. See [run management](docs/run-management.md) and [configuration presets](configs/README.md). Add `--no-plot` for CSV-only runs. This runner uses the existing demo conditions; it does not make hardcoded arm length or molecule counts configurable.

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

## Endpoint candidate search (stage 3-A)

```sh
python scripts/run_demo.py --mode candidates --config configs/demos/candidates_L8_seed3.cfg
```

This saves `candidates.csv` alongside trajectories and draws potential pairs as dashed lines. This candidates-only mode does not select bonds. See [endpoint search](docs/endpoint-search.md) for the index design, complexity, validation, and code reuse.

## Random endpoint bonding (stage 3-B)

```sh
python scripts/run_demo.py --mode bonds --config configs/demos/bonds_L8_seed3.cfg
```

All endpoint IDs are shuffled, then each free endpoint chooses uniformly from its currently free candidate partners. The demo allows separate arms to connect the same molecule pair more than once. It saves both `candidates.csv` and `bonds.csv` and draws selected bonds as solid green lines. See [endpoint bonding](docs/endpoint-bonding.md) for the rule, data structures, validation, and limitations. Union-find cluster analysis comes next.
