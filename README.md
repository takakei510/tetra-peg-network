# tetra-peg-network

Research code for Tetra-PEG type AB-SAW gelation and molecular network analysis.

## Current status

The foundation, a single-molecule demo, and a fixed two-molecule placement demo are implemented: four sequential SAW arms, shared site ownership, full rollback on trapping, trajectory CSV, and Python 3D visualization. Growth currently chooses uniformly among unoccupied nearest neighbors with open boundaries. These are provisional conditions to confirm with the supervisor. Manhattan-radius-2 endpoint candidate search is now implemented (stage 3-A). Endpoint-order random AB bonding is implemented (stage 3-B), with at most one bond per endpoint and optional molecular multibonds. Molecule-level connected-component summaries are implemented with the existing union-find foundation (stage 3-C). General configurable multi-molecule placement and broader graph measurements are not implemented yet.

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

All endpoint IDs are shuffled, then each free endpoint chooses uniformly from its currently free candidate partners. The demo allows separate arms to connect the same molecule pair more than once. It saves both `candidates.csv` and `bonds.csv` and draws selected bonds as solid green lines. See [endpoint bonding](docs/endpoint-bonding.md) for the rule, data structures, validation, and limitations. Use the clusters mode below for union-find summaries.

## Molecular clusters (stage 3-C)

```sh
python scripts/run_demo.py --mode clusters --config configs/demos/clusters_L8_seed3.cfg
```

This keeps the same placement and bonds and adds per-molecule memberships, per-cluster sizes/bond counts, and a network summary. Isolated molecules are included. Multiple bonds between a molecule pair retain their bond count while union-find joins the molecules only once. See [molecular clusters](docs/molecular-clusters.md) for why UF is used, how it differs from bond storage, complexity, and BFS validation. This is a small demonstration, not a gelation measurement.

## Network configuration preparation (next step 1)

```sh
make build/check_network_config
./build/check_network_config configs/network/preview_N20_L16_seed3.cfg
```

A separate `NetworkConfig` adds arm length, A-type probability, and exactly one of `target_molecules` or `target_phi`. Occupancy input is converted to an integer target molecule count by flooring; requested and resolved occupancy are reported separately. Use `configs/network/preview_phi010_L16_seed3.cfg` for the occupancy-input example. This command only validates the configuration and reports `generation_run=0`. The network mode below connects these settings to multi-molecule generation. Existing three-key demo configs and fixed demo conditions remain unchanged. See [network settings](docs/network-settings.md).

## Configurable multi-molecule network (next step 2)

```sh
python scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg
```

This generates the resolved target count with random A/B types, then reuses the endpoint search, bonding and union-find modules and saves a Python figure. The example resolves to 19 molecules. `placement-summary.csv` separates requested/resolved occupancy from actual placement, counts occupied-center and trapped rejections, and records completion. The default attempt budget is 10000 total trials; override with `--max-placement-attempts N` (network mode only). Budget exhaustion saves completed trajectories with `status=incomplete` and exit code 3, skips bonding/UF, and plots partial trajectories only when nonempty. See [multi-molecule generation](docs/multi-molecule-network.md).


### 分子グラフと最大コンポーネント

`--mode network`または`--mode clusters`の正常完了時は、既存の結果に加えて`graph/`へ保存します。

| ファイル | 内容 |
|---|---|
| `network.png` | 分子を頂点、成立結合を枝とする全体図。A/B二列の表示用配置 |
| `largest-network.png` | 代表最大コンポーネントのグラフ図 |
| `largest-molecules.png` | 代表最大コンポーネントの実空間形状 |
| `molecular-edges.csv` | 全成立枝。元の結合・分子・末端IDを保持 |
| `largest-members.csv` / `largest-edges.csv` | 代表最大成分の所属分子と成立枝 |
| `graph-analysis.json` | 代表選択・同率最大ID・集計・正規化・描画状態 |

同率最大が複数ある場合は最小cluster_idの1成分を代表として抽出し、同率の一覧をJSONに残します。
多重結合は別の曲線で描画します。線の交差は接続を表しません。グラフ図の距離は物理的距離ではありません。
目標未達の試料には結合・グラフ解析を実行しません。`--no-plot`でも派生CSV/JSONは保存し、Matplotlibは不要です。

保存済みのcompleted実行も再生成なしで後処理できます。

```bash
python scripts/visualization/analyze_graph.py data/runs/<実行フォルダ>
```

既存の`graph/`は上書きしません。再解析には`--output-dir`で別の出力先を指定してください。
定義・列・処理状態の詳細は[設計書23〜24節](docs/ab-saw-design.md)を参照してください。
