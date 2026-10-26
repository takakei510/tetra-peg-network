# tetra-peg-network

Tetra-PEG型AB-SAWモデルのゲル化と、分子ネットワーク構造を調べる研究用コードです。
分子配置・末端間のA–B結合・連結成分の集計をCで実行し、CSVをPythonまたはブラウザで可視化します。

## 現在できること

- 3次元立方格子上に、4本の自己回避歩行（SAW）の腕を持つ分子を順次生成する。
- 設定ファイルで系サイズ、腕長、A型確率、目標分子数または目標占有率、seedを指定する。
- 行き止まりになった分子の配置を完全に戻し、既存分子を保持して再試行する。
- マンハッタン距離1〜2の末端候補からA–B結合を生成する。1末端1結合とし、同じ分子対の多重結合を許す。
- Union-Find（UF）で、孤立分子を含むクラスタ所属・サイズ・最大クラスタ比を集計する。
- 条件、実績、CSV、ログ、ソースの控えを実行ごとのフォルダに保存する。
- 分子形状・分子グラフ・最大コンポーネントを描画し、枝リストを抽出する。
- [ブラウザ用3Dビューア](web/viewer/)でCSVを読み込み、モデルを回転・拡大・平行移動する。

現在の生成条件は**開境界・空いている最近接点から一様に選択・4本の腕を順番に成長**です。
これは現時点のモデル仕様であり、平衡SAWの一様サンプリングや、分子の運動・反応ダイナミクスを再現するものとは区別します。
周期境界、系サイズ・占有率を変えた本格的な統計実験、一般的な網目解析、電気伝導は未実装です。
測定定義と今後の計画は[設計書](docs/ab-saw-design.md)を参照してください。

## 準備と動作確認

C処理にはC11対応コンパイラとMake、実行管理にはPython 3が必要です。WindowsではWSL上での実行を想定しています。

```bash
make
make test
```

PythonでPNGを描画するときだけ、Matplotlibを追加します。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r scripts/visualization/requirements.txt
```

Matplotlibを使わずCSVと派生グラフCSV/JSONだけ保存する場合は、実行時に`--no-plot`を指定します。
この場合もPythonによるCSV読込み・整合性検査・グラフ抽出は実行します。

## 多分子モデルを生成する

```bash
python3 scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg
```

この例は`L=16`、腕長5、目標占有率0.1、A型確率0.5、seed=3です。
1分子の占有点数は`1 + 4 × arm_length`なので、目標分子数は19に切り下げられます。

設定例：

```ini
dim = 3
L = 16
seed = 3
target_phi = 0.1
arm_length = 5
a_probability = 0.5
```

`target_phi`と`target_molecules`は**どちらか一方**を指定します。目標占有率、分子数に換算した占有率、実際の占有率は区別して記録します。
設定の詳細は[ネットワーク設定](docs/network-settings.md)、生成規則は[多分子生成](docs/multi-molecule-network.md)を参照してください。

配置試行上限は全分子合計で既定10,000回です。大きな系や高占有率では適切に増やします。

```bash
python3 scripts/run_demo.py --mode network \
  --config configs/network/preview_phi010_L16_seed3.cfg \
  --max-placement-attempts 100000 --no-plot
```

上限に達した場合は`status=incomplete`、終了コード3とし、完成分子の軌跡と配置集計を残します。
未達試料では結合・UF・グラフ解析を実行しません。正常完了した目標占有率の試料とは混ぜないでください。

## 結果の保存先

各実行は`data/runs/`の新しいフォルダに保存され、既存の実行結果を上書きしません。

| ファイル | 内容 |
|---|---|
| `config.cfg` / `run.json` / `metadata.log` | 入力条件、実行状態、seed、コード情報など |
| `trajectories.csv` | 分子中心・各腕の格子点と方向 |
| `placement-summary.csv` | 目標と実績、型別分子数、配置試行・棄却回数 |
| `candidates.csv` | 距離条件を満たすA–B末端候補。成立結合とは別 |
| `bonds.csv` | 実際に成立した末端間結合 |
| `cluster-members.csv` | 各分子の所属クラスタと次数 |
| `cluster-summary.csv` | クラスタごとの分子数と結合本数 |
| `network-summary.csv` | 最大・第2位クラスタ、孤立分子数、最大クラスタ比など |
| `molecules.png` | 全体の実空間図。描画を有効にした場合のみ |
| `source.zip` / `build.log` | 再現用ソースの控えとビルドログ |

`network-summary.csv`の`largest_fraction`は**最大クラスタ分子数 / 実生成分子数**です。
先生の資料の`S1 = (1 + 4l) × 最大クラスタ分子数 / L³`とは正規化が異なります。
実占有率を`phi`とすると`S1 = phi × largest_fraction`です。

CSVのID対応と全列の意味は[設計書9〜21節](docs/ab-saw-design.md)、保存方針は[実行管理](docs/run-management.md)を参照してください。

## 分子グラフと最大コンポーネント

`network`または`clusters`モードの正常完了時は`graph/`も保存します。

| ファイル | 内容 |
|---|---|
| `network.png` | 分子を頂点、成立結合を枝とする全体図 |
| `largest-network.png` | 代表最大コンポーネントのグラフ図 |
| `largest-molecules.png` | 代表最大コンポーネントの実空間形状 |
| `molecular-edges.csv` | 全成立枝。元の結合・分子・末端IDを保持 |
| `largest-members.csv` / `largest-edges.csv` | 代表最大成分の所属分子と成立枝 |
| `graph-analysis.json` | 同率最大ID、代表選択、集計値、描画状態 |

同率最大が複数ある場合、最小`cluster_id`の1成分を代表として抽出します。
多重結合は別々の枝として保持します。グラフ図のA/B二列配置は表示用で、距離は物理的距離ではありません。
`--no-plot`ではPNGを作成しません。

保存済みの正常完了結果も再生成せず解析できます。

```bash
python3 scripts/visualization/analyze_graph.py data/runs/<実行フォルダ>
```

既存の`graph/`は上書きしません。別の出力先には`--output-dir`を指定します。
詳細は[設計書23〜24節](docs/ab-saw-design.md)を参照してください。

## CSVを読み込む3Dビューア

[Web版ビューアを開く](https://ab-saw-model-viewer-kei.equalbox128.chatgpt.site)（所有者のみアクセス可能）。操作と制限は[ビューアの説明](web/viewer/README.md)を参照してください。
ローカルでも、リポジトリのルートで次を実行できます。

```bash
python3 -m http.server 8000 --bind 127.0.0.1 --directory web/viewer
```

ブラウザで`http://localhost:8000`を開き、同じ実行フォルダの`trajectories.csv`と、必要に応じて`bonds.csv`を同時に選択します。
CSVはブラウザ内で処理します。通常は表示を優先し、詳細検証は別ボタンで実行します。処理段階の進捗表示・中止・軽量な表示専用JSONの保存と再読込みにも対応しています。ドラッグ回転、ズーム、平行移動、A/B表示切替、最大成分の表示、成分別着色、PNG保存に対応します。
**軌跡だけでも表示できますが、結合・クラスタの判定には`bonds.csv`が必要です。**
大規模データには読込み上限と分子単位の描画間引きがあります。読込み上限内でも、GPU・ブラウザによって動作速度は異なります。

## 段階別のデモ

従来の固定条件デモも残しています。多分子生成には`network`モードを使用してください。

| モード | 内容 | 設定例 |
|---|---|---|
| `single` | 1分子の生成 | `configs/demos/single_L32_seed12345.cfg` |
| `two` | 固定A/Bの2分子配置 | `configs/demos/two_L8_seed12345.cfg` |
| `candidates` | 末端候補の探索 | `configs/demos/candidates_L8_seed3.cfg` |
| `bonds` | 末端間の結合 | `configs/demos/bonds_L8_seed3.cfg` |
| `clusters` | 分子クラスタ集計 | `configs/demos/clusters_L8_seed3.cfg` |

```bash
python3 scripts/run_demo.py --mode clusters --config configs/demos/clusters_L8_seed3.cfg
```

これらの旧デモは腕長や分子数に固定条件を含みます。基盤確認用の`./build/tetra-peg-network configs/default.cfg`は初期化チェックのみで、分子生成は行いません。
詳細：[1分子](docs/molecule-visualization.md) / [2分子](docs/two-molecules.md) / [末端探索](docs/endpoint-search.md) / [結合](docs/endpoint-bonding.md) / [クラスタ](docs/molecular-clusters.md) / [設定一覧](configs/README.md)

## モジュールと計算資源

| モジュール | 役割 |
|---|---|
| `config` / `network_config` | 未知・重複キーなどを検査する設定読込み |
| `rng` | 明示的なseedを持つSplitMix64と、偏りを避ける整数抽選 |
| `lattice` | 格子の添字変換・最近接点・占有分子IDの管理 |
| `molecule` / `placement` | 4腕分子の試行生成、失敗復元、目標数までの反復 |
| `endpoint_search` / `bond` | 末端の空間索引、候補探索、結合の確定 |
| `union_find` / `cluster` | 分子レベルの連結成分と集計 |

格子と末端索引は密な配列を用いるため、3次元では`O(L³)`のメモリが必要です。
固定占有率・腕長では、軌跡点数と保存量も概ね`L³`に比例します。
Python後処理・描画は別にメモリを使います。PC・サーバの使い分けは、今後の性能測定で判断します。

旧`percolation-model`はモジュール分割の参考にしましたが、本リポジトリはそれに依存しません。

## Lを変えた計算時間・メモリ・保存容量の測定

```bash
python3 scripts/benchmark.py --config configs/benchmarks/pc-pilot.json
```

初期設定はL=16,32,64,128、占有率0.1と0.3、各3seedの24試料です。
C処理・Python後処理・描画を独立に測り、未達・省略・資源上限での停止を区別します。
CSVはハッシュ照合済みのgzipで保持し、必要な試料だけ復元できます。
詳しい条件・表の読み方・PC/サーバの判断・容量の上限は[性能測定ガイド](docs/benchmark.md)を参照してください。


## 保存済み試料のH₂・平均クラスタサイズ解析

Gitで更新してから、既存の105試料を後処理できます。実験CSVの再生成やZIPでの解析コード配置は不要です。

```bash
git pull --ff-only
python3 scripts/analysis/analyze_h2_clusters.py
```

`data/`以下で元の105試料を含む保存先が1か所なら自動選択します。見つからない場合や複数ある場合は、保存先を指定します。

```bash
python3 scripts/analysis/analyze_h2_clusters.py \
  --suite data/report_figures/20261019T150945Z_9feb1c/simulation
```

結果はGit管理外の`data/analysis/h2_<UTC日時>_<識別子>/`に保存します。
最後の`UPLOAD THIS FILE:`に表示する`upload_h2_cluster_results.zip`を送れば、報告用の図へ反映できます。
Python標準ライブラリでCSV・JSON・ZIPを作り、Matplotlibがあれば図も保存します。図を省略する場合は`--no-plots`を指定します。

今回のH₂は「ちょうど2本で結ばれた分子ペア数」です。三重・四重結合ペアを含む全閉路数とは区別します。
平均クラスタサイズは、算術平均に加え、最大クラスタ1個を除いた分子数重み付き平均を求めます。
[実行手順と測定定義](docs/h2-cluster-analysis.md)、[解析条件](configs/analysis/h2-cluster-report.json)を参照してください。

```bash
make test-analysis
```

## 占有率0.4〜0.5の詳細計算

既存のΦ=0.40/0.50を再利用し、新たに0.41〜0.49を0.01刻みで調べます。
まずL=64、Φ=0.42/0.44/0.46/0.48、各5seedの20試料を実行します。
今回は保存CSVの再集計に加え、新しい分子ネットワークの生成も行います。

```bash
python3 scripts/run_phi_detail.py pilot
```

正常完了すると、H₂・平均クラスタサイズの集計と送付用ZIPまで作成します。
保存先は`data/refinement/phi040_050/`で、最後の`UPLOAD THIS FILE:`にZIPのパスを表示します。
20試料の結果を確認した後、全3サイズの詳細計算へ広げます。

```bash
python3 scripts/run_phi_detail.py full
```

完了済みのpilot20試料を再利用し、残り115試料を生成します。合計135試料で重複しません。
途中停止したグループは自動でやり直さず、元の出力を保って停止します。
[詳細計算の手順・条件・再利用の範囲](docs/phi-detail-study.md)を参照してください。

```bash
python3 scripts/run_phi_detail.py plan
make test-detail
```
