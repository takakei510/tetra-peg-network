# ③-C：Union-Findによる分子クラスタ集計

## 集計する対象

1分子を1頂点、BondStoreに保存した成立済み末端結合を辺として扱う。
末端IDを4で割れば所属分子IDになる。末端や格子点をUFの要素にはしない。
全登録分子をサイズ1で初期化し、孤立分子も集計に含める。
直接結合していなくても、別の分子を経由してつながれば同じクラスタ。

## UFを使う理由

各結合について2分子を`uf_union`するだけで、接続する集合をまとめられる。
全分子へクラスタ番号を付け直す処理を結合のたびに行う必要がない。
既存の経路圧縮（path halving）とサイズによる統合をそのまま利用できる。
同じクラスタの分子間の結合は統合なしを返すため、多重結合・閉路で分子数を重複計上しない。

今回の連結成分集計だけならBFS/DFSも妥当で、1回の解析はO(N+B)。
UFが必ず速いとは説明しない。今回UFを選ぶ理由は、既存基盤の再利用、
BondStoreを1回読む集計、分子隣接リストを別途作らずO(N)の追加領域で所属を管理できること。
将来、結合を増やしてサイズ変化を追う場合もUFの統合を継続できる。
現在の`cluster_build`は成立済み結合の一括集計で、逐次実験自体は未実装。
結合削除には通常のUFだけでは対応せず、再構築等が必要になる。

## 既存UFの配列と関数

`UnionFind`は`n`、`components`、`parent`、`size`を保持する。
`parent[m]`は親分子ID。`parent[root]==root`が代表。
`uf_find(m)`は親をたどって代表を求め、途中の親を短縮する。
`size[root]`は所属分子数の正本。現行実装では統合された非rootのsizeを0にする。
`uf_union(a,b)`は代表同士を比較し、小さい集合を大きい集合へ接続する。
同じ代表なら0、統合したら1、不正な分子IDなら-1。
`uf_largest_two`で最大と第2位のサイズを取得する。同率も別クラスタとして数える。
第2位が存在しないときは0。基盤の`uf_init`はn>0を要求するので、空ネットワークはcluster側で扱う。

## 新しい派生データ

| ClusterStoreのフィールド | 内容 |
|---|---|
| molecule_count / bond_count / cluster_count | 分子数 / 結合数 / クラスタ数 |
| molecule_cluster[m] | 分子IDから出力用クラスタID |
| cluster_sizes[id] | そのクラスタの分子数 |
| cluster_bonds[id] | そのクラスタ内の結合数。多重結合を含む |
| bond_degree[m] | その分子で結合した腕数。異なる隣接分子数とは別 |
| largest / second_largest | 最大 / 第2位クラスタの分子数 |
| isolated_molecules | 結合次数0の分子数 |

UFのroot番号は結合の入力順で変わる。出力クラスタIDは最小所属分子IDの順に0から付け、
同じ接続関係なら結合順を変えても同じ所属IDを得る。クラスタIDは試料内のラベルで、
結合を追加・変更した場合に以前のIDを保持する仕様ではない。
分子・結合をまとめるのは所属の派生配列であり、元の軌跡・結合配列を移動や上書きしない。
BondStoreのIDと逆引きの整合性を検査する。型・距離は前段の候補探索と結合処理の責務。

## 出力と実行

```sh
python scripts/run_demo.py --mode clusters --config configs/demos/clusters_L8_seed3.cfg
```

既存の軌跡・候補・結合CSVと図に、次の3ファイルを追加する。

| ファイル | 内容 |
|---|---|
| cluster-members.csv | molecule, cluster_id, cluster_size, bond_degree |
| cluster-summary.csv | cluster_id, molecule_count, bond_count |
| network-summary.csv | molecules, bonds, clusters, largest, second_largest, isolated_molecules, largest_fraction |

`largest_fraction=largest/molecules`。N=0ならCSVは空欄（未定義）。
同じ指標とcluster_method、cluster_id_rule、clustering_completedをmetadataからrun.jsonへ記録する。
現行デモは2分子固定、4腕5歩、A/B固定、開境界、距離2、多重結合許可。cfgの入力はdim,L,seedだけ。
L8seed3は結合2本・クラスタ1個・最大2・第2位0・孤立0・最大比1。
L8seed12345は結合0本・クラスタ2個・最大1・第2位1・孤立2・最大比0.5。
小さな動作確認例であり、最大比1をゲル化の証拠とはしない。ゲル化点の推定には多分子・濃度・サイズの統計が必要。

## 計算量・検証・再利用

N分子・B結合でUF構築と所属取得はO((N+B)α(N))、CSV集計はO(N+B)、追加メモリO(N)。
α(N)は逆Ackermann関数。経路圧縮とサイズによる統合に基づく償却計算量で、実測倍率は未評価。
元のUF構造体は集計後に解放し、ClusterStoreの所属・サイズ・結合数・次数を残す。

鎖による間接接続、多重結合、閉路、孤立分子、同率最大2クラスタ、0分子、
結合順による出力ID不変、不正末端ID・逆引きの拒否をテストする。
実際に生成した6分子・50seed、両方の多重結合条件で、独立したBFSと所属・サイズ・次数を比較。
元のBondStoreを変更しないことも確認する。

そのまま再利用したもの：新リポジトリの`uf_init`、`uf_union`、`uf_find`、`uf_largest_two`、`uf_free`。
これらは準備段階で整備したUF基盤で、旧リポジトリの関数を今回直接呼んでいるわけではない。
生成・候補・結合・図・実行保存は既存の③-Bを直接利用。
今回新しく作ったもの：BondStoreから分子IDへの接続、入力検査、ClusterStore、安定した出力ID、3種類のCSV。
UFは結合相手の腕ID・座標・直接の隣接関係を保持しない。それらはBondStoreとMoleculeStoreを参照する。
