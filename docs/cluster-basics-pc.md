# 基本クラスタ量を18,300試料で再作図する

H2の列を使わず、保存済みsamples.csvからネットワークの基本量を再集計します。シミュレーションの再実行や原データの変更は行いません。

## PCで実行

リポジトリを開いて実行します。

```bash
cd ~/tetra-peg-network
git pull --ff-only
python3 scripts/analysis/plot_cluster_basics.py
```

入力の既定は `data/validation/h2_relation_fine_compact` です。そこにある `results_*` のうち、manifestが完了状態で、完了試料数と要求試料数が一致するものを候補にし、最も多い試料数、その次に更新時刻で選びます。選んだ入力パスは実行後に表示します。

入力を明示する場合：

```bash
python3 scripts/analysis/plot_cluster_basics.py --results /path/to/results_XXXXXXXX
```

フォルダから取り出したsamples.csvだけでも実行できます。

```bash
python3 scripts/analysis/plot_cluster_basics.py --results /path/to/samples.csv
```

WindowsのPowerShellならリポジトリ内で `py scripts/analysis/plot_cluster_basics.py` を使えます。別環境ではMatplotlibと日本語フォントが必要です。普段のWSLでMatplotlibを使えていれば、同じ環境で実行してください。

依存を追加する場合は使っているPython環境で `python3 -m pip install matplotlib`。WSLの日本語フォントは `sudo apt install fonts-noto-cjk`。明示指定は `--font /path/to/font.ttf`。

## 出力

`data/figures/cluster_basics_<日時>/` にPNGとSVGを各7種類出力します。日時ごとに別フォルダを作り、既存の結果は上書きしません。

|ファイル|何を確かめるか|
|---|---|
|01_largest_cluster|最大集団の絶対分子数N1と、全分子に占める割合N1/N|
|02_mean_cluster_sizes|平均の取り方と最大クラスタの除外で、何が変わるか|
|03_bonding_and_largest_volume|濃度を増やしたときの結合率と最大クラスタ占有率|
|04_largest_volume_vs_bonding|末端がどれだけ結合した状態で、大きな集団が形成されるか|
|05_cluster_counts|全体が何個の集団に分かれるかと、体積あたりの個数|
|06_sfin_peak_region|最大以外の大きな集団が発達する領域と有限サイズ依存|
|07_sample_count_stability|5〜100試料を使ったときの平均とSEの推移|

`figure_guide.md` に、各図の目的・定義・物理的意味・説明例・限界を出力します。報告前にこの文章と図を対応させて確認してください。

`basic_summary.csv` は1条件1行で、各量の有効数・平均・SD・SEを保存します。`peak_summary.csv`、`sample_count_stability.csv`、入力とコードのSHA-256を含む `run_metadata.json` も保存します。

共有するときは最後に表示される `cluster_basics_share.zip` を送ってください。図・表・説明・実行記録をまとめて確認できます。元のsamples.csvはZIPには入れません。

## 定義と比較上の注意

- Nは孤立分子を含む総分子数。N1は最大クラスタに属する分子数。
- Pmax=N1/N、S1=21N1/L³。分母が違うので同じ割合ではありません。
- A=N/K、Sall=Σs²/N、Sfin=Σ′s²/Σ′s。Sfinは各試料の最大クラスタを1個除外します。
- Sfinが大きいことはクラスタ数Kが多いことではありません。
- p=2E/(4N)は、生成・結合した結果の測定値です。入力で指定した反応確率ではありません。
- 横軸は実占有率、条件の表やピーク注記は目標占有率です。
- 比率は各試料で計算してから平均し、試料をプールした比へ置き換えません。
- SEは平均の不確かさ、SDは試料間のばらつき。欠測Sfinは0とせず、有効数を記録します。
- 誤差棒はSEです。線・フィット・転移点の推定・貫通判定は追加していません。
- 07は今回の100試料の先頭部分集合を使う例示であり、前回の5試料そのものとの直接比較ではありません。各点は入れ子です。
- このCSVには第2クラスタサイズや貫通の情報がないため、推測して作図しません。
- 旧5試料と新100試料は混ぜず、今回の完了データだけを使います。モデル変更の検証や元結合からの再構築とは別の検証です。

## 設定

`configs/analysis/cluster-basics.json` で対象条件・試料数、1分子の占有点数、腕数、拡大窓、出力形式を設定します。今回の既定は3サイズ×61占有率×100試料です。条件不足・重複試料・保存値と定義の不一致があれば停止します。

別の実験に使う場合は設定をコピーし、条件集合と期待試料数を明示して `--config コピーした設定.json` を指定してください。1分子の占有点数や腕数もモデルに合わせる必要があります。

## 確認

```bash
python3 -m unittest discover -s tests -p test_cluster_basics.py -v
```

平均の取り方、未定義値、重複試料、不正なモーメントを確認するテストです。受領済み18,300試料でも全図を生成し、既存summaryの基本量の平均・SD・SEと一致を確認しています。
