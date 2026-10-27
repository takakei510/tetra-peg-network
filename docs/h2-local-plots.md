# PC上で集計CSVから図を作る

再シミュレーション不要。図のために数GiBの原データをコピーする必要もない。
作図入口はscripts/analysis/plot_h2_validation_results.py。元のrun_h2_validation.pyを実行し直さない。
この作図スクリプトは入力CSVやstudy.jsonを変更しない。

## Windows + WSL（今回の環境）

```bash
cd ~/tetra-peg-network
git pull --ff-only
sudo apt update
sudo apt install python3-matplotlib
python3 scripts/analysis/plot_h2_validation_results.py
```

sudoのパスワードはWSLのLinuxユーザーのもの。入力中は文字が表示されない。
デフォルト入力はdata/validation/h2_relation_fine_compact。
完了済みresults_*のうち試料数が最も多いものを選ぶ。同数なら新しいものを使う。
今回の18,300試料が存在すれば18試料capacityを誤って選ばない。
不完全なexport（EXPORT_INCOMPLETE.jsonあり）やpartialは選ばない。

data/figures/h2_<日時>_<ID>/へ、6種類×PNG/SVGと説明・作図条件を出力。
出力先は毎回新規で、既存の図を上書きしない。
最後のUPLOAD THIS FILEに表示するupload_h2_plots.zipをチャットへ添付する。

```bash
explorer.exe "$(wslpath -w ~/tetra-peg-network/data/figures)"
```

このフォルダ内の今回のh2_<日時>_<ID>を開く。
GoogleドライブへZIPを保存すれば、ノートPCでも同じ図を確認できる。
PNGは閲覧とPowerPoint用、SVGは拡大や編集用。
今回の図一式は数MiB程度。元の全実験データはデスクトップに保持する。

## 別の結果・条件を指定

```bash
python3 scripts/analysis/plot_h2_validation_results.py \
  --results data/validation/h2_relation_fine_compact/results_XXXXX
```

指定フォルダにはsummary.csv、size_relation.csv、convergence.csv、manifest.jsonが必要。
他PCでは「今回の結果ZIP」を展開したフォルダを--resultsへ渡せる。Matplotlibが必要。
03等の全体図は全条件、04/05のサイズ関係図は標準でΦ0.44・L128。
別のΦ・Lを見たい場合：

```bash
python3 scripts/analysis/plot_h2_validation_results.py --phi 0.45 --size 128
```

作図だけのため元の実験を再実行しない。今回のGit更新は実験用バイナリを変更しないが、
既存実験runnerは広い範囲のソースハッシュを比較するため、後日seedを追加する際は履歴と互換性の確認が必要。
ハッシュ不一致を無視してstudy.jsonを書き換えない。

## 図の意味

- 01：H₂密度・全末端割合・結合済み末端割合。横軸は実占有率。
- 02：最大1個除外のSfinと、そのピーク付近の拡大。
- 03：生のK、K/L³、格子全体に対する最大クラスタ占有率。
- 04：Φ0.44での3サイズ比較。最大クラスタを含む/1個除外。白抜きは有効seedが全試料の半数未満。
- 05：L128・Φ0.44のサイズとH₂末端割合、寄与seed数。
- 06：先頭n試料での平均・SE。各LでSfin平均が最大の測定点を使用。入れ子で独立な実験ではない。

すべて点とSE。補間・fit・点をつなぐ線なし。背景格子線は値を読むためのもの。
サイズ関係の縦軸は100h2_i/s_i (%)、まず試料内でクラスタ平均、次にseed平均。
有効seed数はサイズ区間ごとに異なる。SE未定義（n=1）なら誤差棒を描かない。
観測ピークは真の転移点ではなく、クラスタサイズとH₂割合の関係は因果の証明ではない。

共有後、各図の「目的・定義・計算・観測・限界」を確認してから来週用スライドへまとめる。

## 確認済み

18,300試料のCSVでPNG/SVG計12枚と共有ZIPの作成・CRC、入力4ファイルが変更されないことを確認。
完了試料数による入力選択とpartial除外はpython3 tests/test_h2_plot_selection.pyで検証。
