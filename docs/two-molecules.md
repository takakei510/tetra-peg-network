# 2分子配置と可視化（結合前）

## 今回の範囲

1個目の分子を残して2個目を生成し、分子間で格子点を重複使用しないことを
確認する段階です。末端間の結合・距離判定・ネットワーク解析はまだ行いません。
開境界、空近傍からの一様選択、腕を1本ずつ完成させる条件は暫定のままです。

`demo_two_molecules`は確認用にA型1個、B型1個を固定します。A/B比率を確率rで
抽選する実験とは区別してください。中心候補は全格子点から一様ランダムに
選びます。生成失敗時はその分子だけを復元し、既存の完成分子は残します。
腕長5、目標2個、全体の試行上限10000回はこのデモの固定条件です。

## 実行

```sh
make test build/demo_two_molecules
mkdir -p data
./build/demo_two_molecules configs/two-molecules.cfg > data/two-molecules.csv 2> data/two-molecules-metadata.txt
```

終了コード0が配置完了、1が未完了またはエラーです。未完了でも完成した分子の
CSVを出力する場合があるため、metadataの`placed=2 completed=1 occupied=42`と
終了コードを確認してください。正常時CSVは見出しを含め43行です。
metadataには中心占有による棄却、成長中の行き止まり、実績占有率も記録します。

Python環境は1分子可視化と同じです（初回準備はmolecule-visualization.md参照）。

```sh
source .venv/bin/activate
python scripts/visualization/plot_molecules.py data/two-molecules.csv --output data/two-molecules.png
```

青：A型、赤：B型、黒丸：中心、四角：末端。絶対座標を共通に使うため、
2分子の実際の位置関係を保持します。各分子を別々に原点へ移しません。
`--show`でGUI環境の回転表示も可能です。

## 再利用したコード

|コード|今回の利用|
|---|---|
|molecule.cのmolecule_try_generate|1分子試行・共通owner・失敗復元をそのまま呼び出す|
|molecule.cのmolecule_write_csv|既存の全登録分子出力をそのまま使う|
|config/rng/lattice|設定読込・中心抽選・占有格子を直接利用|
|plot_molecule.pyのread_molecule|CSV読込・腕の順序・最近接・内部重複検査を直接利用|

旧percolation-modelのSAWの考え方はmolecule.cに引き継がれていますが、今回
旧コードを新たにコピーしていません。追加のデモと複数分子描画は新規です。

## 確認

Cテストでは1個目の腕上の中心候補を棄却し、2分子完成時に42点、全点owner一致、
最近接移動、1個目の中心・型・軌跡・方向の保持を確認します。従来の部分成長後
trappedテストも継続して、既存分子を残した格子の完全復元を確認します。
Pythonでも各分子内と分子間の格子点重複を検査してから画像を保存します。
これは実装の整合性確認であり、物理モデルの条件確定やゲル化結果ではありません。

今回の実行例（L=12、seed=12345）：3回の試行で2個完成、中心占有棄却0回、
trapped棄却1回、42点占有、phi=42/1728。A中心(5,6,8)、B中心(11,8,0)。
試行数はこの1例の値であり、成功率の推定には使いません。
CSVの42点と座標保持、行順を逆転しても同じ読込結果になること、内部軌跡が
正しい2分子でも分子間で重複すれば可視化が拒否することを確認しました。
L=2の未完了例では終了コード1とcompleted=0、occupied=0を確認しています。
