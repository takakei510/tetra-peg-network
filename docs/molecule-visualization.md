# 1分子の軌跡を画像にする

`scripts/visualization/plot_molecule.py`は、Cの生成処理が出力したCSVを読み、
4本の腕をMatplotlibで3D描画します。図を描くために新しく分子を生成する
処理はありません。Cで作った座標をそのまま使います。

## WSLでの初回準備

リポジトリ直下（例：`/home/takak/tetra-peg-network`）で実行します。

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r scripts/visualization/requirements.txt
```

`venv`が利用できない場合は、WSLのUbuntuで`sudo apt install python3-venv`を
実行してからやり直します。以後は`source .venv/bin/activate`で再開できます。

## CSV生成から画像保存まで

```sh
make test build/demo_molecule
mkdir -p data
./build/demo_molecule configs/default.cfg > data/molecule.csv 2> data/molecule-metadata.txt
python scripts/visualization/plot_molecule.py data/molecule.csv --output data/molecule.png
```

`data/molecule.png`が報告スライドに貼る画像です。VS Codeでこのファイルを
開けば確認できます。WSLでは`explorer.exe data`でWindows側からも開けます。
画像は300 dpiで保存します。メタデータにはseedや生成条件が記録されるため、
CSVと併せて残してください。

画面上で回転させたい場合は、GUIが利用できる環境で次を実行します。
画像保存だけならGUIは不要です。

```sh
python scripts/visualization/plot_molecule.py data/molecule.csv --output data/molecule.png --show
```

## 図の読み方と処理

- 黒丸：分子中心。色付き線と丸：腕の各ステップ。四角：末端。
- CSVの中心行は`arm=-1, step=0`。腕は`arm=0..3, step=1..l`。
- 各腕をstep順に並べ、中心を先頭に補ってから座標を線でつなぎます。
- 標準では中心を原点に平行移動します。`--absolute`なら元の座標です。
- 3軸を同じ尺度で表示します。図は視点による投影なので、画面上で線が
  重なって見えても同じ格子点を占有するとは限りません。
- 4腕の長さ・step・最近接移動・格子点重複を確認し、不整合なら停止します。
- `--molecule 0`で対象ID、`--dpi 300`で解像度を指定できます。

現在は開境界の4腕モデル専用です。周期境界のwrapped座標を展開する処理や、
複数分子を一度に表示する処理は含みません。タイトルと凡例は環境による
日本語フォントの違いを避けるため英語です。

## 以前のコードとの関係と報告文

このPython可視化コードは今回新規に作成しました。旧SAWコードの直接利用は
ありません。旧SAWの候補選択・行き止まり判定を参考にしたC側の生成処理が
CSVを出力し、このコードが読み込みます。

報告例：「実装したCコードで4腕分子を生成し、軌跡をCSVとして出力した。
その座標をPythonで可視化し、中心・各腕・末端を区別して表示した。」

画像は生成例の可視化であり、ゲル化やモデルの物理的妥当性を示す結果では
ありません。生成条件はメタデータを参照します。
