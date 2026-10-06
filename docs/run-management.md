# 実行条件と結果をひとまとまりで保存する

## 実行方法

リポジトリのルートで、Matplotlibを入れたPython環境から実行する。
`make`とC11コンパイラが必要。実行対象は自動でビルドされる。

```sh
python scripts/run_demo.py --mode single --config configs/demos/single_L32_seed12345.cfg
python scripts/run_demo.py --mode two --config configs/demos/two_L8_seed12345.cfg
```

CSVだけ必要なら `--no-plot` を付ける。この場合Matplotlibは不要。
設定を変更して試す場合は、例えば次のように見本をコピーする。

```sh
mkdir -p configs/local
cp configs/demos/two_L8_seed12345.cfg configs/local/my-two.cfg
python scripts/run_demo.py --mode two --config configs/local/my-two.cfg
```

## 1回の結果フォルダ

`data/runs/<UTC日時>_<single・two・candidates>_<識別子>/` に毎回新規保存する。
同じ条件で再実行しても上書きしない。フォルダ名はUTCを示すZ付きで、日本時間とは9時間差がある。

| ファイル | 用途 |
|---|---|
| `config.cfg` | 実際にCへ渡した設定のコピー |
| `run.json` | 実行ID、開始・終了時刻、成功/失敗、コマンド、コード情報、実現分子数・型・占有数、ファイルのSHA-256 |
| `trajectories.csv` | 全分子の中心・腕の軌跡。1分子/2分子とも同じ列形式 |
| `metadata.log` | Cが出力したL・seed・固定条件・試行情報 |
| `molecules.png` | CSVから描いた図。描画成功時のみ存在 |
| `build.log` | コンパイルの出力 |
| `plot.log` | 描画の出力。描画を指定した場合のみ存在 |
| `source.zip` | 実行時のsrc・include・Pythonスクリプト・Makefile・requirements.txtの控え |

コードの控えには未コミット変更も含める。GitコミットIDに加え、実行バイナリのハッシュも残す。
コンパイラやMatplotlib等の完全な環境固定はまだ行わないため、異なる環境で画像の見た目まで一致する保証ではない。
描画は既存のCSV検証を通す。1分子は中心相対、2分子は共通絶対座標という既存の表示を使う。

## 成功と失敗の区別

- 正常終了：`run.json` の `status=completed`、終了コード0。
- `--no-plot` の正常終了：生成は成功、`plot_generated=false`。
- 生成・ビルド・描画の失敗：`status=failed`、終了コード1。ログと作成済みデータは残す。
- 描画だけ失敗した場合：`generation_completed=true`。CSVは残るが、一連の処理は失敗として扱う。
- JSONの `reported_conditions` はCのmetadataから読み取った文字列。数値集計には明示的に型変換する。

CSVは成長中の試行履歴ではなく、完成して登録された分子の軌跡である。未完了実行は目標濃度の試料として扱わない。

## Gitとバックアップ

設定の見本・コード・説明書をGitで管理する。`data/`内の実行結果と`configs/local/`はGit対象外。
実験結果をGitへ大量に追加しない。大事な実行フォルダは、フォルダごと別の保存先へバックアップする。
Git対象外は自動バックアップされるという意味ではない。このスクリプトはバックアップ先への転送は行わない。
スライド用の図は該当実行のものを使い、どの実行IDを使ったか記録する。
過去の `data/molecule.csv` 等は自動移動・削除しない。条件情報が不足している古い結果へ推測でmetadataを追加しない。

## 既存コードの活用

生成は既存の `demo_molecule` / `demo_two_molecules` を呼び、図は既存の `plot_molecule.py` / `plot_molecules.py` を直接利用する。
今回はPython標準ライブラリで実行・記録をまとめる層を追加した。旧SAWの生成アルゴリズムの変更はない。

## 今後の拡張

生成条件の設定化、結合CSV、クラスタ集計も同じ実行フォルダへ追加する。
大規模実験ではコードの控えを各試行で複製せず、共通のコード版と試行IDで整理する方式へ拡張する。
並行実験は専用の実行基盤で扱う。現行ラッパーの共有buildディレクトリへの同時ビルドは想定していない。

## ③-Aの候補探索

`--mode candidates --config configs/demos/candidates_L8_seed3.cfg` で
既存の保存内容に `candidates.csv` を追加する。画像には候補を緑点線で表示する。
候補0件も正常。距離規則は固定のmanhattan・半径2で、metadataからrun.jsonへ記録する。
候補は未確定のペアであり、結合数ではない。詳しくはendpoint-search.md参照。
