# Φ=0.10〜0.70、全域0.01刻みの追加設定

目的・定義・解析方法は[h2-validation.md](h2-validation.md)と同じ。
旧設定の「0.40〜0.50だけ0.01刻み」を、全域0.01刻みへ拡張する。

- L=32/64/128、61占有率、計183条件。
- 各条件100試料なら18,300試料、500なら91,500試料。
- 並列単位は(L,Φ,seed)。--workers 4は最大4プロセスで、CPUと検出空きメモリで制限する。
- 新しい図は測定点とSEのみ。線・補間・fitを追加しない。点が増えても統計誤差や偏りが消えるわけではない。
- 出力容量予算は従来どおり既定20GiB。到達時は追加投入を停止して部分結果を保存する。

```bash
cd ~/tetra-peg-network
git pull --ff-only
python3 scripts/run_h2_validation.py run \
  --config configs/experiments/h2-validation-fine.json \
  --output data/validation/h2_relation_fine \
  --samples 100 --workers 4
```

旧設定とは条件一覧が異なるため保存先を分ける。旧結果・途中結果は保持される。
新設定は旧3試料を自動取込しない。同じnamespaceと共通条件は同じseedとなるので、旧結果と新結果を足して独立試料数を増やしたことにしない。
旧計算が実行中なら同時に走らせず、Ctrl-Cで終了を待ってから実行する。
容量で停止した場合、空き容量を確認し、必要なら同じコマンドに--max-storage-gib 50等を付けて再開する。
完了済み試料は再利用される。同じ設定・保存先で--samples 500へ拡張できる。
最後のUPLOAD THIS FILEに表示されるZIPを確認する。manifestのcompleted_samplesとrequested_samplesが一致して初めて計画完了。

計画確認（シミュレーションは実行しない）：

```bash
python3 scripts/run_h2_validation.py plan \
  --config configs/experiments/h2-validation-fine.json \
  --samples 100 --workers 4
```

検証：61占有率の端点・0.01刻み、18,300件のseed重複がないこと、旧共通条件とseedが一致すること、既存runnerでplanが通ることを確認。
生成・並列処理・解析のコードは変更していない。本条件18,300試料の生成はPCで実施する。
