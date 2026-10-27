# 大規模H₂検証の保存容量対策（2026-10-27）

## 保存方針

run_h2_validation.pyの標準は--storage-policy compact。
従来のbenchmark.pyを単独で使う場合は互換性のためfullが標準。
Cの生成・探索・結合・乱数・クラスタ計算は変更しない。出力先だけを切り替える。

|データ|compactの動作|用途|
|---|---|---|
|config.cfg、実行ログ・完了記録・ソースZIP|保持|再現・失敗原因調査|
|bonds.csv.gz|保持|H₂の再計数、4・6分子サイクル解析|
|cluster-members.csv.gz、cluster-summary.csv.gz|保持|所属・次数・クラスタサイズの検証|
|network-summary.csv.gz、placement-summary.csv.gz|保持|全体統計・目標達成の検証|
|軌跡trajectories.csv、結合候補candidates.csv|/dev/nullへ出力し、ディスクに保存しない|形状・全候補の再調査が必要な実験ではfullを明示|
|試料別clusters.csv.gz|直接gzip出力・復号読込み検証後に保持|クラスタ別H₂・サイズの正本|
|全試料連結cluster_samples.csv|作成しない|重複保存を避ける。正本は各試料のgzip|
|送付ZIP|試料統計・サイズ区間統計・収束・図・manifest・storage_reportのみ|チャット等で確認するための小さい結果セット|

compactでは分子軌跡を保存しないので、元の空間配置をCSVから完全復元できない。
同一ソース・条件・seedから再生成するための情報は保存する。結合グラフの再解析は保持CSVだけで行える。
古い実験フォルダは自動削除・変換しない。新しいコード・保存方針を混ぜず、新保存先を使う。
--storage-policy fullなら軌跡・候補も保存するが、試料別クラスタ解析CSVはどちらもgzip、巨大な全試料連結CSVはどちらも作らない。

## 容量監視の範囲と限界

- 親プロセスは、実行中workerのC/圧縮領域・クラスタgzip領域・メタデータと、集計用256MiBの余裕を予約してから新試料を投入する。
- 同じstudy内の完了試料・失敗フォルダ・過去の結果フォルダも使用量に含める。
- Pythonが書くクラスタgzipに試料別上限を付ける。
- 最終CSV・JSON・PNG・ZIPの書込みは、空き容量と--max-storage-gibの残量から割り当てた共通予算を使う。書込み前に予算を確認する。
- ZIPは巨大な連結CSVを含まず、CRC確認後に完了扱いとする。
- 集計途中で制限に達したフォルダにはEXPORT_INCOMPLETE.jsonが残る。そのZIPを完成結果として扱わない。完了済み試料は保持され、予算を増やして再実行できる。
- 再実行のたびにresults_*が増える。古い結果も予算に入る。不要な結果を削除する場合は、ユーザー自身が対象を確認する。
- Cプロセスは既存の定期監視で停止するため、これは厳密なOSディスククォータではない。監視間隔中の増分や他アプリの書込みもあり、ディスク満杯寸前の予算設定はしない。
- 既定予算20GiB。--max-storage-gibは「保存を保証する容量」ではなく停止のためのstudy予算。より大きい指定で再開する場合も実ドライブの空きを確認する。

スケジューラは起動時にstudyの使用量を走査し、その後は完了試料の容量を加算する。全studyを試料投入ごとに再走査する従来の高コスト処理も減らした。

## デスクトップでの実行手順

Windows上の実空き容量はC約201.7GiB、D約154.6GiB、G約79.1GiBというユーザーの表示を確認済み。
WSLに見える約953GiBだけから使える容量を判断しない。i5-12400F、WSL使用可能メモリ約14GiBでは最大6workerから測定する。

```bash
cd ~/tetra-peg-network
git pull --ff-only
python3 tests/test_h2_storage.py
python3 tests/test_h2_relation.py
```

両方OKなら、全61占有率の設定を維持したまま、3サイズ×代表3占有率×2seed=18試料で容量を測る。

```bash
python3 scripts/run_h2_validation.py capacity \
  --config configs/experiments/h2-validation-fine.json \
  --output data/validation/h2_relation_fine_compact \
  --samples 2 --workers 6 --max-storage-gib 30
```

capacityは各Lに対し、最低Φ・0.44に最も近いΦ・最高Φを選ぶ（標準設定で0.10、0.44、0.70）。
旧pilotでは最大L・Φ0.44だけだったため、全体容量の見積もりにはcapacityを使う。
完了後はUPLOAD THIS FILEのZIPを送付する。storage_report.jsonに各Lの完了試料フォルダの平均・最大bytesを含む。
テスト自体も小さな一時データを作るが、TemporaryDirectoryでテスト終了時に片付ける。ユーザーの既存実験フォルダは触らない。

容量・結果を確認した後の本計算：

```bash
python3 scripts/run_h2_validation.py run \
  --config configs/experiments/h2-validation-fine.json \
  --output data/validation/h2_relation_fine_compact \
  --samples 100 --workers 6 --max-storage-gib 30
```

18試料は再利用。183条件×100=18,300試料。30GiBで止まる場合は、結果を見て拡張や実験範囲を検討する。
統計のために失敗seedを別seedで置き換えない。

## こちらの環境での測定（PC上の容量保証ではない）

L=32/64/128、Φ=0.10/0.44/0.70、各2seedをcompactで実行。18試料・集計ZIPまで完成。

|L|Φ0.10の平均MiB/試料|Φ0.44|Φ0.70|
|---:|---:|---:|---:|
|32|0.0114|0.0295|0.0448|
|64|0.0292|0.1826|0.3320|
|128|0.1845|1.5389|2.8988|

試料フォルダは原CSV gzip・クラスタgzip・結果JSON・ログを含む。全studyは約10.63MiB（ソース・最終集計を含む）。送付ZIPは約15KiB。
各Lの観測最大値×各L6100試料を足すと約19.54GiB。これは代表条件に基づく概算で、全密度の最大値を保証しない。
100試料取得前から500/1000へ広げず、実測storage_reportと誤差の収束を確認する。

## 検証

- 同一seed・条件でcompactとfullの結合CSV・所属・クラスタCSVの復号バイト列、および科学的な集計値が一致。
- gzipの読込みとZIPのCRC確認。
- 書込み予算を跨ぐデータは書込み前に拒否。最終集計中の容量不足でも既存試料を保全。
- 小さい全体予算では試料投入前に停止し、予算増加後に同じstudyで実行できる。
- 既存の並列実行・再利用・試料数追加・改変検出テスト、および既存H₂解析と詳細計算のテストも実行。
