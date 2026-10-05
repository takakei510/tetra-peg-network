# 1分子生成：ミニマム実装の記録

## 実装と暫定条件
- 3D立方格子、開境界、4腕、各5歩のデモ。
- 非占有近傍から一様に選択し、腕0から3を順番に完成させる。
- 行き止まりなら中心を含む分子全体の仮占有を復元。
- MoleculeStoreは完成した分子のみ公開する。
- 中心抽選はデモ側、固定中心での1試行は生成関数側に分離。
- 格子サイズとseedは既存configから読み込む。デモの腕長5、A確率0.5、試行上限10000は固定値で、stderrに記録。量産実験前にconfigへ移す。
- 先生の生成規則・境界条件との一致は未確認。ゲル化結果はまだない。

## 設計案からの簡素化
専用GenerationWorkspaceは今回導入しない。未登録分子の保存スロットを仮軌跡領域として利用し、成功時だけcountを増やす。棄却後に未登録スロットに値が残っても、登録済みデータには含めない。仮占有の復元用一覧のみ試行ごとに確保する。各配列を直接参照する際はcount未満だけを読む。

## 旧コードとの対応
|旧コード|活用した内容|今回の変更|分類|
|---|---|---|---|
|percolation-model/src/random_walk.c: run_one_walk() SAW分岐|近傍列挙、訪問済み点の除外、候補ゼロでtrapped|訪問配列を全分子のownerへ、4腕共有、軌跡保存、失敗復元|アルゴリズムを参考に再実装|
|同関数の候補選択|空候補から等確率で選択|rand()%nを既存rng_boundedへ変更|考え方を継承|
|同ファイルの訪問点解除処理|触った点の一覧で解除|中心と4腕の仮占有だけを解除|考え方を参考に再実装|
|新リポジトリのlattice/rng/config|座標変換、占有格子、乱数、設定読込|既存APIを呼ぶ|基盤を直接利用|

旧run_one_walk()をそのまま移植していない。MSD等の集計は生成処理から分離。方向列挙順も旧実装と異なるため、旧seedの軌跡一致は保証しない。

## 確認結果
- foundation testsとmolecule testsが成功。
- 1分子21点、最近接移動、owner一致を検証。
- 占有数と全軌跡点数の一致により中心・腕間の重複を検出。
- 同じseedと条件で軌跡・方向・乱数状態が一致。
- 1歩成長後に必ずtrappedになる条件で格子の完全復元、登録分子不変を確認。
- AddressSanitizer/UndefinedBehaviorSanitizerで検証（リーク検査は環境の制限で実行できず、detect_leaks=0で実施）。

## 実行
```bash
make test build/demo_molecule
./build/demo_molecule configs/default.cfg > data/molecule.csv 2> data/molecule-metadata.txt
```
CSV: molecule,type,arm,step,site,x,y,z,direction。中心はarm=-1,step=0,direction=-1。腕はstep=1..l。
方向0..5は+x,-x,+y,-y,+z,-z。格子ID変換は既存Latticeに従う（xが最も大きい桁）。

## 報告スライド案（2枚）
1. 分子生成：4腕SAWの軌跡を示し、非占有近傍選択と21点占有を説明。生成規則は暫定と明記。
2. 妥当性確認と旧SAWの活用：最近接・重複なし・失敗復元・seed再現性。旧SAWの候補選択を参考に、分子全体の占有判定へ拡張したことを説明。
