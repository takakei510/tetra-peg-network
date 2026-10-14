# ②：設定した分子数で生成・結合・UF集計

## 実行

```sh
python scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg
```

新形式のcfgはnetworkモード専用。既存の2分子デモは固定条件を維持する。
`target_molecules`か`target_phi`の片方を指定し、Cのnetwork_config_loadで目標Nに統一する。
Python側で換算を繰り返さず、Cのplacement-summary.csvから解決済み設定と実績を読む。

## 生成と結合の順序

1. 配列容量を目標N、各腕の軌跡長をarm_lengthとして初期化する。
2. 全格子から中心候補を一様に選ぶ。占有中心は棄却。空中心ならA型確率rで型を抽選する。
3. molecule_try_generateで4腕を順に成長させる。既存分子・同分子の占有点を回避する。
4. 成功時だけ完成分子数が増える。行き止まりでは現在の仮配置だけ取り消して再試行する。
5. 目標Nが完成してから24近傍の末端探索、A–B結合、UF集計を実行する。

中心の提案は全格子からの一様抽選。空点のみの配列を毎回作って抽選する方式ではない。
境界は開境界、成長は空近傍から一様選択、4腕の順序は0..3。
結合距離はマンハッタン1〜2。末端順をシャッフルして、その時点の未結合候補から相手を一様抽選する。
末端は高々1結合、同じ分子ペアの別腕による複数結合は許可する。
これらは現在採用した条件であり、資料の未確定条件を解決したと断定しない。
生成に使った同じRng状態を結合へ渡すため、seedだけでなく生成規則や乱数消費順も再現に影響する。

## 試行上限と未達

networkモードの既定上限は全分子合計10000回。中心占有も行き止まりも1回の試行に数える。
`--max-placement-attempts 1`などで変更でき、cfgとは別にrun.jsonとplacement-summary.csvへ記録する。
追加cfgキーではない。極端な高濃度では、必要点数が格子内に収まっても目標達成を保証しない。

| 状態 | 終了コード | 保存・後続処理 |
|---|---|---|
| completed | 0 | 全軌跡、配置集計、候補、結合、クラスタ集計、図 |
| incomplete | 3 | 完成した分子の軌跡、配置集計、未達状態。結合・UFは実行しない |
| failed | 1 | 設定・メモリ・入出力・描画等の実行エラー。得られたファイルと理由を保持 |

未達でも完成分子は消さない。非空なら未達と明記した軌跡図を作り、0分子なら図を作らない。
失敗した分子の途中の形の履歴は保存しない。center_occupiedとtrappedの回数は保存する。
メモリ不足は物理的な棄却として数えず実行エラーにする。
未達試料を目標占有率での完成ネットワークとして扱わないため、候補・結合・UFをスキップする。

## 保存する値とファイル

`data/runs/`以下に実行ごとのフォルダを作る。過去結果を上書きしない。

| ファイル | 内容 |
|---|---|
| config.cfg | 入力条件の原文 |
| placement-summary.csv | 入力モード、入力占有率、目標N、換算後占有率、完成数、実占有点数・率、A/B数、試行上限・回数、棄却回数、完了・停止理由 |
| trajectories.csv | 完成分子の中心と各腕の軌跡 |
| candidates.csv / bonds.csv | 候補ペアと成立した結合を区別 |
| cluster-members.csv / cluster-summary.csv / network-summary.csv | 分子の所属、クラスタごとの分子数・結合数、ネットワーク集計 |
| molecules.png | 軌跡と成立結合。多分子時はA/B型ごとの凡例 |
| run.json / metadata.log | 実行状態・コマンド・目標と実績・固定モデル条件 |
| source.zip / build.log | 再現用コードとビルド記録 |

N指定ではrequested_phi欄は空。実現A比率はtype_A/placed_moleculesから求める（0分子では未定義）。
結合率は2*bonds/(4*placed_molecules)、最大クラスタ比はlargest/placed_molecules。
実績率と入力率を区別し、0分子の比率を0と決めつけない。

## 19分子の確認結果

L16、腕長5、r0.5、phi0.1、seed3、上限10000の1回実行：
19分子（A11/B8）、399占有点、実占有率0.097412109375。
23試行＝19成功＋占有中心2棄却＋行き止まり2棄却。
76末端、候補7組、結合3本、16クラスタ、最大2分子、孤立13分子。
全19分子が一つにつながるとは限らない。1回の確認でありゲル化点の推定ではない。

## 既存コードの利用と検証

直接呼ぶ既存関数：molecule_try_generate/molecule_store_init（生成と軌跡）、
endpoint_index_build/endpoint_find_candidates（索引と24近傍探索）、
bond_select（末端のランダム結合）、cluster_build（既存UFで分子を集計）。
新規処理はplacement_generateの反復・試行集計とdemo_networkの接続・保存。
Pythonはrun_demo.pyの実行フォルダ保存を拡張し、plot_molecules.pyのCSV検査・3D描画を利用する。

make testでは24seedの目標達成・占有整合・再現性、1回上限の部分保存、L2で0完成の完全復元を検査。
Python統合テストでは実行入口を通し、CSV再現性、全末端対との候補比較、末端の1回使用、
独立なBFSとのクラスタ一致、N直接指定・A型のみで結合0、部分/0完成の結果保存を検査する。
