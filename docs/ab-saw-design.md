# AB-SAW設計書 v0.3

更新日：2026-10-15。状態：設定からの多分子生成、末端探索・結合、UF集計、CSV保存・Python描画まで実装済み。追加として分子グラフと代表最大成分の抽出・描画を実装。
本書は初期設計と現実装の仕様をまとめる。資料の未確定な物理条件と、現在採用している実装条件は区別する。
既存7 CSVの列名・意味・計算式・保存条件は9〜21節、派生グラフ出力は23〜24節を参照。

## 1. モデル条件と未確定事項

ERATO公開シンポ資料「Tetra-PEG型SAWネットワーク（AB-SAWモデル）」の占有・結合パートを参照。
資料に明記：3D立方格子、ランダムな空中心、確率rでA/1-rでB、各腕lステップ、分岐数f、失敗時は中心選択段階から再試行。
例：f=4、l=5、m=21、結合距離lr=2（近傍24点）。各末端は高々1結合。

未確定：境界条件、成長時の方向選択（空近傍から選ぶか、6方向から選んで占有なら失敗か）、腕の順序、距離の定義、結合候補のランダム化手順。
「近傍24点」はマンハッタン距離1～2の格子点数と一致する。一方、ユークリッド距離2以下では32点。したがってユークリッド距離を無断で採用しない。距離規則は要確認。
失敗時の完全ロールバックは安全な設計として採用するが、資料の再試行記述だけで生成分布の全仕様が確定したとはみなさない。
現在の採用条件：開境界、空近傍から一様に成長、腕0〜3を順次生成、マンハッタン距離1〜2、末端順シャッフル後に未結合相手を一様抽選。同分子ペアの別腕による複数結合を許可する。これらを資料の未確定事項が全て確定したという意味にはしない。

## 2. 責務と正本

| 要素 | 正本となる情報 | 方針 |
|---|---|---|
| Lattice | site→molecule ID、境界 | -1空点。座標・近傍変換を提供 |
| MoleculeStore | 分子ID、A/B、中心、腕軌跡 | IDは成功順、失敗で欠番を作らない |
| Arm | 中心を除くl点の順序 | 中心はMoleculeのみ保持 |
| Endpoint | molecule ID、arm ID | 座標とタイプは分子から取得 |
| CandidatePair | 2末端ID | 重複なし、結合情報とは別 |
| BondStore | 成立した2末端ID | endpoint→bond IDを逆引きとして維持 |
| MolecularGraph | 分子対の辺、解析結果 | BondStoreから構築する派生データ |

4腕を初期モデルの固定条件とする。一般分岐数への変更はモデル拡張として扱う。
分子軌跡は連続配列path[(molecule*4+arm)*l+step]で保持し、腕ごとのmallocは避ける。
endpoint ID=4*molecule ID+arm ID。Endpoint配列を重複保存せず、一覧を列挙できるAPIを提供する。
IDはsize_t、未設定はSIZE_MAX。既存ownerはintなので、生成前に分子数がINT_MAX以下か検証する。
配列サイズの乗算・加算はオーバーフロー検査後に確保する。

## 3. 1分子生成API（現在の実装）

```c
GenerationStatus molecule_try_generate(
    Lattice *lat, MoleculeStore *store, Rng *rng,
    size_t center, double a_probability);
```

1回の呼び出し=中心候補1回の試行。内部で無制限再試行しない。
状態：GEN_SUCCESS / GEN_CENTER_OCCUPIED / GEN_TRAPPED / GEN_INVALID_CONFIG / GEN_ALLOCATION_ERROR。
中心候補は呼出側で抽選する。多分子ループの試行統計はPlacementStatsで集計する。
A/Bは空中心確認後に抽選し、方向番号は+x,−x,+y,−y,+z,−zを0〜5で保存する。乱数消費順も再現に影響する。
試行前に必要な格納領域と1+4l点の作業領域を確保する。
中心と腕を仮占有し、作業領域に書き込み点を記録する。trappedならその点のみ全解放する。
成功時のみ分子数を増加し、完成軌跡を登録する。失敗時は格子・登録済み分子を変更しない。乱数状態と試行統計は進む。
メモリ不足は物理的棄却として数えず実行エラーにする。

## 4. 複数分子配置と設定

target_moleculesを配置目標の正本とし、実績占有率は完成数Nを使ってN*(1+4l)/L³で算出する。
目標分子数target_molecules、または入力占有率target_phiの片方をcfgで指定する。
target_phi指定ではfloor(target_phi*L³/(1+4l))へ換算し、requested_phi・resolved_phi・actual_phiを区別する。
現在の共通cfgキーはdim,L,seed,arm_length,a_probability。試行上限はnetworkモードのCLI引数--max-placement-attempts（既定10000）。
boundary、growth_rule、trial_idはcfgキーとして未実装。追加すると未知キーのエラーになる。
未確認のモデル規則に暗黙の既定値を与えず、実験用選択と資料再現を区別する。
max_placement_attempts到達による目標未達は未完了状態として出力し、目標濃度の試料として解析しない。
中心占有による棄却と、成長開始後のtrappedを別集計する。
Aの入力確率rはcfgに、実現数N_A/N_Bは配置CSVに記録する。実現A比率N_A/Nはこれらから算出する（専用CSV列は未出力、N=0では未定義）。分子数を固定したA/B数の割当とは異なる。

## 5. 軌跡と境界

初期1分子テストでは既存の開境界APIを利用できるが、これを先生のモデルの境界と断定しない。
周期境界を採用する場合はwrapped site IDに加え、各歩の方向を保存する。大きな分子の重心半径・MSDをwrapped座標差から計算しない。
同一腕・別腕・既存分子すべてとの点重複を禁止し、共通中心だけは4腕の起点として共有する。

## 6. 後続の結合・グラフ設計

異種末端の距離条件から候補を生成し、末端を1回だけ使うペアリングを独立モジュールにする。
現在の実装は末端ID全体の順をシャッフルし、その時点で未結合の候補相手から一様に1個選ぶ。候補ペア全体をシャッフルする案とは区別する。極大マッチングであり、最大マッチングや全マッチングの一様抽出とは異なる。
同一分子は単一A/B型のためA-B条件で自己結合は除外される。
結合APIは同一分子対の複数結合の可否を引数で指定できる。現在のnetwork実行入口は許可に固定し、cfg/CLIで切り替える機能は未実装。結合次数と異なる隣接分子数を区別する。
Union-Findは全登録分子（孤立分子を含む）に対するクラスタ解析用。直接結合の正本にはしない。
閉路階数E-V+Cは多重グラフの独立閉路数であり、すべての単純閉路の個数ではない。

## 7. 検証・出力・計算量

1分子生成の検証：占有数1+4l、各腕l点、中心から各歩が最近接、中心以外の重複なし、owner一致。
失敗検証：部分成長後のtrappedでも格子を完全復元、既存分子を破損しない、分子数不変。
再現性：同じ設定・seed・実装バージョンで同じ軌跡と状態。成長順や方向列挙順の変更も記録対象。
複数分子検証：occupied=N*(1+4l)、endpoint_count=4N、N_A+N_B=N。
結合検証：bonded_endpoints=2E、p=2E/(4N)、各末端高々1結合、異種条件・距離条件を満たす。N=0ではpを未定義として扱う。
最大クラスタの分子数N_1と全分子数に対する割合N_1/Nを区別する。資料のS_1=m*N_1/L³とは分母が異なるため、CSVではlargest_fractionと命名する（17節）。
実装済みCSVはtrajectories、placement-summary、candidates、bonds、cluster-members、cluster-summary、network-summaryの7種類（9節）。設定・seed・バージョン等はconfig.cfg、metadata.log、run.jsonに保存する。
配置試行の時間O(4l)、ロールバックO(4l)、軌跡メモリO(N*4l)、格子メモリO(L^3)。
固定距離の候補探索は末端の空間索引でO(endpoint_count*K+candidate_count)を目指す。全末端対O(N^2)は小規模検証専用。

## 8. 実装の区切り

第1段階（実装済み）：データ格納・1分子試行・失敗復元・軌跡CSV。
第2段階（実装済み）：複数配置・試行集計と実占有率・3D可視化。成功割合は保存した完成数/試行回数から算出できる。
第3段階（現在の採用規則で実装済み）：候補・結合・UFによる分子クラスタ集計。資料の未確定条件の確認と統計実験は今後の課題。
生成成功率は生成アルゴリズム依存の診断量として扱う。臨界指数やゲル化判定はこの段階の成果に含めない。

## 9. 実行フォルダとCSVファイル一覧（実装済み）

この節以降は2026-10-15時点の実装仕様。主な対象は以下のnetworkモードである。

```sh
python scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg
```

結果は`data/runs/<UTC日時>_network_<識別番号>/`に保存する。各実行は別フォルダにし、過去結果を上書きしない。
フォルダ名の日時末尾の`Z`はUTC。日本時間はUTCに9時間を加える。

| CSVファイル | 1行の単位 | 役割 | 出力する関数 |
|---|---|---|---|
| trajectories.csv | 中心または腕の1点 | 完成分子の軌跡の正本 | molecule_write_csv |
| placement-summary.csv | 実行全体の1行 | 生成目標・実績・試行・停止理由 | demo_network.c内の保存処理 |
| candidates.csv | 候補の末端ペア1組 | 距離・異種条件を満たす候補。結合確定ではない | candidate_write_csv |
| bonds.csv | 成立した結合1本 | 実際に結ばれた末端ペアの正本 | bond_write_csv |
| cluster-members.csv | 分子1個 | 分子IDから所属クラスタを引く | cluster_members_write_csv |
| cluster-summary.csv | クラスタ1個 | 各クラスタの分子数・結合数 | cluster_summary_write_csv |
| network-summary.csv | ネットワーク全体の1行 | 全体の分子数・結合数・クラスタ統計 | network_summary_write_csv |

すべてカンマ区切り、先頭行は列名。整数IDは0始まり。列順は以下のヘッダのとおり。
比率は0〜1の値（例：0.1は10%）で、百分率ではない。距離・座標は格子間隔を1とした単位。
CSVには基本的に単位・型の別行を入れず、この設計書に意味を定義する。

### 9.1 実行状態によるファイルの違い

| 状態 | 終了コード | trajectories.csv / placement-summary.csv | 候補・結合・クラスタの5 CSV | 図 |
|---|---|---|---|---|
| completed | 0 | 保存 | 保存 | 通常はmolecules.pngを保存 |
| incomplete（試行上限で目標未達） | 3 | 完成済み分子と未達実績を保存 | 解析を実行せず、作らない | 完成数が1以上なら未達と明記した軌跡図 |
| failed（設定・確保・入出力・描画等のエラー） | 1 | エラーの段階までに得られたファイルを保持。完全な出力は保証しない | 段階によって有無・完全性が異なる | 段階による |

`--no-plot`指定時は状態にかかわらず図を作らない。0完成の未達では軌跡CSVはヘッダだけとなり、図を作らない。
候補0組・結合0本でも、正常完了なら対応CSVをヘッダだけで保存し、孤立分子のクラスタ集計を行う。
「ファイルがない」と「正常に0件だった」は区別する。最初に`run.json`のstatusと成果物一覧を確認する。
失敗した途中の分子形状は保存しない。完成した分子だけを軌跡に残し、棄却回数を別に保存する。

### 9.2 従来モードとの違い

| run_demo.pyのmode | 主なCSV |
|---|---|
| single / two | trajectories.csv |
| candidates | trajectories.csv、candidates.csv |
| bonds | 上記＋bonds.csv |
| clusters | 上記＋cluster-members.csv、cluster-summary.csv、network-summary.csv |
| network | 上記＋placement-summary.csv。分子数・腕長・A型確率を新形式cfgから読む |

従来デモの未完了は終了コード1であり、networkモードの未達専用コード3とは異なる。
`make test`のPython統合テストも現在は5回の確認実行を`data/runs/`に残す。これらは`--no-plot`で実行する。
通常実験とテストの保存先の分離は未実装で、画像の有無だけで実験条件を判定せずconfig.cfgとrun.jsonを確認する。

## 10. 共通のパラメータ・ID・座標

以降、目標分子数を`N_target`、実際に完成した分子数を`N`、腕長を`l`、成立結合数を`E`、クラスタ数を`C`と表す。
式中のNと、C変数やCSV列名のtarget_moleculesを混同しない。

### 10.1 入力条件をどこに保存するか

| 入力キー・引数 | 意味 | 保存先 |
|---|---|---|
| dim | 次元。現在3のみ | config.cfg、metadata.log |
| L | 各方向の格子点数。全格子点数はL³ | config.cfg、metadata.log |
| seed | 生成・結合に使う乱数の初期seed | config.cfg、metadata.log |
| arm_length | 腕1本の歩数l。中心を除く | config.cfg、metadata.log |
| a_probability | 空中心確認後にA型を選ぶ入力確率r | config.cfg、metadata.log |
| target_molecules | 目標分子数N_targetを直接指定 | config.cfg、placement-summary.csv |
| target_phi | 入力する目標占有率。分子数へ換算して使用 | config.cfg、CSVではrequested_phi |
| --max-placement-attempts | 全分子合計の試行上限。既定10000。cfgのキーではない | run.json、metadata.log、placement-summary.csv |

target_moleculesとtarget_phiは片方のみ指定。共通5キー（dim,L,seed,arm_length,a_probability）は必須。
各CSVにLやseedを重複させていないため、CSVだけを切り離さず実行フォルダ単位で保存する。

### 10.2 IDの定義とファイル間の対応

| 名前 | 決め方・範囲 | 用途・注意 |
|---|---|---|
| molecule / molecule_a / molecule_b | 完成した順に0〜N−1 | 失敗による欠番なし。同じ実行の全CSVで共通 |
| arm / arm_a / arm_b | 腕0〜3 | trajectories.csvの中心行だけarm=−1 |
| endpoint_a / endpoint_b | endpoint_id = 4 × molecule_id + arm_id、0〜4N−1 | 分子ID = endpoint_id / 4、腕ID = endpoint_id % 4（整数演算） |
| site / site_a / site_b | site = (x × L + y) × L + z、0〜L³−1 | 格子点を1次元の番号にしたもの。現在0≤x,y,z<Lで一意 |
| bond_id | BondStoreの保存順に0〜E−1 | 行番号に対応。物理的な永続IDではない |
| cluster_id | クラスタ内の最小分子IDが小さい順に0〜C−1 | UF内部のroot番号とは別。孤立分子にも付ける |

siteからの復元は`z=site%L`、`y=(site/L)%L`、`x=site/(L*L)`（整数除算）。
分子IDからクラスタを引くにはcluster-members.csvを使い、その分子IDでtrajectories.csvの形状を参照する。
クラスタ単位に軌跡を複製した別ファイルは作らない。結合の詳細はbonds.csvを参照する。
IDは実行内で共通だが、seedや条件を変えた別実行の同じIDが同じ物理的分子・クラスタを表すとは限らない。

## 11. trajectories.csv：分子の中心と腕の軌跡

```csv
molecule,type,arm,step,site,x,y,z,direction
```

| 列名 | 型・範囲 | 意味・値の求め方 |
|---|---|---|
| molecule | 非負整数 | 分子ID。成功順に採番 |
| type | A / B | 分子の型。4本の腕と中心は同じ型 |
| arm | −1、0〜3 | 中心なら−1、腕の到達点なら腕番号 |
| step | 0〜l | 中心は0。各腕の到達点は1〜l。step=lが末端 |
| site | 非負整数 | その点の格子点ID。式は10.2参照 |
| x | 整数0〜L−1 | 絶対座標のx成分 |
| y | 整数0〜L−1 | 絶対座標のy成分 |
| z | 整数0〜L−1 | 絶対座標のz成分 |
| direction | −1、0〜5 | 直前の点からこの点へ進んだ方向。中心のみ−1 |

方向番号は0=+x、1=−x、2=+y、3=−y、4=+z、5=−z。
各腕のstep=1の直前は分子中心。step≥2の直前は同じ腕のstep−1。
中心は1分子につき1行だけ出力する。4腕の共通起点を4重に数えない。
データ行数は`N*(1+4*l)`、ヘッダを含む行数はこれに1を加える。
CSVのstepは1始まりの到達番号、Cのpath_sites添字stepは0始まりである点に注意する。
CSVは絶対座標を保存する。1分子図の「中心相対座標」はPythonで中心座標を差し引いた表示で、CSVそのものの座標ではない。

## 12. placement-summary.csv：生成目標と実績

```csv
input_mode,requested_phi,target_molecules,resolved_phi,placed_molecules,occupied_sites,actual_phi,type_A,type_B,max_placement_attempts,attempts,center_occupied,trapped,completed,stop_reason
```

ヘッダの次に実行全体の集計を1行出力する。候補探索・結合・UFより前に保存する。

| 列名 | 型・値 | 意味・計算方法 |
|---|---|---|
| input_mode | molecule_count / occupied_fraction | 分子数直接指定か、占有率指定か |
| requested_phi | 実数、または空欄 | target_phiとして入力した占有率。N直接指定なら未指定なので空欄 |
| target_molecules | 正整数 | 読込後の生成目標N_target。占有率指定ならfloor(requested_phi*L³/(1+4l)) |
| resolved_phi | 実数 | N_target*(1+4l)/L³。目標Nが全て完成した場合の占有率 |
| placed_molecules | 非負整数 | 実際に完成して登録された分子数N |
| occupied_sites | 非負整数 | 実際の格子ownerを走査して数えた占有点数。通常N*(1+4l)と一致 |
| actual_phi | 実数 | occupied_sites/L³。実際の占有率。0完成なら0 |
| type_A | 非負整数 | 完成したA型分子数N_A |
| type_B | 非負整数 | 完成したB型分子数N_B。N_A+N_B=N |
| max_placement_attempts | 正整数 | CLIで指定した全分子合計の試行上限。既定10000 |
| attempts | 非負整数 | 中心候補を選んで1分子生成を試みた総回数。占有中心・trappedも含む |
| center_occupied | 非負整数 | 選んだ中心候補が既に占有されていた棄却回数 |
| trapped | 非負整数 | 空中心から成長を始め、4腕が完成する前に行き止まりになった回数 |
| completed | 0 / 1 | 生成目標を達成したら1。後続の解析・描画の成否までは表さない |
| stop_reason | target_reached / attempt_limit / generation_error | 生成ループの停止理由。目標達成／試行上限／生成の実行エラー |

`floor`はその数以下の最大の整数。今回の正の値では小数部切り捨てで、目標占有率を超えないN_targetを選ぶ。
0分子へ換算されるtarget_phiは入力エラー。黙って1分子へ増やさない。
換算は十進入力の桁演算で扱い、表示用doubleの丸めによる分子数のずれを防ぐ。
表示ではrequested_phiが`0.10000000000000001`となる場合があるが、入力の十進原文はconfig.cfgに残る。
requested_phi・resolved_phi・actual_phiは別の値。未達ならactual_phiはresolved_phiより低くなる。
正常達成・上限未達では`attempts = placed_molecules + center_occupied + trapped`。
メモリ不足等の実行エラーが起きた場合、最後のエラー試行はこの3分類に含まれないことがある。
completed=1の後に描画等が失敗してrun.jsonがfailedになる場合もある。実行全体の成否はrun.jsonで確認する。
A型の入力確率rはconfig.cfgのa_probability。実現A比率はtype_A/placed_moleculesで別途算出し、0完成では未定義。

## 13. candidates.csv：結合候補の末端ペア

```csv
endpoint_a,molecule_a,arm_a,type_a,site_a,x_a,y_a,z_a,endpoint_b,molecule_b,arm_b,type_b,site_b,x_b,y_b,z_b,manhattan_distance
```

| 列名 | 意味・計算方法 |
|---|---|
| endpoint_a | ペアの先頭側の末端ID |
| molecule_a | endpoint_a/4（整数除算）で得る所属分子ID |
| arm_a | endpoint_a%4で得る腕番号 |
| type_a | molecule_aのA/B型 |
| site_a | molecule_a・arm_aの最後の軌跡点のsite |
| x_a | site_aの絶対x座標 |
| y_a | site_aの絶対y座標 |
| z_a | site_aの絶対z座標 |
| endpoint_b | ペアの後側の末端ID |
| molecule_b | endpoint_b/4（整数除算）で得る所属分子ID |
| arm_b | endpoint_b%4で得る腕番号 |
| type_b | molecule_bのA/B型 |
| site_b | molecule_b・arm_bの最後の軌跡点のsite |
| x_b | site_bの絶対x座標 |
| y_b | site_bの絶対y座標 |
| z_b | site_bの絶対z座標 |
| manhattan_distance | abs(x_a−x_b)+abs(y_a−y_b)+abs(z_a−z_b)。現在1または2 |

ID・siteは非負整数、座標は0〜L−1の整数、typeは文字A/B、距離は正整数。
重複除去のためendpoint_a<endpoint_bで保存する。**接尾辞a/bはペアの左右であり、A型/B型という意味ではない。type_aがB、type_bがAの行もある。**
同一分子のペアは除き、異種型で距離1〜2のペアを全て候補として保存する。
候補段階では同じ末端が複数行に現れてよい。まだ結合していないため、末端高々1回の制約はbonds.csvに適用する。
候補ID列は作らない。1行が1組で、候補の並びは抽選順位や結合確率そのものを表さない。
距離列はCSV出力時に座標から計算する。探索本体は24個の座標差とsite索引を使い、全末端対の距離を計算しない。

## 14. bonds.csv：成立した結合

```csv
bond_id,endpoint_a,molecule_a,arm_a,type_a,site_a,x_a,y_a,z_a,endpoint_b,molecule_b,arm_b,type_b,site_b,x_b,y_b,z_b,manhattan_distance
```

| 列名 | 型・意味 |
|---|---|
| bond_id | 非負整数。成立結合の保存順に0〜E−1 |
| endpoint_a | 非負整数。先頭側の末端ID |
| molecule_a | 非負整数。endpoint_a/4 |
| arm_a | 整数0〜3。endpoint_a%4 |
| type_a | 先頭側分子の型A/B |
| site_a | 先頭側の末端site |
| x_a | 先頭側末端の絶対x座標 |
| y_a | 先頭側末端の絶対y座標 |
| z_a | 先頭側末端の絶対z座標 |
| endpoint_b | 非負整数。後側の末端ID |
| molecule_b | 非負整数。endpoint_b/4 |
| arm_b | 整数0〜3。endpoint_b%4 |
| type_b | 後側分子の型A/B |
| site_b | 後側の末端site |
| x_b | 後側末端の絶対x座標 |
| y_b | 後側末端の絶対y座標 |
| z_b | 後側末端の絶対z座標 |
| manhattan_distance | 正整数。座標から計算したマンハッタン距離。現在1または2 |

bond_id以外の列の定義はcandidates.csvと共通。成立ペアは候補の部分集合で、endpoint_a<endpoint_b。
同じ末端IDはendpoint_a列とendpoint_b列を合わせて高々1回だけ現れる。
同じ分子対を別の腕で複数回結んだ場合は別行・別bond_idで保存する。
候補7組から結合3本、のように件数は一致しない。現在の結合処理は最大結合数を保証するものではない。
緑の結合線は2末端の接続表示であり、その間に新しい格子点を占有する腕の軌跡を追加したものではない。
結合前後でoccupied_sitesは変わらない。

## 15. cluster-members.csv：分子からクラスタへの対応

```csv
molecule,cluster_id,cluster_size,bond_degree
```

| 列名 | 型・意味・求め方 |
|---|---|
| molecule | 非負整数。分子ID。0〜N−1の順に1行ずつ保存 |
| cluster_id | 非負整数。結合を使ったUFの連結成分に、10.2の規則で付け直した出力ID |
| cluster_size | 正整数。そのクラスタに属する分子数。格子点数ではない |
| bond_degree | 整数0〜4。その分子につながる成立結合の本数＝結合済み腕の数 |

データ行数はN。孤立分子はcluster_size=1、bond_degree=0で保存する。
bond_degreeは異なる隣接分子の数ではない。同じ相手に2腕を結べばbond_degreeは2になる。
同じcluster_idの分子を集め、trajectories.csvのmolecule列で絞ればクラスタ形状を復元できる。

## 16. cluster-summary.csv：クラスタごとの集計

```csv
cluster_id,molecule_count,bond_count
```

| 列名 | 型・意味・求め方 |
|---|---|
| cluster_id | 非負整数。cluster-members.csvと共通の出力クラスタID |
| molecule_count | 正整数。そのクラスタ内の分子数 |
| bond_count | 非負整数。そのクラスタ内の成立結合の本数。多重結合も1本ずつ数える |

データ行数はC。孤立分子もmolecule_count=1、bond_count=0として1行保存する。
連結している2分子間に2本の結合があってもmolecule_countは2、bond_countは2である。
クラスタ内の独立閉路数を後で求めるなら`bond_count−molecule_count+1`（連結多重グラフ）。この値は現在CSVの独立列としては出力しない。

## 17. network-summary.csv：ネットワーク全体の集計

```csv
molecules,bonds,clusters,largest,second_largest,isolated_molecules,largest_fraction
```

| 列名 | 型・意味・求め方 |
|---|---|
| molecules | 非負整数。集計対象の完成分子数N |
| bonds | 非負整数。成立結合数E。bonds.csvのデータ行数 |
| clusters | 非負整数。孤立分子を含むクラスタ数C |
| largest | 非負整数。最大クラスタに属する分子数N_1 |
| second_largest | 非負整数。2番目のクラスタの分子数。同じサイズのクラスタが複数あれば同値を許す。1クラスタしかなければ0 |
| isolated_molecules | 非負整数。bond_degree=0の分子数 |
| largest_fraction | 実数、または空欄。largest/molecules。0分子の空ネットワークでは未定義なので空欄 |

ヘッダの次に全体の集計を1行保存する。largestとsecond_largestはクラスタIDではなく分子数。
現在のnetworkモードは0完成の未達時にはUFを実行しないため、このCSV自体を作らない。
0分子でlargest_fractionを空欄にする規則は、clusterモジュールの汎用出力仕様である。
先生の資料の`S_1 = m*N_1/L³`と、このCSVの`largest_fraction = N_1/N`は分母が異なるため別の量。
両者は`m=1+4l`として`資料のS_1 = actual_phi * largest_fraction`で対応する（N>0）。

## 18. ファイル間の整合性と確認例

正常完了時には、同じ実行フォルダ内で以下を確認できる。

| 確認項目 | 成立する関係 |
|---|---|
| 完成分子数 | placement-summary.placed_molecules = network-summary.molecules = cluster-membersの行数 |
| 占有点数 | occupied_sites = N*(1+4l) = trajectoriesのデータ行数（中心を1行で保存） |
| 型の集計 | type_A + type_B = N |
| 末端数 | 4N。trajectoriesのstep=lの行数 |
| 結合の集計 | network-summary.bonds = bondsの行数 = cluster-summary.bond_countの総和 |
| 結合次数 | cluster-members.bond_degreeの総和 = 2E |
| 分子の所属 | cluster-summary.molecule_countの総和 = N |
| クラスタ数 | network-summary.clusters = cluster-summaryの行数 |
| 最大クラスタ | largest = cluster-summary.molecule_countの最大値 |
| 孤立分子 | isolated_molecules = bond_degreeが0の分子数 |

ここで「行数」はヘッダを除くデータ行数。同一siteを各腕に重複出力する形式ではないので、軌跡行数と占有点数が一致する。
ID列の対応でファイルを結合する。行位置だけで別CSV同士を対応させない。

### 18.1 L16・腕長5・r0.5・phi0.1・seed3の実行例

目標N_target=19、1分子21点、入力占有率0.1、換算後占有率0.097412109375。
上限10000回の実行ではN=19、N_A=11、N_B=8、occupied_sites=399、attempts=23、center_occupied=2、trapped=2。
候補7組、結合3本、16クラスタ（2分子クラスタ3個＋孤立分子13個）。以下はnetwork-summary.csvの実際の形式。

```csv
molecules,bonds,clusters,largest,second_largest,isolated_molecules,largest_fraction
19,3,16,2,2,13,0.10526315789473684
```

同じ設定で上限を1回にするとN=1、occupied_sites=21、actual_phi=0.005126953125、completed=0、stop_reason=attempt_limit。
未達時は配置集計と完成軌跡だけを残し、候補・結合・クラスタのCSVは作らない。
この例は1seedの動作確認であり、ゲル化点や統計的な結論ではない。

## 19. CSV以外の保存ファイルと再現

| ファイル | 内容・用途 |
|---|---|
| config.cfg | 実行時の入力cfgの原文コピー。元の設定変更後もこの実行の入力を残す |
| run.json | mode、status、開始・終了時刻、コマンド、目標、summary、placement、reported_conditions、files、終了コード、Git情報、ファイルのSHA-256など |
| metadata.log | Cが標準エラーに出す固定条件・生成・探索・結合・クラスタの集計。L、seed、腕長、A確率、境界、成長、距離、結合規則も記録 |
| molecules.png | CSVからPythonで描画。A/B型と結合を表示。全体の実空間図。代表最大成分の抽出・描画はgraph/の追加出力を参照 |
| source.zip | srcの.c、includeの.h、scriptsの.py、Makefile、可視化requirements.txtの控え。build・data・仮想環境・tests・docsは含めない |
| build.log | ビルドコマンドの出力とエラー記録 |
| plot.log | 描画コマンドの出力とエラー記録。描画を実行した場合のみ作成 |

run.jsonのnetworkモードのschema_versionは2。既存モードは1。
run.jsonのplacementはplacement-summary.csvを読み込んだ値で、数値も文字列として格納される。
summaryのplaced_molecules、type_counts、occupied_sites、trajectory_rowsは軌跡CSVから集計したJSON数値。
generation_completedとanalysis_completedはnetworkモードの生成・解析完了の真偽値、plot_generatedは描画完了の真偽値。
描画だけの失敗でもstatusはfailedになるため、生成・解析の各完了フラグとログを併せて確認する。
早い段階の失敗では完了フラグがまだ存在しない場合もある。存在しない値を成功や0として補完しない。
sha256は各成果物、binary_sha256は実行バイナリの指紋。run.json自身の自己ハッシュは含めない。
source.zipは実行時の作業コピーを記録するので、未コミット変更がある場合もGitコミットIDだけに依存しない。
再現にはconfig.cfg、試行上限を含む実行コマンド、実装と乱数消費順が必要。

## 20. 占有率換算式の根拠

参考：水高将吾「Tetra-PEG型SAWネットワークのゲル化転移」、ERATO複素ゲルプロジェクト公開シンポジウム、2026年7月22日。提供PDF「ERATO公開シンポ.pdf」の5・6ページ。
資料では`phi = m*n/L³`、例としてf=4、l=5、m=21点/分子を示す。
現在のモデルでは`m=1+4l`なので、分子数について解くと`n=phi*L³/(1+4l)`。
目標占有率を超えない整数分子数にするため、今回の実装で`N_target=floor(phi*L³/(1+4l))`を採用した。
占有率の定義は資料に基づく。floorによる整数化は実装上の選択であり、資料にその丸め指定が明記されていると主張しない。

## 21. CSV仕様の正本と更新時の照合先

| 仕様 | 照合先 |
|---|---|
| 軌跡の列と中心行 | src/molecule.c：molecule_write_csv |
| 配置実績の列 | src/demo_network.c：placement-summary保存処理 |
| 候補の列と末端の復元 | src/endpoint_search.c：candidate_write_csv |
| 結合の列 | src/bond.c：bond_write_csv |
| クラスタの3 CSV | src/cluster.c：各write_csv関数 |
| ファイル構成と状態管理 | scripts/run_demo.py |

コードで列名や出力条件を変更した際は、この設計書のヘッダ・表・例も同時に更新する。
現段階の未出力指標（独立閉路数、クラスタ別色分け等）は、既存CSV列と混ぜず追加機能として扱う。

## 22. 2026-10-15研究報告後の方針と完成範囲

基本機能（cfg→多分子配置→候補探索→A–B結合→UF→7 CSV・実空間図）は完成済み。
これは1試料の生成・保存機能の完成であり、ゲル化の実証・モデル条件の確定ではない。
UFは所属・サイズの集計、BFSは独立したテストの照合に使用する。
`tests/test_cluster.c`の生成試験は3D、L=8、6分子、腕長5、A/B各3分子、seed=0..49。
各seedで同分子対の複数結合を許可・禁止の両方についてUFとBFSを比較する。
2分子デモの集計とは別の試験。BFSの通常解析APIや最短距離解析は未実装。

| 優先順 | 追加する内容 | 完了判断 |
|---|---|---|
| 1（今回） | 分子グラフの枝リスト、全体図、最大コンポーネントの抽出と実空間・グラフ図 | 元CSVとID・枝本数・所属が一致、多重結合と孤立分子を保持 |
| 2 | 測定量の定義確定 | 先生の資料と記号・分母・平均方法を対応づける |
| 3 | phiを変える複数seed実験 | 平均・誤差・未達率を記録し、未達を目標phiの試料に混ぜない |
| 4 | Lを変える比較 | 有限サイズ効果を調べ、1試料の最大比だけで転移を断定しない |
| 5 | 同一分子対の1/2/3/4重結合の頻度 | 単純なペア別結合本数から始める |
| 6 | 文献の網目分類・サイクル解析 | 必要性と定義を確認して追加 |

周期境界、電気伝導・抵抗、ラプラシアン、一般サイクル列挙、UF/BFSの解析機能拡張は後段。
当面は先生のメモと資料に沿う単純な測定を優先し、独自指標を先行導入しない。

### 22.1 測定量：定義済みと未確定を分離

N=完成分子数、N1=最大クラスタ分子数、m=1+4l、M=成立した枝の総本数、C=孤立分子を含むクラスタ数。
Mは既存コード・CSVではEまたはbondsに対応する。文献ごとに同じ記号の意味が異なるので読み替える。

| 測定量 | 定義・方針 | 現在の出力 |
|---|---|---|
| 実占有率 phi | m*N/L³。入力目標phiとは区別 | actual_phi |
| 最大分子数 N1 | 最大クラスタの分子数 | largest |
| 最大クラスタ比 Pmax | N1/N。N>0 | largest_fraction |
| 資料のS1 | m*N1/L³ = actual_phi*Pmax。資料PDF5・6ページ | 既存CSVから算出可、専用CSV列は未追加 |
| 枝本数 M | 成立結合1本につき1、多重結合を保持 | bonds |
| 末端結合率 p | 2M/(4N)。結合本数ではなく結合済み末端の割合 | 既存CSVから算出可 |
| 平均クラスタサイズ | 算術平均か重み付き平均か、最大を除くか、同率最大をどう扱うかを確定する | 未追加 |
| ペア別結合本数 | 各分子ペアに何本の枝があるか。1〜4本の頻度 | bonds.csvから算出可、頻度表は未追加 |
| M1、m1=M1/M | M1は最小網目の個数という仮理解。三重結合などの数え方と分母を先生に確認 | 未確定・未実装 |
| 独立閉路数 mu | M-N+C。全単純閉路数とは別 | 追加候補、今回は未追加 |

phiを横軸に最大成分指標・p・Mを描き、L別に比較する。Mは示量量なのでpも併記して系サイズを比較する。
L³は当面実験条件であり、縦軸の指定と断定しない。各条件でseedを変え、統計誤差と有限サイズ効果を評価する。
M1を「網目の枝本数の割合」と読み替えない。組合せ数、独立閉路数、局所構造分類を混同しない。

### 22.2 後段の参考文献

以下は測定定義の比較対象。今回の実装の生成規則やM1の定義の根拠と断定しない。
- Schwenke, Lang, Sommer, “On the Structure of Star-Polymer Networks”, https://arxiv.org/abs/2103.16259 （図1.1・第5節の多重結合構造）。
- Lang, Schwenke, Sommer, “Short Cyclic Structures in Polymer Model Networks: A Test of Mean Field Approximation by Monte Carlo Simulations”, Macromolecules 45, 4886–4895 (2012), DOI:10.1021/ma300317z。
- Lang, Müller, “An Analysis of the Gel Point of Polymer Model Networks by Computer Simulations”, Macromolecules 53, 498–512 (2020), DOI:10.1021/acs.macromol.9b02217 （式3のcycle rank）。
分子の運動・反応規則・濃度定義が現モデルと異なるので、数値的一致を無条件に期待しない。

## 23. 分子グラフと最大コンポーネント（実装済み）

Cの生成・結合・UFと既存7 CSVは変更しない。保存済みCSVをPythonで後処理し、新しい乱数を消費しない。
分子中心を頂点、成立結合を無向枝とする多重グラフ。候補CSVを成立枝として扱わない。
N×Nの隣接行列や全サイクル列挙は導入せず、枝リストと分子別の隣接リストを使う。
追加依存は設けず、標準ライブラリと既存Matplotlibを利用する。

### 23.1 最大成分の選択

cluster-members.csvのUF所属を利用し、結合からの独立した探索でも所属・サイズ・次数を照合する。
最大サイズが同率の場合はcluster_idが最小の成分を代表として1個抽出する。
同率最大のcluster_idをすべてgraph-analysis.jsonに保存し、任意の代表選択を隠さない。
全孤立なら分子0を含む成分を代表とし、結合0本のグラフ・形状図を保存する。
空試料・未達試料に結合解析を勝手に追加しない。run.jsonの生成状態を確認する。
IDは元実行と同じ値を保持し、部分抽出後も0始まりに振り直さない。

### 23.2 出力設計

追加成果物は実行フォルダのgraph/配下にまとめ、元のCSVを上書きしない。
- molecular-edges.csv: bond_id,molecule_a,molecule_b,endpoint_a,endpoint_b,cluster_id。成立結合1本1行。多重結合も各1行。
- largest-members.csv: molecule,cluster_id,cluster_size,bond_degree。代表最大成分のみ。列の意味は15節と同じ。
- largest-edges.csv: molecular-edgesと同じ列。代表最大成分の内部枝のみ。
- graph-analysis.json: 全分子数、枝数、クラスタ数、代表cluster_id、同率最大ID、最大分子数、最大枝数、正規化量、描画配置規則。
- network.png: 全分子グラフ。A/Bの列にID順で配置し、多重結合は曲線で分離する。座標軸は実空間距離を表さない。
- largest-network.png: 同じ規則で代表最大成分を描画。
- largest-molecules.png: 代表最大成分を絶対座標で描画。元軌跡・成立枝の位置を利用する。

グラフ図の位置は決定的な表示用配置であり、SAWの実空間位置・距離・物理的な網目面積を表さない。
大規模グラフの可読性は保証せず、IDラベルは小規模時だけ付ける。
--no-plotでは新しいCSV/JSONだけ出力しMatplotlibを要求しない。描画結果の保存有無はJSONで区別する。
通常network/clusters実行の完了時に後処理を接続し、既存実行への単独コマンドも用意する。
検証：三重結合、孤立、同率最大、4分子サイクル、破損所属の拒否、既存19分子実行、未達時の非実行。


## 24. 派生グラフ出力の仕様と実行方法（実装済み）

入力は同じ実行のtrajectories.csv、bonds.csv、cluster-members.csv、network-summary.csv、run.json。
read_moleculesはCSVを1回読み、既存の分子ごとの軌跡検査を再利用する。従来の分子ごとの全CSV再読込を除去。
所属の検証用BFSは通常後処理で1回実行するが、UFの集計実装は置換しない。
グラフ処理の追加メモリ・探索はO(N+M)、軌跡読込はO(N*l)、ID・stepの整列費用は別。
可視化のためにN×N隣接行列や全頂点対距離、乱数によるレイアウトを導入しない。

### 24.1 molecular-edges.csv / largest-edges.csv

```csv
bond_id,molecule_a,molecule_b,endpoint_a,endpoint_b,cluster_id
```

| 列 | 意味 |
|---|---|
| bond_id | 元bonds.csvの成立結合ID。抽出後も採番し直さない |
| molecule_a | 元の先頭末端の所属分子ID。A型とは限らない |
| molecule_b | 元の後側末端の所属分子ID。B型とは限らない |
| endpoint_a | 元の先頭末端ID、4*molecule_a+arm_a |
| endpoint_b | 元の後側末端ID、4*molecule_b+arm_b |
| cluster_id | 両分子共通のUF所属ID。異なる所属ならエラー |

すべて非負整数。全体枝リストはM行、最大成分枝リストは代表成分の枝数行。
元bonds.csvの順序・IDを保持する。多重結合を1枝へ縮約しない。結合0本ならヘッダだけ。
孤立分子は枝リストに現れないため、頂点一覧の正本であるcluster-members.csvも併用する。
largest-members.csvの4列は15節と同義で、代表最大成分に所属する分子行だけを保持する。

### 24.2 graph-analysis.json

| キー | 型・意味 |
|---|---|
| schema_version | 整数1。この派生JSONの形式バージョン |
| molecules / bonds / clusters / largest | 整数。全体のN / M / C / 最大分子数N1 |
| selected_cluster_id | 非負整数。代表最大成分の元ID |
| tied_largest_cluster_ids | 同率最大成分すべてのIDの昇順配列。1個の場合も配列 |
| selection_rule | smallest_cluster_id。代表の選択規則 |
| largest_bonds | 整数。代表最大成分の枝数 |
| largest_fraction | N1/N（全体分子数を分母にする） |
| actual_phi | occupied_sites/L³。run.jsonにLがない旧実行ではnull |
| largest_occupied_fraction | actual_phi*largest_fraction。資料のS1に対応。L不明ならnull |
| layout | deterministic_type_columns_id_order。A/B二列、各列は元分子ID順 |
| spatial_coordinates | original_absolute_open_boundary。実空間図は元の絶対座標 |
| validation | bond_constraints_and_bfs_partition_and_degrees。検証規則の識別子 |
| plots_generated | 3種類の図をすべて保存できたらtrue。--no-plotはfalse |
| files | CSV3種類と、正常に全描画完了した場合のPNG3種類の相対ファイル名 |
| status | completed / failed。後処理の成否。途中失敗時はフォルダの診断用ファイルも確認 |

run.jsonにgraph_analysis_completedとgraph/の各出力・ハッシュを追加する。
既存run.jsonのschema_versionおよび既存7 CSVの列・物理的生成結果は変更しない。
後処理失敗時は元CSVを保持し実行全体をfailedにする。生成完了フラグと混同しない。
単独後処理は元run.jsonを書き換えず、graph-analysis.jsonで成否を記録する。
未達・0完成・failed実行は単独後処理でも拒否する。

### 24.3 コマンド

```sh
# 新しい実行：既存図に加えgraph/を自動保存
python scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg
# 描画なし：既存CSVと派生グラフCSV/JSONを保存。Matplotlib不要
python scripts/run_demo.py --mode network --config configs/network/preview_phi010_L16_seed3.cfg --no-plot
# 古いcompleted実行を後処理：実行フォルダを指定
python scripts/visualization/analyze_graph.py data/runs/<実行フォルダ>
# 既存graph/がある場合、元を残して別フォルダへ出力
python scripts/visualization/analyze_graph.py data/runs/<実行フォルダ> --output-dir data/runs/<実行フォルダ>/graph-v2
```

出力先が既に存在する場合は上書きせずエラーにする。dpiは正整数、既定200。
CLIの--no-plotと--dpiは単独後処理でも指定可能。
今回のseed3・L16・19分子では最大2分子の成分が3個あり、代表はcluster_id=2、分子ID=2と5、内部枝1本。
全体図の線の交差は接続を表さない。頂点として表示された分子だけが接続点。
この抽出例は転移の証拠ではなく、描画・ID対応・保存の動作確認である。


## 25. CSV読込み式の対話的3Dビューア（2026-10-15）

### 目的と実装範囲

静止画に加えて、同じ実行のtrajectories.csvとbonds.csv（任意）をブラウザで読み込み、
実空間形状を回転・拡大・平行移動して確認する。Cの生成・乱数・CSV仕様は変更しない。
ソースはweb/viewer/、実行方法と制限は同ディレクトリのREADME.mdを正本とする。

- 外部描画ライブラリを使わずWebGLの線と中心点をバッファへまとめる。視点・表示変更時だけ再描画。
- CSV処理はModule Workerで行い、ファイルをサーバへ送信しない。
- 基本検査と詳細検証を分離。通常は型・ID・腕長・末端重複使用を検査して表示を優先し、最近接性・占有点重複・CSV間の座標・結合距離は別ボタンで検査。
- 結合から連結成分を再集計する。多重枝を保持。同率最大は所属最小分子IDが最小の1成分を選ぶ。
- 軌跡のみなら形状を表示するが、最大成分・成分別着色・結合数は未判定として表示する。
- A/B別表示、腕・結合切替、最大成分のみ、成分別着色、直交視点、PNG保存。
- config.cfgを読まないので、占有点の最大座標からLやphiを推定しない。境界箱も描かない。

### 負荷と表示の意味

CSV合計32 MiB・軌跡50万行・結合20万行を上限とし、読込み60秒で中断する。
最大描画分子数1,000／5,000／20,000を選択でき、超過時は分子ID順の等間隔間引き。
腕は分子単位で保持し、結合は両端の分子が表示される場合のみ描画する。
統計値・成分所属は全入力で計算し、表示数と間引き状態を画面・PNGに明記する。
この上限は実機性能で保証された最大値ではなく、無制限な読込みを避ける初期設定である。
ブラウザ読込みは文字列・オブジェクト・Worker間転送でもメモリを使う。Cの消費量と別に評価する。

### 検証と今後

make test-viewerで、既存20分子試料・二重結合・孤立・同率最大・表示間引き・不正CSVを確認。
Node.jsはこの検証時のみ必要で、生成処理・ビューア利用自体には不要。
PC上での操作・GPU性能・大規模入力は別途実機確認する。
性能ベンチマークと保存モードの設計は次の課題とし、この変更では自動実験や保存省略を追加しない。


### 25.1 読込み状態・代替経路・表示専用形式

2026-10-15追記。19分子399点・3結合（添付CSV計9,297バイト）で解析と描画用データ生成は
Node.js上で約8 msだった。長時間待ちはこの入力のサイズでは説明できないが、
ブラウザ実行環境の根本原因は未特定。処理段階・エラーを表示して切り分け可能にする。

- Workerはready応答後にファイルを受け取り、段階の進捗・結果・エラーを別メッセージで返す。
- 2.5秒で起動応答なし、またはWorker起動失敗時、合計2 MiB以下なら同じパーサを直接呼ぶ。
- 大きな入力でWorkerを使えなければ案内を表示。中止・再選択時は以前の結果を混入させない。
- 表示の初期状態は詳細未検証。詳細検証は元の選択ファイルから再実行し、視点を保つ。
- 表示専用形式ab-saw-viewer-1は座標・型・ID・末端間結合を配列で保存。集計と成分は再計算。
- ブラウザからの保存、およびscripts/export_viewer.pyによるCSV後処理に対応。
- 元データの条件・seed・候補・ログを省略した派生物なので、研究データの正本として置き換えない。
- 添付試料は表示専用形式で3,366バイト。サイズ縮小は試料に依存し、実ブラウザの速度向上は未計測。
- Worker起動失敗・無応答・中止・正常なメッセージ処理、CSVと表示専用形式の往復を回帰検証。


### 25.2 読込み後に空画面となる場合の処理

2026-10-15追記。「表示準備完了」のまま表示件数・ファイル名が更新されない報告を受け、
結果コールバックの例外が完了済みフラグにより握りつぶされる経路を修正した。
報告環境で発生した元の例外自体は未再現なので、特定済みと混同しない。

- 既定描画をCanvasの互換投影へ変更。回転・ズーム・移動・PNG保存は同じカメラ座標で行う。
- WebGLはrenderer=webgl指定で選択可。初期化・バッファ確保・描画・コンテキスト消失で失敗したらCanvasへ切替。
- 読込み後の結果処理と非同期描画の例外を画面へ出し、失敗を成功のままにしない。
- 表示成功は描画処理終了時に更新する。
- 空表示でもNaNを渡さない。カメラ投影と表示範囲を検査する。
- UI契約テストはサンプル読込み→データ受渡し→描画呼出し→表示件数と成功表示を確認。
- 添付19分子のデータをCanvasでラスター描画し、描画画素が存在することと実際の画像を確認。
- ブラウザ操作・実GPUでの動作はこの環境では未検証。


### 25.3 空画面の根本原因と回帰検証

2026-10-15追記。画面のエラー「Cannot set properties of undefined (setting disabled)」から、
表示範囲selectの2番目のoption開始タグが誤って閉じタグになっていることを特定した。
ブラウザ上では表示範囲selectに1選択肢しかなく、options[1].disabledで結果処理が停止していた。
正しいoption開始タグに修正した。CSV・GPU・ファイル容量が今回の原因ではない。

従来のUI疑似DOMがoptionsを固定で2個作っていたため、この誤りを検出できなかった。
テストを実HTMLから開始タグ・valueを読む方式へ変更し、修正前に報告と同じ例外を再現した。
また、ランタイム例外収集ハンドラをassert前に解除し、assert失敗が握りつぶされず非ゼロ終了するよう修正。
修正後はWebGL呼出し経路・Canvas経路ともに表示件数・選択肢・描画呼出し・成功表示を確認した。

## 26. 系サイズによる性能測定（設計：2026-10-15）

### 目的と条件

固定腕長・占有率・A型確率・境界条件でLを変え、PCで実行可能な範囲とサーバへの移行を判断する。
物理モデル、乱数系列、既存CSVの正本は変更しない。Linux/WSLの逐次実行を対象とする。
初期条件はconfigs/benchmarks/pc-pilot.json。L=16,32,64,128、phi=0.1,0.3、seed=3,17,42、腕長5。
各seedの反復数は設定可能。異なるseedの時間差には構造・棄却回数の差も含まれる。
初期設定は反復1回であり、短時間の精密な速度比較には同一seedの反復を増やす。

### 測定単位

- ビルドは開始時に1回だけ実行し、測定区間に含めない。バイナリSHA256、CC/CFLAGS/CPPFLAGS、build.logを保存。
- 占有率からの分子数換算と入力検査は既存Cのcheck_network_configを正本にする。
- C段階：分子生成・末端探索・結合・UF・7 CSVの書出しを含むプロセス全体。
- 後処理段階：独立Pythonプロセスで既存analyze_graphをno_plot=Trueで実行。CSV読込み・検証・成分抽出を含む。
- 描画段階：独立PythonプロセスでCSVを読み直し、既存描画関数で全実空間・全グラフ・最大グラフ・最大実空間の4 PNGを出力。
  プロセス全体に加え、内部の読込み時間と描画保存時間を分ける。描画時間にはMatplotlib初回importを含む。
- 最大RSSは各プロセスのLinux wait4 ru_maxrss。追加で/procを短周期監視する。段階ごとのRSSは加算しない。
  親PythonやビルドのRSSを混ぜない。小さいCの測定にはプロセス起動・exec前の高水位等の影響がある。
- gzip圧縮は別区間。ファイルごとの元バイト数・gzipバイト数・SHA256を保存し、復号したSHA256を照合する。
- ウォームアップなし、プロセス起動込み、OSページキャッシュを明示的に消さない。測定条件として記録する。

### 完了状態と停止条件

Cの終了0・配置completed=1を成功とし、終了3は目標未達、資源上限での停止・実行エラーと区別。
未達では後処理と描画を行わない。postprocess/plotの分子数上限による省略を明記する。
配置実績、実占有率、目標数、試行上限・実試行数・棄却数、最大成分も記録する。
試行上限はmax(min_attempts, attempts_per_molecule × Cが解決した目標数)。

経過時間、RSS、試料・全体の保存量、ディスク空き容量を監視。上限超過なら子プロセスを停止し、
未完了のファイルも保持して部分出力と識別し、以後の試料を開始しない。割込み時も子を回収する。
RSS・保存量の監視はサンプリングなので一時的な超過を完全には防がない。OSの厳密なハード制限ではない。
事前の格子・軌跡メモリ見積りは下限的な目安。候補・結合・UF・Pythonのオブジェクト領域等を別に使う。

### 保存と集計

data/benchmarks/<UTC+UUID>/へ、設定控え・環境・ソースZIP・実行情報・試料別結果・集計CSVを保存。
run.jsonはanalyze_graphとの互換性に必要な生成状態・条件を持たせる。通常のrun_demo出力とは入口を区別。
測定表は1試料×1段階の生値を保存し、成功した試料だけの中央値・最小・最大と失敗・未達数を集計。
ファイル別表は圧縮前と圧縮後のサイズを分ける。合計保存量にはPNG・ログ・JSON・ソースも含める。
retention=gzipでは照合後に当該実行が作ったCSVのみ圧縮版へ置換。元CSVのバイト列は復元可能。
retention=rawではCSVを残し、圧縮比較に使用したgzipを破棄する。既存データへの自動削除は行わない。
復元は明示コマンドで行い、空き容量を検査し、既存ファイルを上書きせずハッシュを照合する。

CPU・OS・Python・コンパイラ・実メモリ・cgroup/WSL情報を記録し、ここでの結果をユーザーPCの性能と混同しない。
初期結果の傾向からL=192,256へ広げ、必要な場合に384,512へ進む。32GBという搭載量だけで最大Lは決めない。
