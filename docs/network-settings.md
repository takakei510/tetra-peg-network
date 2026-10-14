# 次の①：多分子用設定の読込と検証

## 今回の範囲

多分子生成に向け、分子数・腕長・A型確率をファイルから読めるようにした。
①は設定の読込と入力検査。②として生成・候補探索・結合・UFへの接続も実装済み。
詳しくは[multi-molecule-network.md](multi-molecule-network.md)。専用チェッカー自体は生成しない。
既存のsingle/two/candidates/bonds/clustersデモは以前の固定条件を維持する。

## 2つの設定形式

| 形式 | 読み込み関数 | キー |
|---|---|---|
| 既存デモ用 | config_load | dim, L, seed |
| 多分子用 | network_config_load | 上記3項目＋arm_length, a_probability＋目標の片方 |

`NetworkConfig`構造体に既存の`Config base`を含め、3つの生成条件を追加した。
共通5項目は必須。目標は`target_molecules`または`target_phi`の片方だけ指定する。
暗黙の既定値や優先順位は入れず、両方指定・未知・重複・欠落キーはエラー。
従来のkey=value、#コメント、前後空白の書式を引き継ぐ。
古いデモへ新形式を渡すとエラーになる。追加設定を無視して固定条件で生成させない。

| キー | Cの参照先・型 | 意味・条件 |
|---|---|---|
| dim | config.base.dim : int | 現段階は3のみ |
| L | config.base.L : int | 1辺の格子点数。2以上、既存格子の点数上限内 |
| seed | config.base.seed : uint64_t | 乱数seed。0〜UINT64_MAX |
| target_molecules | config.target_molecules : size_t | 完成させたい分子数。正の整数、既存実装の範囲内 |
| target_phi | config.requested_phi : double | 入力する目標占有率。十進数で0より大きく1以下 |
| arm_length | config.arm_length : size_t | 各腕の歩数。中心を除き1以上 |
| a_probability | config.a_probability : double | A型を選ぶ入力確率r。有限な0〜1、B型は1-r |

A型の確率と実際のA型比率は別。r=0.5でもAとBの数を半数ずつ固定しない。
型の抽選位置は現在の分子生成関数と同じく空中心確認後を予定する。
今回は型の抽選は実行しない。

## 検査すること

整数は符号・小数・後続文字を含む場合に拒否。浮動小数はNaN/無限値や範囲外を拒否。
格子L^3、軌跡配列、1分子の仮占有リスト、探索回数の計算が既存の上限を超えないか検査。
4腕なので1分子の必要占有点数は `1+4*arm_length`。
`target_molecules*(1+4*arm_length) <= L^3` を満たさない設定は事前に拒否する。
乗算前に除算で上限を比較し、計算のオーバーフローを避ける。

この条件は必要条件で、配置成功・生成時間・メモリ確保成功を保証しない。
分子の形状や込み合いによるtrapped・配置試行上限は、②の生成時に扱う。
失敗した読込は呼出元の出力構造体を変更しない。

## 使用例

`configs/network/preview_N20_L16_seed3.cfg`：

```ini
dim = 3
L = 16
seed = 3
target_molecules = 20
arm_length = 5
a_probability = 0.5
```

```sh
make build/check_network_config
./build/check_network_config configs/network/preview_N20_L16_seed3.cfg
```

正常時は終了コード0、`config_validated=1`と設定・換算値を出力する。
例の必要占有点数は420点、目標の占有率は420/4096=0.1025390625。
`generation_run=0`を明記する。`resolved_phi`は分子数と腕長から計算した値で、生成実績ではない。
分子数指定では`input_mode=molecule_count requested_phi=not_set rounding=none`を出力する。
入力エラーは標準エラーに理由を出し、終了コード1。

`configs/network/`は多分子用、`configs/demos/`は従来デモ用、手元の変更は`configs/local/`。
例は読込確認用で、生成済み試料や統計実験ではない。
距離、境界、成長規則、最大試行数、多重結合条件は今回の新しい設定キーに含めない。
未対応のキーを追加するとエラーになる。

## 検証・既存コードとの関係

並び替え・空白・コメント、小数と指数表記、確率0/1、seed0/最大値を検証。
全必須キーの欠落、重複、未知キー、範囲外、数値の桁溢れ、長過ぎる行、
物理的に必要点数が収まらない設定、失敗時の出力不変、旧設定の互換性を検証。

直接使ったもの：`Config`型、`ARM_COUNT`、既存格子・軌跡配列の仕様と上限。
既存Config読込の書式と厳格な入力検査方針を引き継いで、新しい専用読込関数を実装した。
旧config_loadを呼んで追加項目を無視する方式ではない。既存生成やUFは今回は変更しない。
②ではtarget_moleculesを容量と生成目標、arm_lengthをMoleculeStore初期化、
a_probabilityをmolecule_try_generateの引数に接続する。

## 占有率から分子数への換算

生成側が使う値は常に`config.target_molecules`。占有率指定は設定読込時に分子数へ変換する。
4本の腕がすべて完成し、占有点が重複しないモデルなので、1分子は中心込みで`1+4*arm_length`点を占有する。

`N = floor(target_phi * L^3 / (1+4*arm_length))`

切り捨てを使い、入力占有率を超える分子数にはしない。0分子になる入力はエラーにし、黙って1分子へ増やさない。
十進の入力は整数の桁演算で`floor(target_phi*L^3)`を求め、さらに1分子分で整数除算する。
例えばL=10、腕長5、phi=0.105では正確に5分子になる。表示用doubleの丸めに左右されない。

`configs/network/preview_phi010_L16_seed3.cfg`では、N指定の行を`target_phi = 0.1`に置き換える。
専用チェッカーでこのファイルを指定すると、4096点の格子、1分子21点から19分子となる。
必要点数399、`resolved_phi=0.097412109375`。目標の0.1との差は整数分子数への切り捨てによる。

| 保存・出力する値 | 意味 |
|---|---|
| target_mode / input_mode | 分子数指定か占有率指定か |
| requested_phi | 入力した占有率。N指定の場合は出力not_set、構造体0 |
| phi_site_budget | floor(入力占有率×全格子点数)。phi指定時のみ使用 |
| target_molecules | 換算後、または直接指定した生成目標の分子数 |
| resolved_phi | 目標分子数がすべて完成した場合の占有率。実績ではない |

設定ファイル自体に十進の原文を残し、チェック出力には表示用doubleを17桁で出す。
`requested_phi=0.10000000000000001`のような表示はdoubleの表現であり、分子数の換算誤差ではない。

②では目標と実際の完成分子数・実占有率・未達状態・失敗回数を区別して保存する。
生成失敗時もそれまで完成した分子は残す。失敗した途中の形の履歴は保存しない方針。
この生成結果の保存はnetworkモードで実行する。設定チェッカー自体は生成しない。
