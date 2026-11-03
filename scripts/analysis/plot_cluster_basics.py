#!/usr/bin/env python3
"""Replot basic molecular-network observables from saved sample-level CSV.

No simulation is run. H2 columns are not used. Means of sample-level ratios
are calculated before averaging; missing finite-cluster means remain missing.
"""
import argparse
import csv
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import sys
from zipfile import ZipFile, ZIP_DEFLATED

REPO = Path(__file__).resolve().parents[2]
METRICS = ('N', 'N1', 'Pmax', 'S1', 'p', 'K', 'K_density', 'A', 'Sall', 'Sfin')
FIELDS = ('L', 'target_phi', 'actual_phi', 'N', 'N1', 'bonds', 'clusters',
          'cluster_number_mean', 'cluster_weighted_all', 'cluster_weighted_excl_one')


def checked_number(value, name, integer=False):
    if value is None or value == '':
        raise ValueError(f'{name}: 必須値がありません')
    number = float(value)
    if not math.isfinite(number) or (integer and number != int(number)):
        raise ValueError(f'{name}: 不正な数値 {value!r}')
    return int(number) if integer else number


def optional_number(value, name):
    return None if value is None or value == '' else checked_number(value, name)


def assert_close(a, b, label):
    if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11):
        raise ValueError(f'{label}: 保存値と定義が一致しません ({a} != {b})')


def parse_sample(row, sites_per_molecule=21, arms=4):
    L = checked_number(row['L'], 'L', True)
    N = checked_number(row['N'], 'N', True)
    N1 = checked_number(row['N1'], 'N1', True)
    K = checked_number(row['clusters'], 'clusters', True)
    bonds = checked_number(row['bonds'], 'bonds', True)
    if L <= 0 or N < 1 or not 1 <= N1 <= N or not 1 <= K <= N or bonds < 0:
        raise ValueError('L,N,N1,K,bondsの範囲が不正です')
    if 2*bonds > arms*N:
        raise ValueError('結合済み末端数が全末端数を超えています')
    actual_phi = checked_number(row['actual_phi'], 'actual_phi')
    assert_close(actual_phi, sites_per_molecule*N/L**3, 'actual_phi')
    if not 0 < actual_phi <= 1:
        raise ValueError('実占有率が(0,1]の範囲外です')
    pmax, s1, p = N1/N, sites_per_molecule*N1/L**3, 2*bonds/(arms*N)
    for name, value in [('Pmax', pmax), ('S1', s1), ('p', p)]:
        if row.get(name) not in (None, ''):
            assert_close(checked_number(row[name], name), value, name)
    A = checked_number(row['cluster_number_mean'], 'cluster_number_mean')
    Sall = checked_number(row['cluster_weighted_all'], 'cluster_weighted_all')
    assert_close(A, N/K, 'cluster_number_mean')
    if not 1 <= Sall <= N:
        raise ValueError('Sallの範囲が不正です')
    Sfin = optional_number(row['cluster_weighted_excl_one'], 'Sfin')
    if N == N1:
        if Sfin is not None:
            raise ValueError('全分子が最大クラスタに属する試料のSfinは未定義です')
    elif Sfin is None or not 1 <= Sfin <= N-N1:
        raise ValueError('残余クラスタがある試料のSfinが不正です')
    if Sfin is not None:
        assert_close(Sfin, (Sall*N-N1*N1)/(N-N1), 'Sfin')
    target = checked_number(row['target_phi'], 'target_phi')
    if not 0 < target <= 1:
        raise ValueError('target_phiの範囲が不正です')
    seed = row.get('seed', '')
    if not seed:
        raise ValueError('seedがありません。重複試料を判定できません')
    return dict(L=L, target_phi=target, actual_phi=actual_phi, seed=seed,
                case_id=row.get('case_id', ''), N=N, N1=N1, Pmax=pmax, S1=s1,
                p=p, K=K, K_density=K/L**3, A=A, Sall=Sall, Sfin=Sfin)


def load_groups(path, config):
    groups = defaultdict(list)
    seen = set()
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        missing = set(FIELDS)-set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'入力CSVの列が不足しています: {sorted(missing)}')
        for lineno, row in enumerate(reader, 2):
            try:
                sample = parse_sample(row, config['sites_per_molecule'], config['arms'])
            except (ValueError, KeyError) as error:
                raise ValueError(f'{path}:{lineno}: {error}') from error
            key = (sample['L'], sample['target_phi'], sample['seed'])
            if key in seen:
                raise ValueError(f'重複試料: {key}')
            seen.add(key)
            groups[key[:2]].append(sample)
    if not groups:
        raise ValueError('入力CSVが空です')
    expected_sizes = config.get('sizes')
    expected_phis = config.get('target_phis')
    if expected_sizes and expected_phis:
        expected = {(int(L), round(float(phi), 8)) for L in expected_sizes for phi in expected_phis}
        actual = {(L, round(phi, 8)) for L, phi in groups}
        if expected != actual:
            raise ValueError(f'条件集合が設定と異なります。不足={sorted(expected-actual)} 余分={sorted(actual-expected)}')
    n_expected = config.get('expected_samples_per_condition')
    if n_expected is not None:
        bad = {str(key): len(values) for key, values in groups.items() if len(values) != n_expected}
        if bad:
            raise ValueError(f'各条件{n_expected}試料を期待しましたが一致しません: {bad}')
    for samples in groups.values():
        samples.sort(key=lambda s: (s['case_id'], s['seed']))
    return dict(groups)


def describe(values):
    values = [v for v in values if v is not None]
    n = len(values)
    return {'n': n, 'mean': statistics.mean(values) if n else None,
            'sd': statistics.stdev(values) if n > 1 else None,
            'se': statistics.stdev(values)/math.sqrt(n) if n > 1 else None}


def summarize(groups):
    summary = []
    for (L, phi), samples in sorted(groups.items()):
        row = {'L': L, 'target_phi': phi, 'n_samples': len(samples),
               'actual_phi': statistics.mean(s['actual_phi'] for s in samples)}
        for metric in METRICS:
            row.update({f'{metric}_{stat}': value
                        for stat, value in describe(s[metric] for s in samples).items()})
        summary.append(row)
    return summary


def resolve_samples(source):
    if source.is_file():
        return source
    if (source/'EXPORT_INCOMPLETE.json').exists():
        raise ValueError('不完全な出力が指定されています。完了済み結果を指定してください')
    if (source/'samples.csv').is_file() and not (source/'manifest.json').is_file():
        # Explicit standalone CSV export, with scientific validation below.
        return source/'samples.csv'
    folders = [source] + list(source.glob('results_*'))
    valid = []
    for folder in folders:
        if not (folder/'samples.csv').is_file() or not (folder/'manifest.json').is_file():
            continue
        if (folder/'EXPORT_INCOMPLETE.json').exists():
            continue
        m = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        count = m.get('completed_samples', 0)
        if m.get('status') == 'completed' and count > 0 and count == m.get('requested_samples'):
            valid.append((count, folder.stat().st_mtime_ns, folder/'samples.csv'))
    if not valid:
        raise ValueError('完了済みの結果が見つかりません。--resultsでresults_*フォルダまたはsamples.csvを指定してください')
    return max(valid)[2]


def setup_font(explicit):
    from matplotlib import font_manager
    if explicit:
        font_manager.fontManager.addfont(str(explicit))
        return font_manager.FontProperties(fname=str(explicit)).get_name()
    names = {f.name for f in font_manager.fontManager.ttflist}
    for name in ['Noto Sans CJK JP', 'Noto Sans JP', 'IPAexGothic', 'Yu Gothic', 'Meiryo']:
        if name in names:
            return name
    for candidate in ['/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
                      '/mnt/c/Windows/Fonts/meiryo.ttc', 'C:/Windows/Fonts/meiryo.ttc']:
        if Path(candidate).is_file():
            font_manager.fontManager.addfont(candidate)
            return font_manager.FontProperties(fname=candidate).get_name()
    raise ValueError('日本語フォントが必要です。WSL: sudo apt install fonts-noto-cjk、または --font フォントファイル')


def draw(summary, groups, output, config, font_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    family = setup_font(font_path)
    plt.rcParams.update({'font.family': family, 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.unicode_minus': False, 'svg.fonttype': 'path'})
    sizes = sorted({r['L'] for r in summary})
    colors = {L: color for L, color in zip(sizes, ['#2873b8', '#bd681b', '#218477', '#8862ad'])}
    if len(sizes) > len(colors):
        raise ValueError('色設定は最大4サイズに対応しています')
    peaks = [max((r for r in summary if r['L'] == L and r['Sfin_mean'] is not None),
                 key=lambda r: r['Sfin_mean']) for L in sizes]
    fig_names = []

    def save(fig, name, title, note):
        fig.suptitle(title, fontsize=16, y=.98)
        fig.text(.08, .018, note, fontsize=9, color='#45576a', va='bottom')
        fig.tight_layout(rect=(0, .09, 1, .92), pad=1.5)
        for ext in config['formats']:
            if ext not in ('png', 'svg'):
                raise ValueError(f'未対応の出力形式: {ext}')
            destination = output/f'{name}.{ext}'
            temporary = output/f'.{name}.tmp.{ext}'
            fig.savefig(temporary, dpi=config['dpi'], bbox_inches='tight', pad_inches=.2)
            if ext == 'png':
                from PIL import Image
                with Image.open(temporary) as image:
                    image.verify()
            temporary.replace(destination)
        plt.close(fig)
        fig_names.append(name)

    def curve(ax, metric, ylabel, scale=1, log=False, xmetric=None, limits=None):
        for L in sizes:
            rows = sorted((r for r in summary if r['L'] == L and r[metric+'_mean'] is not None), key=lambda r: r['target_phi'])
            if limits:
                rows = [r for r in rows if limits[0] <= r['target_phi'] <= limits[1]]
            ax.errorbar([r[xmetric+'_mean'] if xmetric else r['actual_phi'] for r in rows],
                        [scale*r[metric+'_mean'] for r in rows],
                        xerr=[r[xmetric+'_se'] or 0 for r in rows] if xmetric else None,
                        yerr=[scale*(r[metric+'_se'] or 0) for r in rows], fmt='o',
                        markersize=3, capsize=2, linestyle='none', color=colors[L], label=f'L={L}')
        ax.set_xlabel('末端結合率 p の試料平均' if xmetric else '実占有率 φ')
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale('log')
        ax.grid(axis='y', alpha=.18)
        ax.legend(frameon=False)

    common = '点は条件別の試料平均、誤差棒はSE。各試料で量を計算してから平均。線・補間・フィットなし。'
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    curve(axes[0], 'N1', '最大クラスタの分子数 N1［分子］', log=True)
    curve(axes[1], 'Pmax', '最大クラスタ比 N1/N［%］', 100)
    save(fig, '01_largest_cluster', '最大クラスタの絶対サイズと全分子に占める割合', common+'\n左は縦軸対数。占有率やLによって全分子数Nも変わるため、右の比率と併せて読む。')

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.8))
    curve(axes[0], 'A', 'クラスタ等重み平均 N/K［分子］')
    curve(axes[1], 'Sall', '分子で重み付けした平均 Sall［分子］', log=True)
    curve(axes[2], 'Sfin', '最大クラスタを1個除く Sfin［分子］')
    save(fig, '02_mean_cluster_sizes', '平均の定義と最大クラスタの扱いによる違い', common+'\n中央は縦軸対数。Sfinの除外後に分子が残らない試料は欠測として有効数を記録し、0とはしない。')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    curve(axes[0], 'p', '末端結合率 p = 2E/(4N)［%］', 100)
    curve(axes[1], 'S1', '最大クラスタ占有率 21N1/L³［%］', 100)
    save(fig, '03_bonding_and_largest_volume', '占有率を増やしたときの結合と最大クラスタの成長', common+'\npは入力値ではなく成立結合から測定した量。右は格子全体を分母にした割合で、N1/Nとは異なる。')

    fig, ax = plt.subplots(figsize=(8.5, 6))
    curve(ax, 'S1', '最大クラスタ占有率 S1［%］', 100, xmetric='p')
    save(fig, '04_largest_volume_vs_bonding', '末端結合率と最大クラスタ占有率の対応', '横・縦の誤差棒は各量のSE。各点は1つの占有率条件の平均。両量の共分散は表示していない。\n別々の占有率で生成した系の比較であり、同じ系の反応時間変化や因果効果を表さない。')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    curve(axes[0], 'K', '孤立分子を含むクラスタ数 K［個］', log=True)
    curve(axes[1], 'K_density', '格子点あたりのクラスタ数 K/L³')
    save(fig, '05_cluster_counts', 'クラスタ数と系体積あたりのクラスタ数', common+'\n左は縦軸対数。φを変えるとNも変わるため、個数の増減を同じ分子集合の合体過程とは解釈しない。')

    fig, ax = plt.subplots(figsize=(9, 6))
    curve(ax, 'Sfin', '最大クラスタを1個除く Sfin［分子］', limits=config['zoom_phi'])
    for peak in peaks:
        ax.annotate(f'L={peak["L"]}: 目標φ={peak["target_phi"]:.2f}',
                    (peak['actual_phi'], peak['Sfin_mean']), xytext=(10, 12), textcoords='offset points',
                    fontsize=10, color=colors[peak['L']])
    save(fig, '06_sfin_peak_region', '大きな有限クラスタが発達する占有率領域', common+'\n表示範囲は設定した拡大窓。注記は全占有率で平均最大の測定点であり、厳密な転移点やその信頼区間ではない。')

    fig, axes = plt.subplots(1, len(sizes), figsize=(5*len(sizes), 5.8), squeeze=False)
    stability = []
    for ax, peak in zip(axes[0], peaks):
        samples = groups[peak['L'], peak['target_phi']]
        ns = sorted(set(n for n in config['prefix_counts'] if n <= len(samples)) | {len(samples)})
        stats = [describe(s['Sfin'] for s in samples[:n]) for n in ns]
        ax.errorbar(ns, [v['mean'] for v in stats], yerr=[v['se'] or 0 for v in stats], fmt='o', capsize=3, color=colors[peak['L']])
        ax.set(xlabel='先頭から使用した試料数', ylabel='Sfin 平均 ± SE［分子］',
               title=f'L={peak["L"]}, 目標φ={peak["target_phi"]:.2f}')
        ax.grid(axis='y', alpha=.18)
        for n, v in zip(ns, stats):
            stability.append({'L': peak['L'], 'target_phi': peak['target_phi'], 'prefix_n': n, **v})
    save(fig, '07_sample_count_stability', '試料数を増やしたときの平均とSEの変化', 'case_id順の先頭試料を使用。前回の5試料そのものとの比較ではない。各点の試料集合は入れ子。\n各Lの100試料平均Sfinが最大の条件を例示。全条件の収束や100試料の十分性を保証しない。')
    return peaks, stability, fig_names, matplotlib.__version__


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def guide_text(peaks, summary):
    peak_lines = '\n'.join(f'|{r["L"]}|{r["target_phi"]:.2f}|{r["Sfin_mean"]:.2f} ± {r["Sfin_se"]:.2f}|{r["Sfin_sd"]:.2f}|' for r in peaks)
    return '''# 基本クラスタ量：目的・読み方・報告での説明

## 最初に伝える研究上の問い
占有率を変えると、分子ネットワークのつながり方はどう変化するか。また、その変化は格子サイズによってどう異なるか。
この作図はH2を入力として使わず、生成済みネットワークの基本構造を調べる。異なる占有率の試料は、同じ系を時間発展させたデータではない。

## 01 最大クラスタの絶対サイズと比率
- 目的：最も大きな集団が何分子に成長し、全分子のどれだけを含むかを確認する。
- 定義：N1は最大クラスタ分子数、Pmax=N1/N。Nは孤立分子も含む全分子数。
- 物理的意味：絶対サイズだけだと大きな系が有利なので、比率を併せて比較する。
- 説明例：「占有率の増加に伴って、最大クラスタが全分子のどの程度を占めるようになるかを調べました。系サイズが違うため、絶対分子数と比率を分けています。」
- 限界：有限系でPmaxが大きくなることだけでは、無限系の転移点・貫通・ゲルの弾性を確定できない。

## 02 平均クラスタサイズの3つの定義
- A=N/Kはクラスタを1個選んだときの平均分子数。孤立分子が多いと小さくなる。
- Sall=Σs²/Nは分子を1個選んだとき、それが属するクラスタの平均分子数。最大クラスタも含み、高占有率ではその寄与が大きい。
- Sfin=Σ′s²/Σ′sは最大クラスタ1個を除いた同種の平均。除外後に残る分子から1個選んだときの所属クラスタの平均サイズ。
- 目的：最大クラスタの成長と、それ以外の集団の発達を分けて見る。
- 説明例：「巨大なクラスタができた後は、それを含む平均が支配されます。そこで、最大クラスタを除いた平均も調べ、周囲に残るクラスタの発達を見ました。」
- 限界：Sfinはクラスタの個数ではない。大きな第2・第3クラスタも二乗で強く寄与する。同率最大は1個だけ除く。

## 03 占有率と末端結合率・最大クラスタ占有率
- 目的：分子を多く配置したとき、末端の結合がどれだけ増え、格子全体に対して最大クラスタがどれだけ大きくなるかを見る。
- 定義：p=2E/(4N)、S1=21N1/L³。前者の分母は全末端、後者は全格子点。
- 説明例：「分子濃度を増やした結果として結合率pがどう変化するかと、大きな集団の形成を比較しました。」
- 限界：pは独立に設定した反応確率ではなく測定結果。S1=φPmaxは各試料で成り立つが、別々に平均した量の積が平均S1に一致するとは限らない。

## 04 末端結合率と最大クラスタ占有率
- 目的：どれだけの末端が結合した状態で、大きな集団が形成されるかを記述する。
- 説明例：「横軸を結合済み末端の割合に変え、結合量と大きなクラスタの形成の対応を調べました。」
- 限界：結合を増やした因果効果や時間発展ではない。pとS1は同じ試料から測定され、共分散を持つ。誤差棒は周辺SEのみ。

## 05 クラスタ数とクラスタ数密度
- 目的：集団が何個に分かれているかと、体積をそろえたときの数を確認する。
- 定義：Kは孤立分子を含むクラスタ数、K/L³は格子点あたりの数。
- 説明例：「クラスタの大きさに加え、全体が何個の集団に分かれるかを調べました。Lによる体積の違いを除くため、単位格子体積あたりの数も示します。」
- 限界：φを変えるとNも変わる。Kの増減だけで、同じ集団の合体・分裂を直接観測したとは言えない。

## 06 Sfinのピーク領域
- 目的：最大クラスタ以外にも大きな集団が発達する領域を、複数Lで比較する。
- 説明例：「Sfinの平均が最大となる占有率と、その系サイズ依存を確認しました。測定点のピークを転移候補の手がかりとして扱います。」
- 限界：ピークの位置は有限サイズ、占有率刻み、統計誤差に依存する。注記は平均の最大点で、転移点や位置の信頼区間ではない。

## 07 試料数による平均・SEの変化
- 目的：少ない試料で見えた平均が、試料追加でどの程度変わるかを見る。
- 説明例：「各Lのピーク条件を例に、5試料から100試料まで増やしたときの平均とSEを比較しました。」
- 限界：先頭試料の部分集合は互いに重なる。前回の5試料そのものとの直接比較ではなく、今回の100試料の部分集合を使う。ケースの順序にも依存する。SEが小さくても偏りは消えず、この図だけで十分な試料数とは宣言しない。

## 平均、SD、SEの区別
SDは試料間のばらつき。SE=SD/√nは試料平均の不確かさ。有効数n=100ならSE=SD/10だが、試料間のばらつき自体が10分の1になったわけではない。
すべての比率は各試料で計算してから平均する。試料のクラスタを全て合算して比を作る集計とは異なる。

## 今回のピークの観測
|L|平均Sfinが最大の目標占有率|平均 ± SE [分子]|試料SD [分子]|
|---:|---:|---:|---:|
''' + peak_lines + '''

この表の最大点は全占有率から選択した記述的な結果。隣接点との差の有意性や、無限系の臨界点を示す表ではない。
生成・結合規則を変えず、保存済みsamples.csvから再集計した。元の結合・クラスタCSVからネットワーク構造を独立に再検証したという意味ではない。
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=REPO/'configs/analysis/cluster-basics.json')
    parser.add_argument('--results', type=Path, help='完了済みresults_*、実験フォルダ、またはsamples.csv')
    parser.add_argument('--output', type=Path, help='新しい出力フォルダ（既存の内容は上書きしない）')
    parser.add_argument('--font', type=Path, help='日本語フォントファイル')
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding='utf-8'))
    source = args.results or REPO/config['results']
    input_path = resolve_samples(source).resolve()
    groups = load_groups(input_path, config)
    summary = summarize(groups)
    # Confirm the environment before creating the output directory.
    setup_font(args.font)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output = (args.output or REPO/'data/figures'/f'cluster_basics_{stamp}').resolve()
    output.mkdir(parents=True, exist_ok=False)
    peaks, stability, names, matplotlib_version = draw(summary, groups, output, config, args.font)
    write_csv(output/'basic_summary.csv', summary)
    write_csv(output/'peak_summary.csv', peaks)
    write_csv(output/'sample_count_stability.csv', stability)
    (output/'figure_guide.md').write_text(guide_text(peaks, summary), encoding='utf-8')
    metadata = {'created_at_utc': datetime.now(timezone.utc).isoformat(),
                'input': str(input_path), 'input_sha256': hashlib.sha256(input_path.read_bytes()).hexdigest(),
                'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'config': config, 'config_path': str(args.config.resolve()),
                'python': platform.python_version(), 'matplotlib': matplotlib_version,
                'conditions': len(groups), 'samples': sum(map(len, groups.values())),
                'figure_names': names, 'aggregation': 'sample-level ratios, then equal-weight mean',
                'missing_policy': 'exclude undefined values; report metric-specific n',
                'se_definition': 'sample SD with ddof=1 / sqrt(valid n)',
                'note': 'No simulation or H2 calculation; no independent reconstruction from raw bonds.'}
    (output/'run_metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    archive = output/'cluster_basics_share.zip'
    temporary_archive = output/'.cluster_basics_share.tmp.zip'
    with ZipFile(temporary_archive, 'w', ZIP_DEFLATED) as zip_file:
        for item in sorted(output.iterdir()):
            if item.is_file() and item not in (archive, temporary_archive):
                zip_file.write(item, item.name)
    with ZipFile(temporary_archive) as zip_file:
        if zip_file.testzip() is not None:
            raise ValueError('共有ZIPの検証に失敗しました')
    temporary_archive.replace(archive)
    print(f'入力: {input_path}')
    print(f'{metadata["conditions"]}条件・{metadata["samples"]:,}試料を再集計しました')
    print(f'図と説明: {output}')
    print(f'共有用ZIP: {archive}')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(2)
