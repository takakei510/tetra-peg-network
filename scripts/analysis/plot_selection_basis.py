"""Plot concentration-selection evidence from existing 100-sample summaries.

python plot_selection_basis.py --summary summary.csv --output figures --font FONT
Selection uses Pmax and S_fin only, independently of the H2-size curves.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, default=Path('summary.csv'))
    parser.add_argument('--output', type=Path, default=Path('figures'))
    parser.add_argument('--font', type=Path)
    parser.add_argument('--L', type=int, default=128)
    parser.add_argument('--low-threshold', type=float, default=0.05)
    parser.add_argument('--high-threshold', type=float, default=0.80)
    args = parser.parse_args()
    if args.font:
        font_manager.fontManager.addfont(str(args.font))
        family = font_manager.FontProperties(fname=args.font).get_name()
    else:
        family = 'Noto Sans CJK JP'
    plt.rcParams.update({'font.family': family, 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.unicode_minus': False, 'svg.fonttype': 'path'})
    with args.summary.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    groups = {}
    for row in rows:
        groups.setdefault(int(row['L']), []).append(row)
    for group in groups.values():
        group.sort(key=lambda row: float(row['target_phi']))
    data = groups[args.L]
    assert len({float(r['target_phi']) for r in data}) == len(data)
    peak = max(data, key=lambda row: float(row['cluster_weighted_excl_one_mean']))
    peak_phi = float(peak['target_phi'])
    low = max((r for r in data if float(r['target_phi']) < peak_phi and
               float(r['Pmax_mean']) <= args.low_threshold),
              key=lambda row: float(row['target_phi']))
    high = min((r for r in data if float(r['target_phi']) > peak_phi and
                float(r['Pmax_mean']) >= args.high_threshold),
               key=lambda row: float(row['target_phi']))
    selected = [low, peak, high]
    colors = ['#2873b8', '#bd681b', '#218477']
    labels = ['最大クラスタが小さい', '平均クラスタサイズが最大', '最大クラスタに8割が所属']
    args.output.mkdir(parents=True, exist_ok=True)

    def curve(ax, group, metric, color, scale=1, label=None):
        x = [float(r['target_phi']) for r in group]
        y = [scale * float(r[metric + '_mean']) for r in group]
        se = [scale * float(r[metric + '_se']) for r in group]
        ax.errorbar(x, y, yerr=se, color=color, lw=1.6, marker='o', ms=3,
                    elinewidth=.8, capsize=2, label=label)
        ax.grid(axis='y', alpha=.20)
        ax.set_xlabel(r'目標占有率 $\phi_{\rm target}$')
        ax.set_xlim(.09, .71)
        ax.xaxis.set_major_locator(MultipleLocator(.1))
        ax.xaxis.set_minor_locator(MultipleLocator(.01))

    fig, axes = plt.subplots(1, 2, figsize=(14.4, 8.6))
    fig.subplots_adjust(left=.08, right=.97, top=.79, bottom=.37, wspace=.30)
    fig.text(.06, .94, '代表占有率の選定根拠', fontsize=23, weight='bold')
    fig.text(.06, .886,
             f'L = {args.L}  ｜  占有率0.10〜0.70（0.01刻み・61条件）  ｜  各100試料の平均 ± SE',
             fontsize=12, color='#445567')
    curve(axes[0], data, 'Pmax', '#536b80', 100)
    curve(axes[1], data, 'cluster_weighted_excl_one', '#536b80')
    axes[0].set_title('① 最大クラスタが全分子に占める割合', loc='left', pad=14, fontsize=13)
    axes[1].set_title('② 最大クラスタを除いた平均クラスタサイズ', loc='left', pad=14, fontsize=13)
    axes[0].set_ylabel(r'最大クラスタ比 $P_{\max}=N_1/N$ ［%］')
    axes[0].set_ylim(-2, 102)
    axes[1].set_ylabel(r'$S_{\rm fin}$ ［分子］')
    axes[1].set_ylim(bottom=0)
    for threshold in [args.low_threshold, args.high_threshold]:
        axes[0].axhline(100*threshold, color='#aab5bf', ls='--', lw=1)
        axes[0].text(.105, 100*threshold+2, f'基準 {100*threshold:g}%', fontsize=10, color='#657586')
    for row, color in zip(selected, colors):
        phi = float(row['target_phi'])
        for ax, metric, scale in [(axes[0], 'Pmax', 100),
                                   (axes[1], 'cluster_weighted_excl_one', 1)]:
            ax.axvline(phi, color=color, lw=1.1, ls='--', alpha=.75)
            ax.scatter([phi], [scale*float(row[metric+'_mean'])], color=color,
                       edgecolor='white', linewidth=1, s=90, zorder=5)
            ax.text(phi, 1.025, f'{phi:.2f}', transform=ax.get_xaxis_transform(),
                    ha='center', fontsize=10, color=color)
    # The first two concentrations are close; stagger their labels for legibility.
    for ax in axes:
        for t in ax.texts:
            if t.get_text() == f"{float(low['target_phi']):.2f}":
                t.set_position((float(low['target_phi'])-.013, 1.025))
                t.set_ha('right')
            elif t.get_text() == f'{peak_phi:.2f}':
                t.set_position((peak_phi+.007, 1.025))
                t.set_ha('left')
    rules = [f'ピークより低い側で、平均Pmax ≤ {100*args.low_threshold:g}%となる最大の占有率',
             '平均Sfinが最大となる測定点',
             f'ピークより高い側で、平均Pmax ≥ {100*args.high_threshold:g}%となる最小の占有率']
    for i, (row, color, label, rule) in enumerate(zip(selected, colors, labels, rules)):
        x = .07 + i*.31
        fig.text(x, .265, f'φ = {float(row["target_phi"]):.2f}  |  {label}', color=color, fontsize=11, weight='bold')
        # Wrap criterion to fit its column.
        wrapped = rule.replace('平均Pmax', '\n平均Pmax').replace('となる最大', '\nとなる最大').replace('となる最小', '\nとなる最小')
        fig.text(x, .223, wrapped, fontsize=10, va='top', linespacing=1.55)
        fig.text(x, .105, f'Pmax = {100*float(row["Pmax_mean"]):.2f}%   /   Sfin = {float(row["cluster_weighted_excl_one_mean"]):.2f}', fontsize=10)
    fig.text(.07, .055, r'$S_{\rm fin}=\sum_i s_i^2 / \sum_i s_i$（各試料の最大クラスタを1個除外して算出し、その後に試料間平均）', fontsize=10)
    fig.text(.07, .025, '5%・80%は状態比較のための操作的基準。ピークは測定点中の最大であり、転移点の確定ではない。', fontsize=10, color='#586776')
    for ext in ['png', 'svg']:
        fig.savefig(args.output / f'selection_basis_L{args.L}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.3))
    fig.subplots_adjust(left=.08, right=.97, top=.76, bottom=.19, wspace=.28)
    fig.text(.06, .92, '格子サイズによる成長曲線とピーク位置の違い', fontsize=20, weight='bold')
    fig.text(.06, .855, '全183条件・18,300試料  ｜  各点：100試料平均 ± SE  ｜  同じ占有率での比較', fontsize=11, color='#445567')
    peaks = []
    for L, color in zip(sorted(groups), colors):
        group = groups[L]
        curve(axes[0], group, 'Pmax', color, 100, f'L = {L}')
        curve(axes[1], group, 'cluster_weighted_excl_one', color, 1, f'L = {L}')
        row = max(group, key=lambda r: float(r['cluster_weighted_excl_one_mean']))
        x, y = float(row['target_phi']), float(row['cluster_weighted_excl_one_mean'])
        axes[1].scatter([x], [y], color=color, s=65, edgecolor='white', zorder=5)
        axes[1].annotate(f'L={L}: φ={x:.2f}', (x, y), xytext=(12, 12),
                         textcoords='offset points', fontsize=10, color=color)
        peaks.append({'L': L, 'phi': x, 'Sfin_mean': y,
                      'Sfin_se': float(row['cluster_weighted_excl_one_se'])})
    axes[0].set_ylabel(r'最大クラスタ比 $P_{\max}$ ［%］')
    axes[0].set_ylim(0, 100)
    axes[1].set_ylabel(r'最大クラスタを除いた平均サイズ $S_{\rm fin}$ ［分子］')
    axes[1].set_ylim(0, 235)
    axes[0].legend(frameon=False, loc='upper left')
    fig.text(.08, .075, 'ピーク位置はLによって異なる。L=128の代表条件を、全サイズ共通の転移点としては扱わない。', fontsize=11)
    fig.text(.08, .035, '横軸は目標占有率。SEは試料間平均の標準誤差であり、ピーク位置の信頼区間ではない。', fontsize=10, color='#586776')
    for ext in ['png', 'svg']:
        fig.savefig(args.output / f'finite_size_overview.{ext}', dpi=180, facecolor='white')
    plt.close(fig)

    result = {'summary_sha256': hashlib.sha256(args.summary.read_bytes()).hexdigest(),
              'L': args.L, 'low_threshold': args.low_threshold,
              'high_threshold': args.high_threshold, 'selection_rules': rules,
              'selected': selected, 'grid_peaks': peaks,
              'note': 'Descriptive selection from existing data; not a preregistered test or a critical-point estimate.'}
    (args.output / 'selection.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print('Selected concentrations:', [float(r['target_phi']) for r in selected])
    print('Saved:', args.output.resolve())


if __name__ == '__main__':
    main()
