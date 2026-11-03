"""All-concentration H2 maps and detail at the measured Sfin maximum."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--summary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font', type=Path, required=True)
    p.add_argument('--L', type=int, default=128)
    p.add_argument('--min-samples', type=int, default=50)
    args = p.parse_args()
    font_manager.fontManager.addfont(str(args.font))
    name = font_manager.FontProperties(fname=args.font).get_name()
    plt.rcParams.update({'font.family': name, 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'path', 'axes.unicode_minus': False})
    with args.input.open(encoding='utf-8-sig', newline='') as f:
        rows = [r for r in csv.DictReader(f) if int(r['L']) == args.L]
    with args.summary.open(encoding='utf-8-sig', newline='') as f:
        summary = [r for r in csv.DictReader(f) if int(r['L']) == args.L]
    if not rows or not summary:
        raise ValueError(f'L={args.L}: データなし')
    peak = max(summary, key=lambda r: float(r['cluster_weighted_excl_one_mean']))
    peak_phi = round(float(peak['target_phi']), 2)
    phis = sorted({round(float(r['target_phi']), 2) for r in rows})
    if set(phis) != {round(float(r['target_phi']), 2) for r in summary}:
        raise ValueError('集計CSV間で占有率の集合が一致しません')
    sample_counts = {int(r['total_samples']) for r in rows}
    if len(sample_counts) != 1:
        raise ValueError('条件間で総試料数が異なります')
    total_samples = sample_counts.pop()
    bins = sorted({(int(r['size_low']), int(r['size_high'])) for r in rows if int(r['size_low']) > 1})
    data = {(round(float(r['target_phi']), 2), r['scope'], int(r['size_low'])): r for r in rows}
    if len(data) != len(rows):
        raise ValueError('重複する集計行があります')
    labels = [str(lo) if lo == hi else f'{lo:,}〜{hi:,}' for lo, hi in bins]
    n = np.zeros((len(bins), len(phis)), dtype=int)
    q = np.full(n.shape, np.nan)
    for i, (lo, _) in enumerate(bins):
        for j, phi in enumerate(phis):
            r = data.get((phi, 'exclude_one_largest', lo))
            if r:
                n[i, j] = int(r['fraction_n'])
                if n[i, j] > 0:
                    q[i, j] = 100 * float(r['fraction_mean'])
    args.output.mkdir(parents=True, exist_ok=True)

    def save(fig, stem):
        for ext in ['png', 'svg']:
            fig.savefig(args.output / f'{stem}.{ext}', dpi=180, bbox_inches='tight', pad_inches=.2)
        plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(13.2, 9.8), sharex=True)
    fig.subplots_adjust(left=.15, right=.88, top=.87, bottom=.17, hspace=.25)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#e4e7eb')
    grid_x = np.arange(len(phis)+1)-.5
    grid_y = np.arange(len(bins)+1)-.5
    masked = np.ma.masked_where((n < args.min_samples) | ~np.isfinite(q), q)
    im = axes[0].pcolormesh(grid_x, grid_y, masked, cmap=cmap, shading='flat', vmin=0)
    imn = axes[1].pcolormesh(grid_x, grid_y, n, cmap='Blues', shading='flat', vmin=0, vmax=total_samples)
    for ax in axes:
        ax.set_yticks(range(len(bins)), labels)
        ax.set_ylabel('クラスタサイズ区間［分子］')
        ax.set_ylim(-.5, len(bins)-.5)
        ax.axvline(phis.index(peak_phi), color='#e27525', linestyle='--', linewidth=1.7)
    ticks = [i for i, phi in enumerate(phis) if round(phi*100) % 10 == 0]
    axes[1].set_xticks(ticks, [f'{phis[i]:.2f}' for i in ticks])
    axes[1].set_xlabel(r'目標占有率 $\phi_{\rm target}$')
    axes[0].set_title(f'$H_2$末端割合：寄与試料数 n ≥ {args.min_samples} の区間を表示', loc='left')
    axes[1].set_title('寄与試料数：そのサイズ区間が存在した試料数', loc='left')
    fig.colorbar(im, ax=axes[0], fraction=.035, pad=.025, label='平均$H_2$末端割合［%］')
    fig.colorbar(imn, ax=axes[1], fraction=.035, pad=.025, label='寄与試料数 n')
    fig.suptitle(f'全占有率での$H_2$末端割合と比較可能なサイズ範囲\nL={args.L} ｜ 各条件{total_samples}試料 ｜ 各試料の最大クラスタを1個除外', fontsize=16, y=.97)
    fig.text(.15, .10, f'橙の破線：平均Sfinが最大の測定点 φ={peak_phi:.2f}。厳密な転移点ではない。', fontsize=10)
    fig.text(.15, .067, f'上段の灰色：寄与試料数{args.min_samples}未満（未観測を含む）。50試料は表示上の目安で、精度の保証ではない。'.replace('50試料', f'{args.min_samples}試料'), fontsize=9)
    fig.text(.15, .034, '$H_2$割合は試料内でクラスタを等重み平均し、次に寄与試料を等重み平均。$H_2$=0も含む。区間は等間隔配置。', fontsize=9)
    save(fig, 'all_phi_h2_and_support')

    # Aligned detailed panels: total cluster counts above, finite-cluster H2 below.
    detail_bins = [(1, 1)] + bins
    relevant = [i for i, (lo, hi) in enumerate(detail_bins) if (peak_phi, 'all', lo) in data]
    detail_bins = detail_bins[:max(relevant)+1]
    detail_labels = [str(lo) if lo == hi else f'{lo:,}〜{hi:,}' for lo, hi in detail_bins]
    fig, (ac, aq, an) = plt.subplots(3, 1, figsize=(13.2, 8.5), sharex=True,
                                     gridspec_kw={'height_ratios': [1, 1, .14]})
    fig.subplots_adjust(left=.12, right=.97, top=.84, bottom=.25, hspace=.22)
    color = '#2873b8'
    for i, (lo, hi) in enumerate(detail_bins):
        r = data.get((peak_phi, 'all', lo))
        count = int(r['clusters'])/int(r['total_samples']) if r else 0
        if count > 0:
            ac.scatter(i, count, color=color, s=28)
        else:
            ac.text(i, .025, '0', transform=ac.get_xaxis_transform(), ha='center', fontsize=9)
        r = data.get((peak_phi, 'exclude_one_largest', lo))
        if lo == 1:
            ns_label = '省略'
        elif r and int(r['fraction_n']) > 0:
            ns = int(r['fraction_n'])
            y = 100*float(r['fraction_mean'])
            ns_label = str(ns)
            if ns >= args.min_samples:
                aq.errorbar(i, y, yerr=100*float(r['fraction_se']), fmt='o', color=color, capsize=3)
            else:
                aq.scatter(i, y, facecolors='none', edgecolors='#888888', s=30)
        else:
            ns_label = '—'
        an.text(i, .5, ns_label, ha='center', va='center', fontsize=9)
    ac.set_yscale('log')
    ac.set_ylabel('平均クラスタ数\n［個／試料］')
    ac.set_title('全クラスタを含む個数（区間幅では割らない）', loc='left')
    aq.set_ylabel('平均$H_2$末端割合［%］')
    aq.set_ylim(0, max(11, aq.get_ylim()[1]))
    aq.set_title('各試料の最大クラスタを1個除外した$H_2$末端割合', loc='left')
    for ax in [ac, aq]:
        ax.grid(axis='y', alpha=.2)
        ax.tick_params(axis='x', labelbottom=False)
    an.set_ylim(0, 1)
    an.set_yticks([])
    an.set_ylabel('寄与試料数', rotation=0, ha='right', va='center', labelpad=12)
    for spine in an.spines.values():
        spine.set_visible(False)
    an.set_xticks(range(len(detail_bins)), detail_labels, rotation=40, ha='right')
    an.set_xlim(-.6, len(detail_bins)-.4)
    an.set_xlabel('クラスタサイズ区間［分子］')
    fig.suptitle(f'平均Sfinが最大の条件でのクラスタサイズと$H_2$末端割合\nL={args.L} ｜ φtarget={peak_phi:.2f} ｜ 各条件{total_samples}試料', fontsize=16, y=.97)
    fig.text(.12, .075, f'色付き点：平均 ± SE（n≥{args.min_samples}）。白抜き：n<{args.min_samples}。—：該当なし。1分子の$H_2$割合は常に0のため省略。', fontsize=9)
    fig.text(.12, .04, 'クラスタ数は全試料平均。$H_2$割合は区間のある試料のみの平均。横軸はカテゴリ軸。相関から因果は判断しない。', fontsize=9)
    save(fig, 'peak_phi_detail')
    (args.output / 'analysis_metadata.json').write_text(json.dumps({
        'L': args.L, 'peak_phi': peak_phi, 'min_samples_display': args.min_samples,
        'peak_rule': 'maximum mean Sfin among measured concentrations',
        'total_samples': total_samples, 'phis': phis,
        'notes': ['Descriptive selection, not a critical-point estimate.',
                  'H2 means are conditional on the size interval being observed.']}, indent=2)+'\n', encoding='utf-8')
    print(f'L={args.L}: measured Sfin maximum at phi={peak_phi:.2f}')
    print('Output:', args.output)


if __name__ == '__main__':
    main()
