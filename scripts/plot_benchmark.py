#!/usr/bin/env python3
"""性能測定の集計CSVから、実行環境比較用の図を作る。"""
import argparse
import csv
from pathlib import Path
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('suite', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('出力先が既に存在します。別名を指定してください')
    with (args.suite / 'summary.csv').open(newline='') as stream:
        summary = list(csv.DictReader(stream))
    with (args.suite / 'cases.csv').open(newline='') as stream:
        cases = list(csv.DictReader(stream))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    phis = sorted({row['target_phi'] for row in summary}, key=float)
    colors = plt.get_cmap('tab10')
    for index, phi in enumerate(phis):
        color = colors(index)
        for stage, style in [('c', '-'), ('postprocess', '--')]:
            rows = sorted((r for r in summary if r['target_phi'] == phi and r['stage'] == stage
                           and r['wall_seconds_median']), key=lambda r: int(r['L']))
            if rows:
                x = [int(r['L']) for r in rows]
                y = [float(r['wall_seconds_median']) for r in rows]
                axes[0, 0].plot(x, y, style+'o', color=color, label=f'phi={phi}, {stage}')
                axes[0, 1].plot(x, [float(r['peak_rss_bytes_median'])/2**20 for r in rows],
                                style+'o', color=color, label=f'phi={phi}, {stage}')
        rows = sorted((r for r in summary if r['target_phi'] == phi and r['stage'] == 'plot'
                       and r['render_save_seconds_median']), key=lambda r: int(r['L']))
        if rows:
            axes[1, 0].plot([int(r['L']) for r in rows],
                            [float(r['render_save_seconds_median']) for r in rows],
                            '-o', color=color, label=f'phi={phi}')
        groups = {}
        for row in cases:
            if row['target_phi'] == phi and row['generation_status'] == 'completed':
                groups.setdefault(int(row['L']), []).append(row)
        for metric, style, name in [('raw_csv_bytes', '-', 'raw'), ('gzip_csv_bytes', '--', 'gzip')]:
            if groups:
                x = sorted(groups)
                y = [statistics.median(float(r[metric]) for r in groups[L])/2**20 for L in x]
                axes[1, 1].plot(x, y, style+'o', color=color, label=f'phi={phi}, {name}')
    labels = [('Process wall time (startup + I/O)', 'Seconds'),
              ('Process peak RSS (do not add stages)', 'MiB'),
              ('4 PNGs: render + save incl. Matplotlib import', 'Seconds'),
              ('Original 7 CSV files per sample', 'MiB')]
    for axis, (title, ylabel) in zip(axes.flat, labels):
        axis.set_title(title, fontsize=10)
        axis.set_xlabel('L')
        axis.set_ylabel(ylabel)
        axis.set_xscale('log', base=2)
        axis.set_yscale('log')
        axis.grid(True, which='both', alpha=.2)
        axis.legend(fontsize=7)
    fig.suptitle('AB-SAW resource benchmark — medians across seeds / repeats', fontsize=13)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    print(args.output)


if __name__ == '__main__':
    main()
