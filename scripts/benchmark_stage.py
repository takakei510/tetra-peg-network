#!/usr/bin/env python3
"""ベンチマークの後処理・描画を独立プロセスで行う内部入口。"""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts' / 'visualization'))


def plot(run, dpi):
    from analyze_graph import draw_graph, read_rows
    from plot_molecules import read_molecules, read_bonds, draw_molecules
    started = time.perf_counter()
    molecules, occupied = read_molecules(run / 'trajectories.csv')
    points = read_bonds(run / 'bonds.csv', molecules)
    edges = read_rows(run / 'bonds.csv', ('molecule_a', 'molecule_b'))
    info = json.loads((run / 'graph' / 'graph-analysis.json').read_text())
    members = read_rows(run / 'graph' / 'largest-members.csv', ('molecule',))
    chosen = {int(row['molecule']) for row in members}
    selected = {m: item for m, item in molecules.items() if m in chosen}
    selected_edges = [row for row in edges if int(row['molecule_a']) in chosen]
    selected_points = [p for row, p in zip(edges, points) if int(row['molecule_a']) in chosen]
    read_finished = time.perf_counter()
    output = run / 'plots'
    output.mkdir(exist_ok=False)
    draw_molecules(molecules, occupied, output / 'molecules.png', dpi, bonds=points)
    draw_graph(molecules, edges, output / 'network.png', 'Molecular multigraph', dpi)
    draw_graph(selected, selected_edges, output / 'largest-network.png',
               f"Largest component {info['selected_cluster_id']}", dpi)
    selected_occupied = sum(1 + sum(len(path) - 1 for path in item[2].values())
                            for item in selected.values())
    draw_molecules(selected, selected_occupied, output / 'largest-molecules.png',
                   dpi, bonds=selected_points)
    (run / 'plot-timing.json').write_text(json.dumps({
        'read_validate_seconds': read_finished - started,
        'render_save_seconds': time.perf_counter() - read_finished,
        'note': 'render includes lazy Matplotlib import and PNG writes; peak RSS covers whole process'
    }, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['postprocess', 'plot'])
    parser.add_argument('run', type=Path)
    parser.add_argument('--dpi', type=int, default=120)
    args = parser.parse_args()
    if args.stage == 'postprocess':
        from analyze_graph import analyze
        analyze(args.run, args.run / 'graph', no_plot=True)
    else:
        plot(args.run, args.dpi)


if __name__ == '__main__':
    main()
