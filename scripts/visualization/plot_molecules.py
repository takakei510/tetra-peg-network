#!/usr/bin/env python3
"""複数の4腕分子を同じ絶対座標系に描画する（現在は開境界専用）。

例: python scripts/visualization/plot_molecules.py data/two-molecules.csv \
        --output data/two-molecules.png

1分子用plot_molecule.pyのCSV読込と整合性検査を直接再利用する。
分子ごとに原点を移すと位置関係が失われるので、ここでは移動しない。
"""
import argparse
import csv
import sys
from pathlib import Path

from plot_molecule import read_molecule


def read_molecules(csv_path):
    """CSV内の分子IDを列挙し、個々の軌跡と分子間重複を確認する。"""
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if "molecule" not in (reader.fieldnames or []):
            raise ValueError("CSVにmolecule列がありません")
        ids = sorted({int(row["molecule"]) for row in reader})
    if not ids or ids[0] < 0:
        raise ValueError("CSVに有効な分子がありません")
    molecules = {}
    occupied = set()
    for molecule_id in ids:
        center, kind, paths = read_molecule(csv_path, molecule_id)
        # 中心は4腕の共通起点なので1回だけ数える。
        sites = {center}
        for path in paths.values():
            sites.update(path[1:])
        if occupied.intersection(sites):
            raise ValueError(f"分子{molecule_id}と別分子の格子点が重複しています")
        occupied.update(sites)
        molecules[molecule_id] = (center, kind, paths)
    return molecules, len(occupied)


def read_candidates(csv_path, molecules):
    """候補CSVを軌跡と照合する。点線は潜在ペアであり結合ではない。"""
    pairs, seen = [], set()
    with csv_path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {"endpoint_a", "endpoint_b", "molecule_a", "molecule_b",
                    "arm_a", "arm_b", "type_a", "type_b", "manhattan_distance"}
        required.update(f"{axis}_{side}" for axis in "xyz" for side in "ab")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("候補CSVの列が不足しています")
        for row in reader:
            ids = tuple(int(row[f"endpoint_{side}"]) for side in "ab")
            if ids[0] >= ids[1] or ids in seen:
                raise ValueError("候補IDの順序または重複が不正です")
            points, kinds = [], []
            for side, endpoint in zip("ab", ids):
                m, arm = int(row[f"molecule_{side}"]), int(row[f"arm_{side}"])
                if endpoint != 4*m + arm or not 0 <= arm < 4:
                    raise ValueError("末端IDと分子・腕IDが一致しません")
                point = tuple(int(row[f"{axis}_{side}"]) for axis in "xyz")
                _, kind, paths = molecules[m]
                if point != paths[arm][-1] or kind != row[f"type_{side}"]:
                    raise ValueError("候補の末端位置・型と軌跡が一致しません")
                points.append(point)
                kinds.append(kind)
            distance = sum(abs(a-b) for a, b in zip(*points))
            if kinds[0] == kinds[1] or not 1 <= distance <= 2 or distance != int(row["manhattan_distance"]):
                raise ValueError("A-B条件または距離条件に違反しています")
            seen.add(ids)
            pairs.append(tuple(points))
    return pairs


def draw_molecules(molecules, occupied_count, output, dpi, show=False, candidates=None):
    """A/Bを色、中心と末端を形で区別し、絶対座標で保存する。"""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    colors = {"A": "#4472C4", "B": "#C0504D"}
    all_points = []
    for molecule_id, (center, kind, paths) in molecules.items():
        for arm, path in paths.items():
            all_points.extend(path)
            xs, ys, zs = zip(*path)
            # 凡例は分子ごとに1項目。4本の腕は同じ型の色で描く。
            ax.plot(xs, ys, zs, "-o", color=colors[kind], markersize=3.5,
                    linewidth=2.2, label=f"Molecule {molecule_id} ({kind})" if arm == 0 else None)
            ax.scatter(xs[-1], ys[-1], zs[-1], color=colors[kind], marker="s",
                       s=55, depthshade=False)
        ax.scatter(*center, color="black", s=95, depthshade=False)
        ax.text(*center, f"  {kind}{molecule_id}", color="black", fontsize=11)
    # 候補は点線で描き、確定した結合と誤認しない凡例を付ける。
    if candidates is not None:
        for i, (start, end) in enumerate(candidates):
            xs, ys, zs = zip(start, end)
            ax.plot(xs, ys, zs, "--", color="#23854B", linewidth=2.5,
                    label="Candidate (not bond)" if i == 0 else None)
    # 黒丸と四角の説明も凡例に加える。
    ax.scatter([], [], [], color="black", s=70, label="Center")
    ax.scatter([], [], [], color="#777777", marker="s", s=55, label="Endpoint")
    ranges = [(min(p[i] for p in all_points), max(p[i] for p in all_points))
              for i in range(3)]
    span = max(high - low for low, high in ranges) + 1
    for setter, (low, high) in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), ranges):
        middle = (low + high) / 2
        setter(middle - span / 2, middle + span / 2)
    ax.set_box_aspect((1, 1, 1))
    for axis, setter in zip("xyz", (ax.set_xlabel, ax.set_ylabel, ax.set_zlabel)):
        setter(axis, labelpad=2)
    detail = "Absolute coordinates, no bonds" if candidates is None else f"{len(candidates)} candidate pairs, no bonds"
    ax.set_title(f"{len(molecules)} molecules, {occupied_count} occupied sites\n" + detail)
    ax.view_init(elev=24, azim=-55)
    ax.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output}")
    print(f"Molecules: {len(molecules)}, occupied sites: {occupied_count}, no overlap")
    if show:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, default=Path("data/two-molecules.png"))
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--candidates", type=Path, help="③-Aの候補CSV。確定結合ではない")
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("dpiは正にしてください")
    try:
        molecules, occupied_count = read_molecules(args.csv)
        candidates = read_candidates(args.candidates, molecules) if args.candidates else None
        draw_molecules(molecules, occupied_count, args.output, args.dpi, args.show, candidates)
    except (ValueError, OSError, ImportError, TypeError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
