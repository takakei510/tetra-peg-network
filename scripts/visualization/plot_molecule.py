#!/usr/bin/env python3
"""demo_moleculeのCSVから、1分子の4腕軌跡を3D描画する。

使用例（リポジトリ直下で実行）:
    python scripts/visualization/plot_molecule.py data/molecule.csv \
        --output data/molecule.png

CSVの読み込みには標準ライブラリ、描画にはMatplotlibだけを使う。
現在の開境界で生成した座標を対象とする。周期境界の座標を連続的な
軌跡に戻す処理は含めない。モデルを変更した場合はここも見直すこと。
"""

import argparse
import csv
import sys
from pathlib import Path


def read_molecule(csv_path, molecule_id):
    """指定分子を読み、中心・分子タイプ・腕ごとの座標列を返す。"""
    selected = []
    required = {"molecule", "type", "arm", "step", "x", "y", "z"}
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CSVに必要な列がありません: " + ", ".join(sorted(required)))
        for line, row in enumerate(reader, start=2):
            try:
                if int(row["molecule"]) != molecule_id:
                    continue
                selected.append((
                    row["type"], int(row["arm"]), int(row["step"]),
                    tuple(int(row[axis]) for axis in ("x", "y", "z")),
                ))
            except (ValueError, TypeError, KeyError) as error:
                raise ValueError(f"CSVの{line}行目を読み込めません") from error

    if not selected:
        raise ValueError(f"分子ID {molecule_id} はCSVにありません")
    types = {row[0] for row in selected}
    if len(types) != 1 or not types.issubset({"A", "B"}):
        raise ValueError("1分子のtypeはAまたはBで統一してください")

    # C側の出力規則では、中心だけがarm=-1, step=0。
    # 各腕のCSV行は中心を含まず、step=1から始まる。
    centers = [row for row in selected if row[1] == -1]
    if len(centers) != 1 or centers[0][2] != 0:
        raise ValueError("中心行（arm=-1, step=0）が1行必要です")
    center = centers[0][3]
    arms = {arm: [] for arm in range(4)}
    for _, arm, step, xyz in selected:
        if arm == -1:
            continue
        if arm not in arms:
            raise ValueError(f"4腕モデルで使用できないarm ID: {arm}")
        arms[arm].append((step, xyz))

    # 行の並びに依存せず、step順に結ぶ。欠損や重複があるCSVを
    # そのまま描くと誤った線になるため、最低限の整合性を検査する。
    paths = {}
    occupied = {center}
    lengths = set()
    for arm, points in arms.items():
        points.sort(key=lambda item: item[0])
        steps = [step for step, _ in points]
        if not steps or steps != list(range(1, len(points) + 1)):
            raise ValueError(f"arm {arm}: stepは1から連続する必要があります")
        lengths.add(len(points))
        previous = center
        path = [center]  # 中心を先頭に補い、最初の1歩も線で表示する。
        for _, xyz in points:
            # 開境界の立方格子では、最近接移動のマンハッタン距離は1。
            if sum(abs(a - b) for a, b in zip(previous, xyz)) != 1:
                raise ValueError(f"arm {arm}: 最近接でない移動があります")
            if xyz in occupied:
                raise ValueError("中心以外の軌跡に格子点の重複があります")
            occupied.add(xyz)
            path.append(xyz)
            previous = xyz
        paths[arm] = path
    if len(lengths) != 1:
        raise ValueError("4本の腕のステップ数が一致していません")
    return center, next(iter(types)), paths


def draw_molecule(center, molecule_type, paths, molecule_id, output, dpi,
                  absolute=False, show=False):
    """軌跡をPNG等に保存し、必要なら画面にも表示する。"""
    import matplotlib

    # 画像保存だけの場合は、画面のないサーバーでも動くAggを選ぶ。
    # --showの場合はユーザー環境のGUIバックエンドを利用する。
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    colors = ("#4472C4", "#C0504D", "#70AD47", "#8064A2")
    all_points = []
    for arm, path in paths.items():
        # 相対座標なら分子中心を(0,0,0)にする。平行移動なので
        # 分子の形・距離・自己回避の性質は変わらない。
        points = [p if absolute else tuple(a - b for a, b in zip(p, center))
                  for p in path]
        all_points.extend(points)
        xs, ys, zs = zip(*points)
        ax.plot(xs, ys, zs, "-o", color=colors[arm], linewidth=2.4,
                markersize=4, label=f"arm {arm}")
        ax.scatter(xs[-1], ys[-1], zs[-1], color=colors[arm],
                   marker="s", s=75, depthshade=False)

    origin = center if absolute else (0, 0, 0)
    ax.scatter(*origin, color="black", s=100, depthshade=False, label="center")
    # 3軸の表示幅と描画箱の比率を揃え、1格子間隔を同じ尺度で示す。
    ranges = [(min(p[i] for p in all_points), max(p[i] for p in all_points))
              for i in range(3)]
    span = max(high - low for low, high in ranges) + 1
    for setter, (low, high) in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), ranges):
        middle = (low + high) / 2
        setter(middle - span / 2, middle + span / 2)
    ax.set_box_aspect((1, 1, 1))
    for axis, setter in zip("xyz", (ax.set_xlabel, ax.set_ylabel, ax.set_zlabel)):
        setter(axis, labelpad=2)
    arm_length = len(paths[0]) - 1
    coordinate_label = "Absolute coordinates" if absolute else "Center-relative coordinates"
    ax.set_title(f"Molecule {molecule_id} (type {molecule_type}), "
                 f"4 arms, {arm_length} steps/arm\n{coordinate_label}")
    ax.view_init(elev=24, azim=-55)
    ax.legend(loc="upper left")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"Saved: {output}")
    print(f"Center: {center}, type: {molecule_type}, arm length: {arm_length}")
    if show:
        plt.show()  # マウスで回転できる。GUIが使える環境でのみ指定する。
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="demo_moleculeの軌跡CSV")
    parser.add_argument("--output", type=Path, default=Path("data/molecule.png"))
    parser.add_argument("--molecule", type=int, default=0, help="描画対象の分子ID")
    parser.add_argument("--dpi", type=int, default=300, help="画像の解像度")
    parser.add_argument("--absolute", action="store_true", help="絶対座標で表示")
    parser.add_argument("--show", action="store_true", help="保存後に3D表示を開く")
    args = parser.parse_args()
    if args.dpi <= 0 or args.molecule < 0:
        parser.error("dpiは正、moleculeは0以上にしてください")
    try:
        center, kind, paths = read_molecule(args.csv, args.molecule)
        draw_molecule(center, kind, paths, args.molecule, args.output,
                      args.dpi, args.absolute, args.show)
    except (ValueError, OSError, ImportError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
