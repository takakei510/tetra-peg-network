#!/usr/bin/env python3
"""保存済みの完了試料から分子多重グラフと代表最大成分を抽出する。

python scripts/visualization/analyze_graph.py data/runs/<run> [--no-plot]
図の頂点位置は表示用のA/B二列配置であり、実空間距離ではない。
"""
import argparse
from collections import Counter, defaultdict, deque
import csv
import json
from pathlib import Path
import sys

from plot_molecules import draw_molecules, read_bonds, read_molecules

EDGE_FIELDS = ("bond_id", "molecule_a", "molecule_b", "endpoint_a", "endpoint_b", "cluster_id")
MEMBER_FIELDS = ("molecule", "cluster_id", "cluster_size", "bond_degree")


def read_rows(path, fields):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not set(fields).issubset(reader.fieldnames or []):
            raise ValueError(f"{path.name}: 必要な列がありません")
        return list(reader)


def validate_graph(ids, bonds, members):
    """枝リストからBFSでUFの分割・次数を照合。多重枝を捨てない。"""
    ids = set(ids)
    if not ids:
        raise ValueError("空試料は解析できません")
    assignments, groups = {}, defaultdict(set)
    for row in members:
        m, cid = int(row["molecule"]), int(row["cluster_id"])
        if m not in ids or m in assignments or cid < 0:
            raise ValueError("所属CSVの分子ID・クラスタIDが不正です")
        assignments[m] = cid
        groups[cid].add(m)
    if set(assignments) != ids:
        raise ValueError("所属CSVと軌跡CSVの分子一覧が一致しません")
    adjacency = {m: [] for m in ids}
    for row in bonds:
        u, v = int(row["molecule_a"]), int(row["molecule_b"])
        if u not in ids or v not in ids or u == v:
            raise ValueError("枝の分子IDが不正です")
        if assignments[u] != assignments[v]:
            raise ValueError("成立枝が異なるクラスタを結んでいます")
        adjacency[u].append(v)
        adjacency[v].append(u)
    for row in members:
        m = int(row["molecule"])
        if int(row["cluster_size"]) != len(groups[assignments[m]]):
            raise ValueError("cluster_sizeと所属一覧が一致しません")
        if int(row["bond_degree"]) != len(adjacency[m]):
            raise ValueError("bond_degreeと枝リストが一致しません")
    visited = set()
    for start in sorted(ids):
        if start in visited:
            continue
        reached, pending = {start}, deque([start])
        while pending:
            for other in adjacency[pending.popleft()]:
                if other not in reached:
                    reached.add(other)
                    pending.append(other)
        if reached != groups[assignments[start]]:
            raise ValueError("UFのクラスタ所属と枝リストのBFSが一致しません")
        visited.update(reached)
    largest_size = max(map(len, groups.values()))
    tied = sorted(cid for cid, group in groups.items() if len(group) == largest_size)
    return assignments, groups, tied


def write_rows(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in rows)


def draw_graph(molecules, edges, output, title, dpi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch

    fig, ax = plt.subplots(figsize=(9, 7))
    positions = {}
    colors = {"A": "#4472C4", "B": "#C0504D"}
    for kind, x in (("A", -1), ("B", 1)):
        ids = sorted(m for m, item in molecules.items() if item[1] == kind)
        for i, m in enumerate(ids):
            positions[m] = (x, 0 if len(ids) == 1 else 1 - 2 * i / (len(ids) - 1))
    pairs = Counter(tuple(sorted((int(row["molecule_a"]), int(row["molecule_b"]))))
                    for row in edges)
    for (u, v), count in sorted(pairs.items()):
        for i in range(count):
            curve = (i - (count - 1) / 2) * 0.22
            ax.add_patch(FancyArrowPatch(positions[u], positions[v], arrowstyle="-",
                         connectionstyle=f"arc3,rad={curve}", linewidth=1.6,
                         color="#23854B", alpha=0.85, shrinkA=9, shrinkB=9, zorder=1))
    for kind in ("A", "B"):
        ids = sorted(m for m, item in molecules.items() if item[1] == kind)
        if ids:
            ax.scatter([positions[m][0] for m in ids], [positions[m][1] for m in ids],
                       s=180 if len(molecules) <= 50 else 35, c=colors[kind],
                       edgecolors="white", zorder=2, label=f"Type {kind} ({len(ids)})")
    if len(molecules) <= 50:
        for m, (x, y) in positions.items():
            ax.annotate(str(m), (x, y), xytext=(-18 if x < 0 else 12, 0),
                        textcoords="offset points", va="center", fontsize=9)
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.5, 1.5)
    ax.set_axis_off()
    ax.set_title(f"{title}\n{len(molecules)} molecules, {len(edges)} bonds")
    ax.legend(loc="upper right")
    fig.text(0.5, 0.025, "A/B columns; positions are schematic, not spatial coordinates.\n"
             "Parallel curves show separate bonds. Node labels are original molecule IDs.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(output, dpi=dpi, facecolor="white")
    plt.close(fig)


def analyze(run, output, no_plot=False, dpi=200, allow_running=False):
    if output.exists():
        raise ValueError("出力先が既に存在します。別の--output-dirを指定してください")
    manifest_path = run / "run.json"
    if not manifest_path.is_file():
        raise ValueError("run.jsonのある実行フォルダを指定してください")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    allowed = {"completed", "running"} if allow_running else {"completed"}
    if (manifest.get("status") not in allowed or
            manifest.get("generation_completed") is not True or
            manifest.get("mode") not in ("clusters", "network")):
        raise ValueError("生成・クラスタ解析が完了したclusters/network実行だけを対象にします")
    if manifest.get("mode") == "network" and manifest.get("analysis_completed") is not True:
        raise ValueError("クラスタ解析が完了していません")
    molecules, occupied = read_molecules(run / "trajectories.csv")
    bond_points = read_bonds(run / "bonds.csv", molecules)
    bonds = read_rows(run / "bonds.csv", EDGE_FIELDS[:-1])
    members = read_rows(run / "cluster-members.csv", MEMBER_FIELDS)
    assignments, groups, tied = validate_graph(molecules, bonds, members)
    chosen = tied[0]
    selected = groups[chosen]
    edges = [{**row, "cluster_id": assignments[int(row["molecule_a"])]} for row in bonds]
    largest_edges = [row for row in edges if int(row["cluster_id"]) == chosen]
    chosen_molecules = {m: item for m, item in molecules.items() if m in selected}
    selected_points = [points for row, points in zip(edges, bond_points)
                       if int(row["cluster_id"]) == chosen]
    net = read_rows(run / "network-summary.csv", ("molecules", "bonds", "clusters", "largest"))
    expected = {"molecules": len(molecules), "bonds": len(edges),
                "clusters": len(groups), "largest": len(selected)}
    if len(net) != 1 or any(int(net[0][key]) != value for key, value in expected.items()):
        raise ValueError("network-summary.csvと枝・所属の集計が一致しません")
    conditions = manifest.get("reported_conditions", {})
    phi = s1 = None
    if conditions.get("L"):
        L = int(conditions["L"])
        if L <= 0:
            raise ValueError("Lは正にしてください")
        phi = occupied / L**3
        s1 = phi * len(selected) / len(molecules)
    summary = {"schema_version": 1, **expected, "selected_cluster_id": chosen,
               "tied_largest_cluster_ids": tied, "selection_rule": "smallest_cluster_id",
               "largest_bonds": len(largest_edges), "largest_fraction": len(selected) / len(molecules),
               "actual_phi": phi, "largest_occupied_fraction": s1,
               "layout": "deterministic_type_columns_id_order",
               "spatial_coordinates": "original_absolute_open_boundary",
               "validation": "bond_constraints_and_bfs_partition_and_degrees",
               "plots_generated": False,
               "files": ["molecular-edges.csv", "largest-members.csv", "largest-edges.csv"]}
    output.mkdir(parents=True)
    summary_path = output / "graph-analysis.json"
    try:
        write_rows(output / "molecular-edges.csv", EDGE_FIELDS, edges)
        write_rows(output / "largest-members.csv", MEMBER_FIELDS,
                   [row for row in members if int(row["molecule"]) in selected])
        write_rows(output / "largest-edges.csv", EDGE_FIELDS, largest_edges)
        if not no_plot:
            draw_graph(molecules, edges, output / "network.png", "Molecular multigraph", dpi)
            draw_graph(chosen_molecules, largest_edges, output / "largest-network.png",
                       f"Representative largest component (ID {chosen}, {len(tied)} tied)", dpi)
            chosen_occupied = sum(1 + sum(len(path)-1 for path in item[2].values())
                                  for item in chosen_molecules.values())
            draw_molecules(chosen_molecules, chosen_occupied, output / "largest-molecules.png",
                           dpi, bonds=selected_points)
            summary["plots_generated"] = True
            summary["files"].extend(["network.png", "largest-network.png", "largest-molecules.png"])
        summary["status"] = "completed"
    except Exception:
        summary["status"] = "failed"
        raise
    finally:
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output-dir", type=Path, help="既定は実行フォルダのgraph/。既存先の上書きは禁止")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("dpiは正にしてください")
    try:
        summary = analyze(args.run, args.output_dir or args.run / "graph", args.no_plot, args.dpi)
        print(f"Graph saved: {args.output_dir or args.run / 'graph'}")
        print(f"Largest: {summary['largest']} molecules, {summary['largest_bonds']} bonds, "
              f"cluster_id={summary['selected_cluster_id']}, ties={summary['tied_largest_cluster_ids']}")
    except (ValueError, OSError, ImportError, TypeError, KeyError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
