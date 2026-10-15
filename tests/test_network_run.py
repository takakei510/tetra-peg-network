"""実行入口を通し、CSV・未達保存・結合条件・UFを独立に検証する。Matplotlib不要。"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def run(config, limit=10000):
    result = subprocess.run([sys.executable, "scripts/run_demo.py", "--mode", "network",
                             "--config", str(config), "--max-placement-attempts", str(limit), "--no-plot"],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode in (0, 3), result.stdout + result.stderr
    directory = Path(result.stdout.split("保存先: ", 1)[1].strip())
    manifest = json.loads((directory / "run.json").read_text())
    assert manifest["status"] == ("completed" if result.returncode == 0 else "incomplete")
    assert all((directory / file).is_file() for file in manifest["files"].values())
    return directory, manifest


config = ROOT / "configs/network/preview_phi010_L16_seed3.cfg"
a, manifest = run(config)
b, _ = run(config)
# 実行時刻等ではなく、物理結果のCSVが同じことを確認する。
for name in ("trajectories.csv", "candidates.csv", "bonds.csv", "cluster-members.csv",
             "cluster-summary.csv", "network-summary.csv", "placement-summary.csv"):
    assert (a / name).read_bytes() == (b / name).read_bytes()
# 枝リストは1成立結合1行。ID・末端IDを保持し、代表最大成分だけ抽出する。
assert manifest["graph_analysis_completed"]
graph_summary = json.loads((a / "graph/graph-analysis.json").read_text())
assert graph_summary["largest"] == 2 and len(graph_summary["tied_largest_cluster_ids"]) == 3
assert graph_summary["selected_cluster_id"] == min(graph_summary["tied_largest_cluster_ids"])
assert not graph_summary["plots_generated"]
assert len(rows(a / "graph/largest-members.csv")) == 2
original_bonds = rows(a / "bonds.csv")
graph_edges = rows(a / "graph/molecular-edges.csv")
assert len(graph_edges) == len(original_bonds)
for original, edge in zip(original_bonds, graph_edges):
    assert all(original[k] == edge[k] for k in ("bond_id", "molecule_a", "molecule_b", "endpoint_a", "endpoint_b"))
for file in ("molecular-edges.csv", "largest-members.csv", "largest-edges.csv", "graph-analysis.json"):
    assert (a / "graph" / file).read_bytes() == (b / "graph" / file).read_bytes()
placement = rows(a / "placement-summary.csv")[0]
assert int(placement["placed_molecules"]) == 19 and int(placement["occupied_sites"]) == 399
assert int(placement["attempts"]) == 19 + int(placement["center_occupied"]) + int(placement["trapped"])
assert int(placement["type_A"]) + int(placement["type_B"]) == 19
paths = rows(a / "trajectories.csv")
assert len({r["site"] for r in paths}) == 399
endpoints = {4 * int(r["molecule"]) + int(r["arm"]): r for r in paths if r["step"] == "5"}
assert len(endpoints) == 76
# 全末端対との比較はテスト専用。高速探索の候補漏れ・余分な候補を検出する。
expected = set()
for u, x in endpoints.items():
    for v, y in endpoints.items():
        distance = sum(abs(int(x[k]) - int(y[k])) for k in ("x", "y", "z"))
        if u < v and x["type"] != y["type"] and 1 <= distance <= 2:
            expected.add((u, v))
candidates = rows(a / "candidates.csv")
actual = {(int(r["endpoint_a"]), int(r["endpoint_b"])) for r in candidates}
assert actual == expected
used = set()
graph = [set() for _ in range(19)]
for r in rows(a / "bonds.csv"):
    u, v = int(r["endpoint_a"]), int(r["endpoint_b"])
    assert (u, v) in expected and u not in used and v not in used
    used.update((u, v)); graph[u // 4].add(v // 4); graph[v // 4].add(u // 4)
# UFとは別の探索でクラスタ分割を照合。孤立分子も含める。
members = rows(a / "cluster-members.csv")
assignment = {int(r["molecule"]): int(r["cluster_id"]) for r in members}
for start in range(19):
    reached, pending = {start}, [start]
    while pending:
        for other in graph[pending.pop()]:
            if other not in reached:
                reached.add(other); pending.append(other)
    assert reached == {m for m in range(19) if assignment[m] == assignment[start]}
partial, m = run(config, 1)
assert not m["generation_completed"] and not m["analysis_completed"]
assert m["summary"]["placed_molecules"] == 1 and m["summary"]["occupied_sites"] == 21
assert not (partial / "graph").exists()
assert not (partial / "bonds.csv").exists() and not (partial / "network-summary.csv").exists()
with tempfile.TemporaryDirectory(dir=ROOT / "build") as folder:
    cfg = Path(folder) / "count.cfg"
    cfg.write_text("dim=3\nL=16\nseed=3\ntarget_molecules=3\narm_length=5\na_probability=1\n")
    count, m = run(cfg)
    assert m["target_molecules"] == 3 and m["summary"]["type_counts"] == {"A": 3, "B": 0}
    assert not rows(count / "bonds.csv") and not rows(count / "candidates.csv")
    assert len(rows(count / "cluster-summary.csv")) == 3
    isolated_graph = json.loads((count / "graph/graph-analysis.json").read_text())
    assert isolated_graph["tied_largest_cluster_ids"] == [0, 1, 2]
    assert isolated_graph["largest_bonds"] == 0
    assert rows(count / "graph/largest-members.csv")[0]["molecule"] == "0"
    assert not rows(count / "graph/largest-edges.csv")
    assert rows(count / "placement-summary.csv")[0]["requested_phi"] == ""
    cfg = Path(folder) / "zero.cfg"
    cfg.write_text("dim=3\nL=2\nseed=3\ntarget_molecules=1\narm_length=1\na_probability=1\n")
    zero, m = run(cfg, 7)
    assert m["summary"]["placed_molecules"] == 0 and not m["plot_generated"]
    assert not rows(zero / "trajectories.csv") and not (zero / "bonds.csv").exists()
    assert rows(zero / "placement-summary.csv")[0]["trapped"] == "7"
print("network integration tests passed (reproducibility, candidate oracle, bonds, BFS, partial/zero retention)")
