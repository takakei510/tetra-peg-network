#!/usr/bin/env python3
"""デモ1回の設定・結果・コードを同じフォルダに保存する実行入口。"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import uuid
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_output(*args):
    # Gitがない環境でも実行できる。再現に必要なコードは別途ZIPで残す。
    try:
        result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                                text=True, check=False)
        return result.stdout.strip() if result.returncode == 0 else None
    except OSError:
        return None


def snapshot_sources(destination):
    # build/data/.venvを含めず、コンパイル・描画・実行に使う入力だけを控える。
    # 未コミットの変更も残るので、コミットIDだけより確実にコードを追跡できる。
    files = [ROOT / "Makefile", ROOT / "scripts" / "visualization" / "requirements.txt"]
    for directory, suffixes in [("src", {".c"}), ("include", {".h"}),
                                ("scripts", {".py"})]:
        files.extend(p for p in (ROOT / directory).rglob("*")
                     if p.is_file() and p.suffix in suffixes)
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            archive.write(path, path.relative_to(ROOT).as_posix())


def csv_summary(path):
    # CSVから実現した分子数・型・占有点数を集計する。目標値と混同しない。
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"molecule", "type", "arm", "step", "site", "x", "y", "z"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("軌跡CSVの列が不足しています")
        rows = list(reader)
    molecules = {row["molecule"]: row["type"] for row in rows}
    return {"placed_molecules": len(molecules),
            "type_counts": {kind: list(molecules.values()).count(kind)
                            for kind in ("A", "B")},
            "occupied_sites": len({row["site"] for row in rows}),
            "trajectory_rows": len(rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=["single", "two", "candidates", "bonds", "clusters", "network"])
    parser.add_argument("--config", required=True, type=Path,
                        help="実行条件のcfg。相対パスは現在の作業ディレクトリ基準")
    parser.add_argument("--no-plot", action="store_true",
                        help="CSVのみ保存する（Matplotlib不要）")
    parser.add_argument("--max-placement-attempts", type=int, default=None,
                        help="network専用の全分子合計試行上限（既定10000回）")
    args = parser.parse_args()
    if args.max_placement_attempts is not None and (args.mode != "network" or args.max_placement_attempts <= 0):
        parser.error("試行上限はnetworkモードで正の整数を指定してください")
    max_attempts = args.max_placement_attempts if args.max_placement_attempts is not None else 10000
    config = args.config.resolve()
    if not config.is_file():
        parser.error(f"設定ファイルがありません: {config}")

    # UTCを明記し、同時実行も衝突しないIDにする。既存結果は一切上書きしない。
    started = datetime.now(timezone.utc)
    run_id = started.strftime("%Y%m%dT%H%M%S.%fZ") + "_" + args.mode + "_" + uuid.uuid4().hex[:8]
    run = ROOT / "data" / "runs" / run_id
    run.mkdir(parents=True, exist_ok=False)
    (run / "config.cfg").write_bytes(config.read_bytes())
    binary_name = {"single": "demo_molecule", "two": "demo_two_molecules",
                   "candidates": "demo_candidates", "bonds": "demo_bonds", "clusters": "demo_clusters", "network": "demo_network"}[args.mode]
    plot_name = "plot_molecule.py" if args.mode == "single" else "plot_molecules.py"
    binary = ROOT / "build" / binary_name
    command = [str(binary), str(run / "config.cfg")]
    if args.mode in ("candidates", "bonds", "clusters", "network"):
        command.append(str(run / "candidates.csv"))
    if args.mode in ("bonds", "clusters", "network"):
        command.append(str(run / "bonds.csv"))
    if args.mode in ("clusters", "network"):
        command.extend(str(run / name) for name in
                       ("cluster-members.csv", "cluster-summary.csv", "network-summary.csv"))
    if args.mode == "network":
        command.extend([str(run / "placement-summary.csv"), str(max_attempts)])
    manifest = {"schema_version": 2 if args.mode == "network" else 1, "run_id": run_id,
                "started_at_utc": started.isoformat(), "mode": args.mode,
                "config_source": str(config), "status": "running",
                "target_molecules": 1 if args.mode == "single" else 2,
                "generation_command": command, "python_version": sys.version,
                "git_commit": git_output("rev-parse", "HEAD"),
                "git_status": git_output("status", "--porcelain"),
                "conditions_source": "metadata.log emitted by the C demo",
                "files": {"config": "config.cfg", "trajectories": "trajectories.csv",
                          "metadata": "metadata.log", "source": "source.zip",
                          "build_log": "build.log"}}

    if args.mode in ("candidates", "bonds", "clusters", "network"):
        manifest["files"]["candidates"] = "candidates.csv"
    if args.mode in ("bonds", "clusters", "network"):
        manifest["files"]["bonds"] = "bonds.csv"
    if args.mode in ("clusters", "network"):
        manifest["files"].update({"cluster_members": "cluster-members.csv",
                                  "cluster_summary": "cluster-summary.csv",
                                  "network_summary": "network-summary.csv"})

    if args.mode == "network":
        # NはCの専用読込を正本にする。Pythonで占有率の換算を重複させない。
        manifest["target_molecules"] = None
        manifest["max_placement_attempts"] = max_attempts
        manifest["files"]["placement_summary"] = "placement-summary.csv"

    def save_manifest():
        # 保存途中で途切れても、直前に完成したJSONを残す。
        temporary = run / "run.json.tmp"
        temporary.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                             encoding="utf-8")
        temporary.replace(run / "run.json")

    save_manifest()
    try:
        snapshot_sources(run / "source.zip")
        with (run / "build.log").open("w") as log:
            # -Bで現行ソースから確実に再ビルドし、古い実行ファイルの混入を防ぐ。
            build = subprocess.run(["make", "-B", f"build/{binary_name}"], cwd=ROOT,
                                   stdout=log, stderr=subprocess.STDOUT, check=False)
        manifest["build_exit_code"] = build.returncode
        if build.returncode:
            raise RuntimeError("ビルドに失敗しました（build.log参照）")
        manifest["binary_sha256"] = sha256(binary)
        with (run / "trajectories.csv").open("w") as output, \
                (run / "metadata.log").open("w") as errors:
            result = subprocess.run(command, cwd=ROOT, stdout=output, stderr=errors,
                                    check=False)
        manifest["generation_exit_code"] = result.returncode
        # 固定条件もCの出力から記録する。Python側に腕長等の定数を重複させない。
        metadata = (run / "metadata.log").read_text(encoding="utf-8")
        manifest["reported_conditions"] = dict(token.split("=", 1)
                                                for token in metadata.split() if "=" in token)
        incomplete = args.mode == "network" and result.returncode == 3
        if result.returncode and not incomplete:
            raise RuntimeError("生成・解析が完了しませんでした（metadata.log参照）")
        if args.mode == "network":
            with (run / "placement-summary.csv").open(newline="", encoding="utf-8") as stream:
                placement = next(csv.DictReader(stream))
            manifest["placement"] = placement
            manifest["target_molecules"] = int(placement["target_molecules"])
            # 未達は実行エラーと区別。未実行の解析CSVを成果物一覧に含めない。
            if incomplete:
                for key in ("candidates", "bonds", "cluster_members", "cluster_summary", "network_summary"):
                    manifest["files"].pop(key, None)
        manifest["summary"] = csv_summary(run / "trajectories.csv")
        if not incomplete and manifest["summary"]["placed_molecules"] != manifest["target_molecules"]:
            raise RuntimeError("CSVの分子数が目標と一致しません")
        manifest["generation_completed"] = not incomplete
        if args.mode == "network":
            manifest["analysis_completed"] = not incomplete
        manifest["graph_analysis_completed"] = False
        if args.mode in ("clusters", "network") and not incomplete:
            # 保存済みの物理結果からの後処理。乱数とCの生成処理は変更しない。
            save_manifest()
            sys.path.insert(0, str(ROOT / "scripts" / "visualization"))
            from analyze_graph import analyze
            try:
                analyze(run, run / "graph", no_plot=args.no_plot, allow_running=True)
                manifest["graph_analysis_completed"] = True
            finally:
                if (run / "graph").is_dir():
                    for path in sorted((run / "graph").iterdir()):
                        manifest["files"]["graph_" + path.stem.replace("-", "_")] = path.relative_to(run).as_posix()
        manifest["plot_generated"] = False
        if not args.no_plot and manifest["summary"]["placed_molecules"]:
            plot_command = [sys.executable, str(ROOT / "scripts" / "visualization" / plot_name),
                            str(run / "trajectories.csv"), "--output", str(run / "molecules.png")]
            if args.mode == "candidates":
                plot_command.extend(["--candidates", str(run / "candidates.csv")])
            if args.mode in ("bonds", "clusters", "network") and not incomplete:
                plot_command.extend(["--bonds", str(run / "bonds.csv")])
            if incomplete:
                plot_command.append("--incomplete")
            manifest["plot_command"] = plot_command
            manifest["files"]["plot_log"] = "plot.log"
            with (run / "plot.log").open("w") as log:
                plotted = subprocess.run(plot_command, cwd=ROOT, stdout=log,
                                         stderr=subprocess.STDOUT, check=False)
            manifest["plot_exit_code"] = plotted.returncode
            if plotted.returncode:
                raise RuntimeError("描画に失敗しました。CSVは保存済みです（plot.log参照）")
            manifest["plot_generated"] = True
            manifest["files"]["plot"] = "molecules.png"
        # 未達の図は完成済み軌跡のみ。結合済みの完成試料と区別する。
        manifest["status"] = "incomplete" if incomplete else "completed"
        exit_code = 3 if incomplete else 0
        if incomplete:
            print(f"試行上限に到達: 目標{manifest['target_molecules']}分子、完成{manifest['summary']['placed_molecules']}分子。完成済み軌跡を保存しました。", file=sys.stderr)
    except (OSError, ValueError, RuntimeError, ImportError, TypeError, KeyError, csv.Error) as error:
        # 不完全な実行も残し、成功データと識別できるようにする。
        manifest["status"] = "failed"
        manifest["error"] = str(error)
        print(error, file=sys.stderr)
        exit_code = 1
    finally:
        manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["sha256"] = {name: sha256(run / name)
                              for name in manifest["files"].values() if (run / name).is_file()}
        save_manifest()
        print(f"status={manifest['status']}\n保存先: {run}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
