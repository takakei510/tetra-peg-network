#!/usr/bin/env python3
"""Compare C4 and C6 in completed tetra-peg-network samples."""
import argparse
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
import csv
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path, PurePosixPath
import re
import statistics
import tarfile
import tempfile
from zipfile import ZipFile

LENGTHS = (4, 6)
PHI_TICKS = range(10, 71)
FIELDS = ("bonds.csv.gz", "cluster-members.csv.gz")
ENGINE_SHA256 = "d16eede57d87d37ae3578454847ea6655b4124bb6fdd19be35804f7dbbdd9faa"
CASE_ID = re.compile(r"L\d+_phi\d+\.\d+_i\d{6}")


class SplitArchive:
    def __init__(self, directory):
        parts = sorted(directory.glob("h2.tar.part-*"))
        if not parts or [p.name for p in parts] != [
            f"h2.tar.part-{i:03d}" for i in range(len(parts))
        ]:
            raise ValueError("Missing or misnumbered archive parts")
        self.parts = iter(parts)
        self.current = None

    def read(self, size):
        chunks = []
        while size:
            if self.current is None:
                part = next(self.parts, None)
                if part is None:
                    break
                self.current = part.open("rb")
            chunk = self.current.read(size)
            if chunk:
                chunks.append(chunk)
                size -= len(chunk)
            else:
                self.close()
        return b"".join(chunks)

    def close(self):
        if self.current is not None:
            self.current.close()
            self.current = None


def read_samples(path, phis, replicates):
    if path.suffix.lower() == ".zip":
        with ZipFile(path) as bundle:
            content = bundle.read("samples.csv").decode("utf-8-sig")
    else:
        content = path.read_text(encoding="utf-8-sig")
    selected = [
        row for row in csv.DictReader(content.splitlines())
        if int(row["L"]) == 128
        or (int(row["L"]) in (32, 64) and round(100 * float(row["target_phi"])) in phis)
    ]
    groups = defaultdict(list)
    for row in selected:
        if not CASE_ID.fullmatch(row["case_id"]):
            raise ValueError("Invalid case ID")
        groups[int(row["L"]), round(100 * float(row["target_phi"]))].append(row)
    required = {(128, phi) for phi in PHI_TICKS}
    required.update((L, phi) for L in (32, 64) for phi in phis)
    if groups.keys() != required:
        raise ValueError("The sample manifest does not contain all requested conditions")
    result = {}
    for rows in groups.values():
        if len(rows) != 100 or len({r["seed"] for r in rows}) != 100:
            raise ValueError("Expected 100 distinct seeds per condition in the manifest")
        for row in sorted(rows, key=lambda sample: sample["case_id"])[:replicates]:
            if row["case_id"] in result:
                raise ValueError("Duplicate case ID")
            result[row["case_id"]] = row
    return result


def count_sample(engine, study, sample):
    rec = engine.inspect_case(study, sample)
    clusters = rec["clusters"]
    largest, = (row for row in clusters if row["is_largest"])
    if largest["cluster_id"] != int(sample["selected_largest_id"]):
        raise ValueError(f"{sample['case_id']}: inconsistent largest component")
    result = {key: sample[key] for key in ("case_id", "L", "target_phi", "actual_phi", "seed", "N", "N1")}
    n, n1 = int(result["N"]), int(result["N1"])
    for k in LENGTHS:
        name = f"C{k}_mol"
        total = sum(row[name] for row in clusters)
        inside = largest[name]
        outside = total - inside
        in100 = 100 * inside / n1
        out100 = 100 * outside / (n - n1) if n > n1 else None
        result.update({
            name: total, f"C{k}_largest": inside, f"C{k}_outside": outside,
            f"C{k}_all_per_100": 100 * total / n,
            f"C{k}_largest_per_100": in100, f"C{k}_outside_per_100": out100,
            f"C{k}_paired_difference": in100 - out100 if out100 is not None else None,
            f"C{k}_largest_share": inside / total if total else None,
        })
    result["largest_molecule_share"] = n1 / n
    result.update({f"{key}_sha256": value for key, value in rec["input_sha256"].items()})
    return result


def archive_samples(archive, wanted, temporary):
    pending = {}
    completed = set()
    stream = SplitArchive(archive)
    try:
        with tarfile.open(fileobj=stream, mode="r|") as tar:
            for member in tar:
                parts = PurePosixPath(member.name).parts
                if len(parts) < 5 or not member.isfile():
                    continue
                marker, case, run, run_case, filename = parts[-5:]
                if marker != "cases" or run != "runs" or run_case != case:
                    continue
                if case not in wanted or filename not in FIELDS:
                    continue
                if case in completed:
                    raise ValueError(f"{case}: duplicate {filename} in archive")
                folder = temporary / "cases" / case / "runs" / case
                folder.mkdir(parents=True, exist_ok=True)
                seen = pending.setdefault(case, set())
                if filename in seen:
                    raise ValueError(f"{case}: duplicate {filename} in archive")
                with tar.extractfile(member) as src, (folder / filename).open("wb") as target:
                    while True:
                        data = src.read(1024 * 1024)
                        if not data:
                            break
                        target.write(data)
                seen.add(filename)
                if set(FIELDS) == seen:
                    yield case
                    for name in FIELDS:
                        (folder / name).unlink()
                    for directory in (folder, folder.parent, folder.parent.parent):
                        directory.rmdir()
                    del pending[case]
                    completed.add(case)
                    if completed == wanted:
                        return
    finally:
        stream.close()


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[int(row["L"]), float(row["target_phi"])].append(row)
    metrics = ["largest_molecule_share"]
    for k in LENGTHS:
        metrics.extend(f"C{k}_{suffix}" for suffix in (
            "mol", "largest", "outside", "all_per_100", "largest_per_100",
            "outside_per_100", "paired_difference", "largest_share",
        ))
    for (L, phi), cases in sorted(groups.items()):
        item = {"L": L, "target_phi": phi, "actual_phi": float(cases[0]["actual_phi"])}
        for metric in metrics:
            values = [row[metric] for row in cases if row[metric] is not None]
            item[f"{metric}_n"] = len(values)
            item[f"{metric}_mean"] = statistics.mean(values) if values else None
            item[f"{metric}_se"] = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else None
        yield item


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(output, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def series(ax, rows, key, label):
        x = [r["target_phi"] for r in rows]
        y = [np.nan if r[f"{key}_mean"] is None else r[f"{key}_mean"] for r in rows]
        err = [np.nan if r[f"{key}_se"] is None else r[f"{key}_se"] for r in rows]
        ax.errorbar(x, y, yerr=err, marker="o", markersize=3, lw=1.2, label=label)

    l128 = [row for row in summary if row["L"] == 128]
    shared = set.intersection(*({row["target_phi"] for row in summary if row["L"] == L}
                                for L in (32, 64, 128)))
    for name, suffixes, ylabels in (
        ("01_occupancy", ("mol", "all_per_100"), ("Cycles per sample", "Cycles per 100 molecules")),
        ("02_largest_component", ("largest", "largest_per_100"), ("Cycles per sample", "Cycles per 100 molecules")),
        ("03_lattice_size", ("all_per_100",), ("Cycles per 100 molecules",)),
    ):
        fig, axes = plt.subplots(2, len(suffixes), squeeze=False, figsize=(9, 6), layout="constrained")
        for row, k in zip(axes, LENGTHS):
            for ax, suffix, ylabel in zip(row, suffixes, ylabels):
                if name == "03_lattice_size":
                    for L in (32, 64, 128):
                        selected = [r for r in summary if r["L"] == L and r["target_phi"] in shared]
                        series(ax, selected, f"C{k}_{suffix}", f"L = {L}")
                elif name == "02_largest_component":
                    series(ax, l128, f"C{k}_{suffix}", "Inside largest component")
                    series(ax, l128, f"C{k}_{suffix.replace('largest', 'outside')}", "Outside largest component")
                else:
                    series(ax, l128, f"C{k}_{suffix}", "L = 128")
                ax.set(xlabel="Target occupancy", ylabel=ylabel, title=f"C{k}")
                ax.legend(frameon=False)
                ax.grid(alpha=.22)
        fig.savefig(output / f"{name}.png", dpi=190)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True, help="samples.csv or upload_h2_validation.zip")
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--study", type=Path, help="Study directory containing cases/")
    sources.add_argument("--archive", type=Path, help="Directory containing split h2.tar.part-*")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--phis", default="0.10,0.43,0.70",
                        help="Occupancies for L=32 and 64, or all for 0.10 to 0.70")
    parser.add_argument("--plan", action="store_true",
                        help="Validate the sample manifest and print counts without opening bond data")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--case", help="Select one manifest case for an input test")
    args = parser.parse_args()
    if not 1 <= args.replicates <= 100:
        parser.error("--replicates must be 1–100")
    if args.phis.strip().lower() == "all":
        phis = set(PHI_TICKS)
    else:
        try:
            ticks = [Decimal(item.strip()) * 100 for item in args.phis.split(",")]
            phis = {int(tick) for tick in ticks if tick in PHI_TICKS}
            if len(phis) != len(set(ticks)) or not phis:
                raise ValueError("Occupancies must match the sampled grid")
        except (InvalidOperation, ValueError, OverflowError):
            parser.error("--phis must be all or a comma-separated list on the 0.01 grid (0.10 to 0.70)")
    selected = read_samples(args.samples, phis, args.replicates)
    if args.case:
        if args.case not in selected:
            parser.error("--case is outside the selected conditions or replicates")
        selected = {args.case: selected[args.case]}
    if args.plan:
        counts = Counter(int(row["L"]) for row in selected.values())
        for L, count in sorted(counts.items()):
            conditions = {row["target_phi"] for row in selected.values() if int(row["L"]) == L}
            print(f"L={L}: {len(conditions)} occupancies, {count:,} samples")
        print(f"Total: {len(selected):,} samples")
        return
    engine_path = args.repo / "scripts" / "analysis" / "analyze_cycles_46.py"
    if sha256(engine_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest() != ENGINE_SHA256:
        parser.error("Use the verified repository revision 11e6d471bd3b1bf5ec6b9f477aed8db8108e7624")
    spec = importlib.util.spec_from_file_location("cycle46_engine", engine_path)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    args.output.mkdir(parents=True, exist_ok=True)
    cache = args.output / "checkpoint.jsonl"
    identity = {
        "engine": ENGINE_SHA256,
        "runner": sha256(Path(__file__).read_bytes()).hexdigest(),
        "samples": sha256(args.samples.read_bytes()).hexdigest(),
        "source": str((args.study or args.archive).resolve()),
        "parts": sorted((p.name, p.stat().st_size, p.stat().st_mtime_ns)
                        for p in args.archive.glob("h2.tar.part-*")) if args.archive else [],
    }
    identity = json.loads(json.dumps(identity))
    rows = {}
    if cache.exists():
        with cache.open(encoding="utf-8") as handle:
            if json.loads(next(handle)) != identity:
                parser.error("Checkpoint source changed. Use a new output directory.")
            for line in handle:
                row = json.loads(line)
                rows[row["case_id"]] = row
    else:
        cache.write_text(json.dumps(identity) + "\n", encoding="utf-8")
    wanted = selected.keys() - rows.keys()
    with cache.open("a", encoding="utf-8") as handle, tempfile.TemporaryDirectory(dir=args.output) as root:
        study = args.study or Path(root)
        if args.study:
            for case in wanted:
                run = study / "cases" / case / "runs" / case
                for filename in ("bonds.csv", "cluster-members.csv"):
                    engine.raw_path(run, filename)
            for case in rows.keys() & selected.keys():
                run = study / "cases" / case / "runs" / case
                for key, filename in (("bonds", "bonds.csv"), ("members", "cluster-members.csv")):
                    if engine.digest(engine.raw_path(run, filename)) != rows[case][f"{key}_sha256"]:
                        parser.error(f"{case}: source changed. Use a new output directory.")
        cases = sorted(wanted) if args.study else archive_samples(args.archive, wanted, study) if wanted else ()
        for progress, case in enumerate(cases, 1):
            if case in rows:
                raise ValueError(f"Duplicate case: {case}")
            row = count_sample(engine, study, selected[case])
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            rows[case] = row
            if progress == 1 or progress % 100 == 0:
                print(f"{progress}/{len(wanted)} new samples", flush=True)
    missing = selected.keys() - rows.keys()
    if missing:
        raise ValueError(f"{len(missing)} cases missing bonds or members, e.g. {sorted(missing)[:3]}")
    records = [rows[case] for case in sorted(selected)]
    report = list(summarize(records))
    write_csv(args.output / "cycle46_samples.csv", records)
    write_csv(args.output / "cycle46_summary.csv", report)
    if not args.case:
        plot(args.output, report)
    print(f"Completed {len(records)} samples and {len(report)} conditions", flush=True)


if __name__ == "__main__":
    main()
