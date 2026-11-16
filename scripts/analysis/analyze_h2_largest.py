"""Measure exact double-bonded molecule pairs in the largest component."""
import argparse
from collections import defaultdict
from contextlib import closing
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import statistics
import tarfile

from compare_cycles_46 import PHI_TICKS, SplitArchive, read_samples, write_csv

FILENAME = "clusters.csv.gz"
METRICS = ("H2_largest", "H2_outside", "largest_endpoint_pct",
           "outside_endpoint_pct", "largest_H2_share_pct", "largest_molecule_share_pct")


def measure(zipped, sample):
    case_id = sample["case_id"]
    n, n1, total_h2 = (int(sample[k]) for k in ("N", "N1", "H2"))
    largest = None
    sums = [0, 0, 0]
    ids = set()
    with gzip.GzipFile(fileobj=zipped) as data, io.TextIOWrapper(data, encoding="utf-8-sig") as text:
        for row in csv.DictReader(text):
            cid, size, bonds, h2 = (int(row[k]) for k in ("cluster_id", "size", "bonds", "h2"))
            if cid < 0 or cid in ids or not (size > 0 and 0 <= 2*h2 <= bonds <= 2*size):
                raise ValueError(f"{case_id}: invalid cluster row {cid}")
            ids.add(cid)
            sums[0] += size
            sums[1] += bonds
            sums[2] += h2
            if largest is None or (-size, cid) < (-largest[1], largest[0]):
                largest = cid, size, h2
    if (sums != [n, int(sample["bonds"]), total_h2] or len(ids) != int(sample["clusters"])
            or largest is None or largest[:2] != (int(sample["selected_largest_id"]), n1)):
        raise ValueError(f"{case_id}: cluster totals or largest component disagree with samples.csv")
    inside = largest[2]
    outside = total_h2 - inside
    return dict(
        case_id=case_id, L=int(sample["L"]), target_phi=float(sample["target_phi"]),
        actual_phi=float(sample["actual_phi"]), seed=int(sample["seed"]),
        N=n, N1=n1, H2=total_h2, largest_id=largest[0],
        H2_largest=inside, H2_outside=outside,
        largest_endpoint_pct=100*inside/n1,
        outside_endpoint_pct=100*outside/(n-n1) if n>n1 else None,
        largest_H2_share_pct=100*inside/total_h2 if total_h2 else None,
        largest_molecule_share_pct=100*n1/n,
    )


def from_archive(directory, selected):
    with closing(SplitArchive(directory)) as stream, tarfile.open(fileobj=stream, mode="r|") as archive:
        for member in archive:
            parts = PurePosixPath(member.name).parts
            if not member.isfile() or len(parts)<3 or parts[-3] != "cases" or parts[-1] != FILENAME:
                continue
            case_id = parts[-2]
            if case_id in selected:
                with archive.extractfile(member) as data:
                    yield case_id, measure(data, selected[case_id])


def from_study(study, selected):
    for case_id, sample in sorted(selected.items()):
        with (study/"cases"/case_id/FILENAME).open("rb") as data:
            yield case_id, measure(data, sample)


def summary(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["L"], row["target_phi"]].append(row)
    result = []
    for (size, phi), cases in sorted(groups.items()):
        item = dict(L=size, target_phi=phi, n_samples=len(cases))
        for metric in METRICS:
            values = [row[metric] for row in cases if row[metric] is not None]
            n = len(values)
            item.update({f"{metric}_n":n, f"{metric}_mean":statistics.mean(values) if n else None,
                         f"{metric}_se":statistics.stdev(values)/math.sqrt(n) if n>1 else None})
        result.append(item)
    return result


def plot(output, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [
        ("H2_largest", "H₂ pairs in largest component"),
        ("largest_endpoint_pct", "All ends used by H₂ inside [%]"),
        ("largest_H2_share_pct", "Share of all H₂ pairs inside [%]"),
        ("largest_molecule_share_pct", "Share of all molecules inside [%]"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for ax, (key, title) in zip(axes.flat, panels):
        for size in sorted({r["L"] for r in rows}):
            group = [r for r in rows if r["L"]==size]
            means = [float("nan") if r[f"{key}_mean"] is None else r[f"{key}_mean"] for r in group]
            errors = [float("nan") if r[f"{key}_se"] is None else r[f"{key}_se"] for r in group]
            ax.errorbar([r["target_phi"] for r in group], means, yerr=errors,
                        label=f"L = {size}", marker=".", lw=1.2)
        ax.set(xlabel="Target occupancy", title=title)
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    fig.savefig(output/"h2_largest_vs_occupancy.png", dpi=200)
    plt.close(fig)


def save_checkpoint(path, identity, rows):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps({"identity":identity, "rows":rows}), encoding="utf-8")
    tmp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--study", type=Path)
    source.add_argument("--archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=100, choices=range(1,101), metavar="1-100")
    parser.add_argument("--case")
    args = parser.parse_args()
    selected = read_samples(args.samples, set(PHI_TICKS), args.replicates)
    if args.case:
        if args.case not in selected:
            parser.error("--case is outside the selected sample set")
        selected = {args.case:selected[args.case]}
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    identity = dict(samples=hashlib.sha256(args.samples.read_bytes()).hexdigest(),
                    code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    directory=str((args.archive or args.study).resolve()))
    parts = sorted(args.archive.glob("h2.tar.part-*")) if args.archive else []
    if args.archive:
        identity["parts"] = [[p.name,p.stat().st_size,p.stat().st_mtime_ns] for p in parts]
    checkpoint = output/"checkpoint.json"
    rows = {}
    if checkpoint.is_file():
        saved = json.loads(checkpoint.read_text(encoding="utf-8"))
        if saved["identity"] != identity:
            parser.error("Checkpoint inputs changed. Specify a new output directory.")
        rows = saved["rows"]
    pending = {case:sample for case,sample in selected.items() if case not in rows}
    if pending:
        source = from_archive(args.archive, pending) if args.archive else from_study(args.study, pending)
        try:
            done = len(selected) - len(pending)
            with closing(source):
                for case, row in source:
                    if case in rows:
                        raise ValueError(f"{case}: duplicate cluster table")
                    rows[case] = row
                    done += 1
                    if done % 500 == 0 or done == len(selected):
                        save_checkpoint(checkpoint, identity, rows)
                    if done == 1 or done % 100 == 0 or done == len(selected):
                        print(f"{done}/{len(selected)} samples", flush=True)
                    if done == len(selected):
                        break
        finally:
            if rows:
                save_checkpoint(checkpoint, identity, rows)
    missing = selected.keys()-rows.keys()
    if missing:
        raise ValueError(f"Missing cluster tables: {len(missing)}, e.g. {sorted(missing)[:3]}")
    records = [rows[case] for case in sorted(selected)]
    grouped = summary(records)
    write_csv(output/"h2_largest_samples.csv", records)
    write_csv(output/"h2_largest_summary.csv", grouped)
    if not args.case:
        plot(output, grouped)
    print(f"Completed {len(records)} samples and {len(grouped)} conditions", flush=True)


if __name__ == "__main__":
    main()
