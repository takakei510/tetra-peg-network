#!/usr/bin/env python3
"""Generate repeated AB-SAW samples and summarize their resistor response."""

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
from uuid import uuid4

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import analyze_resistance as engine


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


@contextmanager
def output_lock(root):
    """Use an OS lock that is released even when a notebook process stops."""
    with (root / ".study.lock").open("a+b") as stream:
        stream.write(b"0")
        stream.flush()
        stream.seek(0)
        if os.name == "nt":
            import msvcrt
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            stream.seek(0)
            if os.name == "nt":
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def load_config(path):
    config = json.loads(path.read_text(encoding="utf-8-sig"))
    expected = {"schema_version", "sizes", "phi_percent", "seeds", "arm_length", "a_probability",
                "attempts_per_molecule", "minimum_attempts", "analysis"}
    if set(config) != expected or config["schema_version"] != 1:
        raise ValueError("Unexpected study configuration keys/version")
    for field, low, high in (("sizes", 2, 1290), ("phi_percent", 1, 99), ("seeds", 0, 2**64 - 1)):
        values = config[field]
        if not values or any(type(v) is not int or not low <= v <= high for v in values):
            raise ValueError(f"Invalid {field}")
        if len(set(values)) != len(values):
            raise ValueError(f"Duplicate {field}")
    for field in ("arm_length", "attempts_per_molecule", "minimum_attempts"):
        if type(config[field]) is not int or config[field] <= 0:
            raise ValueError(f"Invalid {field}")
    if not 0 <= config["a_probability"] <= 1:
        raise ValueError("Invalid A probability")
    return config


def cases(config):
    for size in sorted(config["sizes"]):
        for percent in sorted(config["phi_percent"]):
            for seed in config["seeds"]:
                yield {"case_id": f"L{size}_phi{percent / 100:.2f}_seed{seed}",
                       "L": size, "target_phi": percent / 100, "seed": seed}


def extras(cores, bonds, branches, settings, summary):
    result = {"active_fraction_1e12": None, "active_fraction_1e10": None,
              "active_fraction_1e8": None, "power_participation_fraction": None,
              "largest_fraction": None, "terminal_relative_error": None}
    indices = {molecule: i for i, molecule in enumerate(cores)}
    left = np.array([indices[a] for _, a, _ in bonds], dtype=int)
    right = np.array([indices[b] for _, _, b in bonds], dtype=int)
    adjacency = coo_matrix((np.ones(2 * len(bonds)), (np.r_[left, right], np.r_[right, left])),
                           shape=(len(cores), len(cores))).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    result["largest_fraction"] = float(np.bincount(labels).max() / len(cores))
    if summary["status"] != "connected":
        return result
    current_scale = (settings["voltage_high"] - settings["voltage_low"]) / settings["branch_resistance"]
    normalized = np.array([abs(branch[3]) / current_scale for branch in branches])
    for exponent in (12, 10, 8):
        result[f"active_fraction_1e{exponent}"] = float(np.mean(normalized > 10.0**-exponent))
    # Scaling the powers first avoids overflow and leaves the participation ratio unchanged.
    powers = np.square(normalized)
    weights = powers / powers.sum()
    result["power_participation_fraction"] = float(1 / (len(bonds) * np.dot(weights, weights)))
    result["terminal_relative_error"] = abs(summary["current_in"] - summary["current_out"]) / summary["current_in"]
    return result


def compress_input(source, destination, cores_only=False):
    with source.open(encoding="utf-8", newline="") as incoming, gzip.open(destination, "wt", encoding="utf-8", newline="") as outgoing:
        if not cores_only:
            shutil.copyfileobj(incoming, outgoing)
            return
        reader = csv.DictReader(incoming)
        writer = csv.DictWriter(outgoing, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(row for row in reader if row["arm"] == "-1")


def process_case(case, config, generator, root, identity):
    attempt = root / "cases" / case["case_id"] / ("attempt_" + uuid4().hex[:8])
    attempt.mkdir(parents=True)
    sample = attempt / "sample"
    sample.mkdir()
    result = dict(case, attempt_dir=attempt.relative_to(root).as_posix(), status="generation_error")
    started = time.perf_counter()
    with TemporaryDirectory(prefix="absaw_resistance_") as temporary:
        raw = Path(temporary)
        config_text = (f"dim=3\nL={case['L']}\nseed={case['seed']}\ntarget_phi={case['target_phi']:.2f}\n"
                       f"arm_length={config['arm_length']}\na_probability={config['a_probability']}\n")
        (raw / "config.cfg").write_text(config_text, encoding="utf-8")
        (sample / "config.cfg").write_text(config_text, encoding="utf-8")
        target_n = (int(case["L"]**3 * case["target_phi"]) // (1 + 4 * config["arm_length"]))
        attempts = max(config["minimum_attempts"], target_n * config["attempts_per_molecule"])
        filenames = ("candidates", "bonds", "cluster-members", "cluster-summary", "network-summary", "placement-summary")
        command = [str(generator), str(raw / "config.cfg"),
                   *(str(raw / (name + ".csv")) for name in filenames), str(attempts)]
        with (raw / "trajectories.csv").open("wb") as coordinates, (attempt / "generation.log").open("wb") as log:
            completed = subprocess.run(command, stdout=coordinates, stderr=log)
        placement_path = raw / "placement-summary.csv"
        placement = {}
        if placement_path.is_file():
            shutil.copyfile(placement_path, sample / placement_path.name)
            with placement_path.open(encoding="utf-8", newline="") as stream:
                placement, = csv.DictReader(stream)
        result.update(actual_phi=float(placement["actual_phi"]) if placement else None,
                      placement_attempts=int(placement["attempts"]) if placement else None,
                      generation_exit=completed.returncode)
        if completed.returncode == 3 and placement.get("completed") == "0":
            result["status"] = "generation_incomplete"
        elif completed.returncode == 0 and placement.get("completed") == "1":
            result["status"] = "analysis_error"
            atomic_json(raw / "run.json", {"status": "completed", "run_id": case["case_id"], "git_commit": identity["code_commit"]})
            size, cores, bonds, paths, _, _ = engine.read_sample(raw)
            settings_path = root / "settings.json"
            settings = engine.read_settings(settings_path, size)
            summary, nodes, branches = engine.solve_resistance(size, cores, bonds, settings)
            nodes, branches = list(nodes), list(branches)
            result.update(summary)
            result.update(extras(cores, bonds, branches, settings, summary))
            # A core-only trajectory is sufficient to reproduce this molecular resistor graph.
            compress_input(raw / "trajectories.csv", sample / "trajectories.csv.gz", cores_only=True)
            compress_input(raw / "bonds.csv", sample / "bonds.csv.gz")
            atomic_json(sample / "run.json", {"status": "completed", "run_id": case["case_id"],
                        "git_commit": identity["code_commit"], "trajectory_representation": "cores_only"})
            for filename, header, values in (
                ("node-potentials.csv.gz", ("molecule", "electrode", "potential"), nodes),
                ("bond-currents.csv.gz", ("bond_id", "molecule_a", "molecule_b", "current_a_to_b", "dissipated_power"), branches),
            ):
                with gzip.open(attempt / filename, "wt", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(header)
                    writer.writerows(values)
            result["input_sha256"] = {name: engine.filehash(sample / name) for name in ("config.cfg", "trajectories.csv.gz", "bonds.csv.gz")}
    result["seconds"] = time.perf_counter() - started
    return result


def mean_se(values):
    if not values:
        return None, None
    return statistics.mean(values), statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else None


def wilson(successes, total):
    if not total:
        return None, None
    z = statistics.NormalDist().inv_cdf(0.975)
    fraction = successes / total
    denominator = 1 + z*z / total
    center = (fraction + z*z / (2 * total)) / denominator
    half = z * math.sqrt(fraction*(1-fraction)/total + z*z/(4*total*total)) / denominator
    return max(0, center-half), min(1, center+half)


def summarize(rows, config):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["L"], row["target_phi"]].append(row)
    output = []
    for size in sorted(config["sizes"]):
        for percent in sorted(config["phi_percent"]):
            group = grouped[size, percent/100]
            valid = [r for r in group if r["status"] in {"connected", "disconnected", "missing_electrode_contact"}]
            connected = [r for r in valid if r["status"] == "connected"]
            measured = [r for r in valid if r["status"] != "missing_electrode_contact"]
            row = {"L": size, "target_phi": percent/100, "n_planned": len(config["seeds"]),
                   "n_attempted": len(group), "n_valid": len(valid), "n_connected": len(connected),
                   "n_contact_missing": len(valid)-len(measured), "n_failed": len(group)-len(valid),
                   "n_pending": len(config["seeds"])-len(group),
                   "spanning_fraction": len(connected)/len(valid) if valid else None}
            row["spanning_ci_low"], row["spanning_ci_high"] = wilson(len(connected), len(valid))
            # With the specified electrode rule, absent contact carries no terminal current.
            # Keep the engine's undefined raw value, but count that realization as G=0 here.
            row["n_conductance"] = len(measured)
            row["conductance_mean"], row["conductance_se"] = mean_se([
                0.0 if r["status"] == "missing_electrode_contact" else r["effective_conductance"] for r in valid
            ])
            row["conductance_contact_mean"], row["conductance_contact_se"] = mean_se([r["effective_conductance"] for r in measured])
            ratios = [r["resistance_ratio"] for r in connected]
            for name, quantile in (("resistance_q25", 0.25), ("resistance_median", 0.5), ("resistance_q75", 0.75)):
                row[name] = float(np.quantile(ratios, quantile)) if ratios else None
            for metric in ("active_fraction_1e12", "active_fraction_1e10", "active_fraction_1e8", "power_participation_fraction"):
                row[metric+"_mean"], row[metric+"_se"] = mean_se([r[metric] for r in connected])
            row["contacts_high_mean"], _ = mean_se([r["high_count"] for r in valid])
            row["contacts_low_mean"], _ = mean_se([r["low_count"] for r in valid])
            row["terminal_relative_error_max"] = max((r["terminal_relative_error"] for r in connected), default=None)
            output.append(row)
    return output


def export(root, rows, config):
    for filename, values in (("resistance_samples.csv", rows), ("resistance_summary.csv", summarize(rows, config))):
        if not values:
            continue
        temporary = root / (filename + ".tmp")
        fields = list(dict.fromkeys(key for row in values for key in row if key != "input_sha256"))
        with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(values)
        temporary.replace(root / filename)


def run(config, root, generator, retry_failed=False, analyze_only=False):
    root.mkdir(parents=True, exist_ok=True)
    source = Path(__file__)
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source.parent, stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    identity = {"schema_version": 1, "config": config, "code_commit": revision,
                "driver_sha256": engine.filehash(source), "solver_sha256": engine.filehash(Path(engine.__file__)),
                "generator_sha256": engine.filehash(generator), "python": platform.python_version(),
                "numpy": np.__version__, "scipy": scipy.__version__}
    marker = root / "study.json"
    with output_lock(root):
        if marker.exists():
            recorded = json.loads(marker.read_text(encoding="utf-8"))
            for key in ("config", "driver_sha256", "solver_sha256", "generator_sha256"):
                if recorded[key] != identity[key]:
                    raise ValueError(f"Study {key} changed; use the original environment or a new output")
        else:
            if analyze_only:
                raise FileNotFoundError(marker)
            atomic_json(marker, identity)
        atomic_json(root / "settings.json", config["analysis"])
        rows = []
        for case in cases(config):
            destination = root / "cases" / case["case_id"] / "result.json"
            previous = json.loads(destination.read_text(encoding="utf-8")) if destination.is_file() else None
            failure = previous and previous["status"] in {"generation_incomplete", "generation_error", "analysis_error"}
            if previous and not (failure and retry_failed and not analyze_only):
                rows.append(previous)
                continue
            if analyze_only:
                continue
            try:
                result = process_case(case, config, generator, root, identity)
            except Exception as exc:
                # Preserve the failed attempt and make its exclusion explicit in sample counts.
                result = dict(case, status="analysis_error", error=f"{type(exc).__name__}: {exc}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            atomic_json(destination, {k: ("inf" if isinstance(v, float) and v == math.inf else v) for k, v in result.items()})
            # Use the serialized form consistently across initial runs and resumes.
            rows.append(json.loads(destination.read_text(encoding="utf-8")))
            if len(rows) % len(config["seeds"]) == 0:
                export(root, rows, config)
            print(f"{len(rows)} {case['case_id']} {result['status']}", flush=True)
        export(root, rows, config)
        counts = dict(Counter(r["status"] for r in rows))
        print(json.dumps(counts), flush=True)
        return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--generator", required=True, type=Path)
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    print(f"Planned samples: {sum(1 for _ in cases(config))}", flush=True)
    counts = run(config, args.output.resolve(), args.generator.resolve(), args.retry_failed, args.analyze_only)
    if any(counts.get(key, 0) for key in ("generation_incomplete", "generation_error", "analysis_error")):
        parser.exit(1, "Some samples failed. Inspect counts/logs before using the figures.\n")


if __name__ == "__main__":
    main()
