#!/usr/bin/env python3
"""Analyze one completed AB-SAW sample as a DC resistor network."""

import argparse
import csv
import json
import math
import subprocess
import sys
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import scipy
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from h2_cluster.raw import filehash, locate


DEFAULT_SETTINGS = {
    "axis": "x",
    "electrode_layers": 1,
    "branch_resistance": 1.0,
    "voltage_high": 1.0,
    "voltage_low": 0.0,
}
INPUT_NAMES = ("config.cfg", "trajectories.csv", "bonds.csv")


def read_settings(path, box_size):
    settings = DEFAULT_SETTINGS.copy()
    if path is not None:
        supplied = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(supplied, dict):
            raise ValueError("Analysis settings must be a JSON object.")
        unknown = supplied.keys() - settings.keys()
        if unknown:
            raise ValueError("Unknown analysis setting: " + ", ".join(sorted(unknown)))
        settings.update(supplied)

    axis = settings["axis"]
    layers = settings["electrode_layers"]
    if axis not in ("x", "y", "z"):
        raise ValueError("axis must be x, y, or z.")
    if type(layers) is not int or not 1 <= layers <= box_size // 2:
        raise ValueError("electrode_layers must satisfy 1 <= electrode_layers <= L/2.")
    for key in ("branch_resistance", "voltage_high", "voltage_low"):
        value = settings[key]
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError(f"{key} must be a finite number.")
        settings[key] = float(value)
    delta_v = settings["voltage_high"] - settings["voltage_low"]
    if settings["branch_resistance"] <= 0 or not 0 < delta_v < math.inf:
        raise ValueError("branch_resistance and voltage_high - voltage_low must be positive and finite.")
    return settings


def csv_rows(path, required):
    import gzip

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames or []
        missing = set(required) - set(header)
        if missing:
            raise ValueError(f"{path.name}: missing columns {', '.join(sorted(missing))}.")
        if len(header) != len(set(header)):
            raise ValueError(f"{path.name}: duplicate column names.")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path.name}:{reader.line_num}: wrong number of columns.")
            try:
                yield tuple(int(row[key]) for key in required)
            except ValueError as exc:
                raise ValueError(f"{path.name}:{reader.line_num}: expected integer values.") from exc


def read_sample(folder):
    manifest_path = folder / "run.json"
    manifest = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        if not isinstance(manifest, dict) or manifest.get("status") != "completed":
            raise ValueError("run.json does not record a completed run.")

    config_path = folder / INPUT_NAMES[0]
    config = {}
    for line_number, line in enumerate(config_path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.partition("#")[0].strip()
        if not line:
            continue
        key, sep, value = (part.strip() for part in line.partition("="))
        if not sep or not key or not value or key in config:
            raise ValueError(f"config.cfg:{line_number}: invalid or duplicate setting.")
        config[key] = value
    if "L" not in config or int(config["L"]) < 2 or int(config.get("dim", "3")) != 3:
        raise ValueError("config.cfg must describe a three-dimensional box with L >= 2.")
    box_size = int(config["L"])

    paths = [config_path, *(locate(folder, name) for name in INPUT_NAMES[1:])]
    cores = {}
    seen_coordinates = set()
    for molecule, arm, x, y, z in csv_rows(paths[1], ("molecule", "arm", "x", "y", "z")):
        if arm != -1:
            continue
        coordinates = (x, y, z)
        if molecule < 0 or molecule in cores or coordinates in seen_coordinates:
            raise ValueError(f"{paths[1].name}: invalid molecule ID or duplicate core for molecule {molecule}.")
        if any(not 0 <= q < box_size for q in coordinates):
            raise ValueError(f"{paths[1].name}: core {molecule} lies outside the box.")
        cores[molecule] = coordinates
        seen_coordinates.add(coordinates)
    if not cores:
        raise ValueError("The trajectories file has no core rows (arm=-1).")

    bonds = []
    seen_bonds = set()
    for bond_id, a, b in csv_rows(paths[2], ("bond_id", "molecule_a", "molecule_b")):
        if bond_id < 0 or bond_id in seen_bonds:
            raise ValueError(f"{paths[2].name}: negative or duplicate bond ID {bond_id}.")
        if a == b or a not in cores or b not in cores:
            raise ValueError(f"{paths[2].name}: self-bond or unknown molecule in bond {bond_id}.")
        bonds.append((bond_id, a, b))
        seen_bonds.add(bond_id)
    return box_size, cores, bonds, paths, config, manifest


def solve_resistance(box_size, cores, bonds, settings):
    molecule_ids = list(cores)
    index_by_molecule = {molecule: index for index, molecule in enumerate(molecule_ids)}
    node_count = len(cores)
    branch_a = np.array([index_by_molecule[a] for _, a, _ in bonds], dtype=np.intp)
    branch_b = np.array([index_by_molecule[b] for _, _, b in bonds], dtype=np.intp)
    axis = settings["axis"]
    layers = settings["electrode_layers"]
    axis_index = ("x", "y", "z").index(axis)
    coordinates = np.array([cores[molecule][axis_index] for molecule in molecule_ids])
    high_mask = coordinates < layers
    low_mask = coordinates >= box_size - layers

    # Sum parallel branches in the adjacency matrix while keeping each bond for output.
    adjacency = coo_matrix(
        (np.ones(2 * len(bonds)), (np.r_[branch_a, branch_b], np.r_[branch_b, branch_a])),
        shape=(node_count, node_count),
    ).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    high_components = np.unique(labels[high_mask])
    low_components = np.unique(labels[low_mask])
    spanning = np.intersect1d(high_components, low_components)
    free_indices = np.flatnonzero(np.isin(labels, spanning) & ~(high_mask | low_mask))

    # Solve with unit voltage and unit branch conductance, then rescale the observables.
    unit_potential = np.full(node_count, np.nan)
    unit_potential[np.isin(labels, low_components)] = 0.0
    unit_potential[np.isin(labels, high_components)] = 1.0
    unit_potential[low_mask] = 0.0
    if free_indices.size:
        laplacian = diags(np.asarray(adjacency.sum(axis=1)).ravel()) - adjacency
        rhs = np.asarray(adjacency[free_indices][:, high_mask].sum(axis=1)).ravel()
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            unit_potential[free_indices] = spsolve(laplacian[free_indices][:, free_indices], rhs)
        solved = unit_potential[free_indices]
        if not np.all(np.isfinite(solved)) or np.any((solved < -1e-9) | (solved > 1 + 1e-9)):
            raise ValueError("The sparse solve produced invalid internal potentials.")

    unit_currents = unit_potential[branch_a] - unit_potential[branch_b]
    unit_currents[~np.isin(labels[branch_a], spanning)] = 0.0
    injection = (
        np.bincount(branch_a, weights=unit_currents, minlength=node_count)
        - np.bincount(branch_b, weights=unit_currents, minlength=node_count)
    )
    unit_in = float(injection[high_mask].sum())
    unit_out = float(-injection[low_mask].sum())
    unit_power = float(unit_currents @ unit_currents)
    unit_kcl = float(np.max(np.abs(injection[free_indices]), initial=0.0))
    if spanning.size:
        if unit_in <= 0 or unit_out <= 0:
            raise ValueError("A spanning network produced nonpositive terminal current.")
        if (unit_kcl > 1e-10 + 1e-8 * np.max(np.abs(unit_currents), initial=0.0)
                or not math.isclose(unit_in, unit_out, rel_tol=1e-8, abs_tol=1e-10)
                or not math.isclose(unit_power, unit_in, rel_tol=1e-8, abs_tol=1e-10)):
            raise ValueError("Current or dissipated power did not balance.")

    delta_v = settings["voltage_high"] - settings["voltage_low"]
    branch_resistance = settings["branch_resistance"]
    current_scale = delta_v / branch_resistance
    power_scale = delta_v * current_scale
    if not (0 < current_scale < math.inf and 0 < power_scale < math.inf):
        raise ValueError("Voltage and resistance scales exceed floating-point range.")
    potentials = settings["voltage_low"] + delta_v * unit_potential
    potentials[unit_potential == 1.0] = settings["voltage_high"]
    currents = unit_currents * current_scale
    powers = np.square(unit_currents) * power_scale

    if spanning.size:
        status = "connected"
        conductance = unit_in / branch_resistance
        resistance_ratio = 1.0 / unit_in
        resistance = branch_resistance * resistance_ratio
        if not all(0 < value < math.inf for value in (conductance, resistance, unit_in * current_scale, unit_power * power_scale)):
            raise ValueError("The terminal observables exceed floating-point range.")
    elif high_mask.any() and low_mask.any():
        status, conductance, resistance, resistance_ratio = "disconnected", 0.0, math.inf, math.inf
    else:
        status, conductance, resistance, resistance_ratio = "missing_electrode_contact", None, None, None

    summary = {
        "status": status, "axis": axis, "electrode_layers": layers, "box_size": box_size,
        "node_count": node_count, "bond_count": len(bonds),
        "high_count": int(high_mask.sum()), "low_count": int(low_mask.sum()),
        "spanning_component_count": int(spanning.size),
        "current_in": unit_in * current_scale, "current_out": unit_out * current_scale,
        "effective_conductance": conductance, "effective_resistance": resistance,
        "resistance_ratio": resistance_ratio, "kcl_residual_max": unit_kcl * current_scale,
        "current_balance_error": abs(unit_in - unit_out) * current_scale,
        "power_balance_error": abs(unit_power - unit_in) * power_scale,
    }
    node_rows = (
        (molecule, "high" if high_mask[i] else "low" if low_mask[i] else "",
         "" if np.isnan(potential) else float(potential))
        for i, (molecule, potential) in enumerate(zip(molecule_ids, potentials))
    )
    bond_rows = (
        (*bond, float(current), float(power))
        for bond, current, power in zip(bonds, currents, powers)
    )
    return summary, node_rows, bond_rows


def write_csv(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def write_results(output, summary, nodes, bonds, metadata):
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".resistance-", dir=output.parent) as scratch:
        staged = Path(scratch)
        write_csv(staged / "resistance-summary.csv", summary.keys(), [summary.values()])
        write_csv(staged / "node-potentials.csv", ("molecule", "electrode", "potential"), nodes)
        write_csv(staged / "bond-currents.csv",
                  ("bond_id", "molecule_a", "molecule_b", "current_a_to_b", "dissipated_power"), bonds)
        marker = "analysis-run.json"
        (staged / marker).write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        # The exclusive mkdir prevents concurrent writers from using the same destination.
        output.mkdir()
        moved = []
        try:
            for name in ("resistance-summary.csv", "node-potentials.csv", "bond-currents.csv", marker):
                (staged / name).rename(output / name)
                moved.append(output / name)
        except BaseException:
            for file in moved:
                file.unlink()
            output.rmdir()
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_dir", type=Path, help="Completed run with config.cfg and both CSVs.")
    parser.add_argument("--output-dir", type=Path, required=True, help="New result directory.")
    parser.add_argument("--settings", type=Path, help="Optional JSON containing electrode and resistance settings.")
    args = parser.parse_args()
    try:
        folder = args.sample_dir.resolve()
        output = args.output_dir.resolve()
        if output.exists():
            raise FileExistsError(f"Output directory already exists: {output}")
        box_size, cores, bonds, paths, config, manifest = read_sample(folder)
        settings = read_settings(args.settings, box_size)
        summary, node_rows, bond_rows = solve_resistance(box_size, cores, bonds, settings)
        script = Path(__file__).resolve()
        try:
            revision = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=script.parent, stderr=subprocess.DEVNULL, text=True,
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            revision = None
        metadata = {
            "schema_version": 1, "sample_id": manifest.get("run_id", folder.name),
            "sample_dir": str(folder),
            "inputs": {name: {"path": str(path), "sha256": filehash(path)} for name, path in zip(INPUT_NAMES, paths)},
            "generation_config": config, "generation_commit": manifest.get("git_commit"),
            "analysis_settings": settings, "analysis_commit": revision, "analysis_script_sha256": filehash(script),
            "solver": "scipy.sparse.linalg.spsolve", "numpy_version": np.__version__, "scipy_version": scipy.__version__,
        }
        write_results(output, summary, node_rows, bond_rows, metadata)
    except (OSError, ValueError, ArithmeticError, RuntimeError, MatrixRankWarning) as exc:
        parser.exit(1, f"Resistance analysis failed: {exc}\n")

    ratio = summary["resistance_ratio"]
    print(f"status={summary['status']} current_in={summary['current_in']:.12g} "
          f"R_eff/R_0={ratio if ratio is not None else 'undefined'}")
    print(f"Results: {output}")


if __name__ == "__main__":
    main()
