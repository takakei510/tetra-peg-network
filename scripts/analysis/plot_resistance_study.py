#!/usr/bin/env python3
"""Draw the initial observables with explicit sample counts and uncertainty."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def number(row, key):
    return float(row[key]) if row.get(key) else np.nan


def plot(root):
    with (root / "resistance_summary.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    config = json.loads((root / "study.json").read_text(encoding="utf-8"))["config"]
    groups = defaultdict(list)
    for row in rows:
        groups[int(row["L"])].append(row)
    counts = [int(r["n_valid"]) for r in rows]
    eligible = [int(r["n_conductance"]) for r in rows]
    connected = [int(r["n_connected"]) for r in rows]
    planned = len(config["seeds"])
    prefix = "" if counts and min(counts) == planned else "PARTIAL DATA. "
    note = f"{prefix}Planned n={planned}/condition; valid n={min(counts)}-{max(counts)}, contact-valid n={min(eligible)}-{max(eligible)}."
    note_conditional = f"{prefix}Connected samples only; n={min(connected)}-{max(connected)} per condition (planned {planned})."
    output = root / "figures"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

    def save(fig, name):
        fig.savefig(output / (name + ".png"), dpi=200)
        fig.savefig(output / (name + ".svg"))
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for size, group in sorted(groups.items()):
        x = [number(r, "target_phi") for r in group]
        y = np.array([number(r, "spanning_fraction") for r in group])
        interval = np.array([[number(r, "spanning_ci_low") for r in group], [number(r, "spanning_ci_high") for r in group]])
        ax.errorbar(x, y, yerr=np.maximum(0, np.array([y-interval[0], interval[1]-y])), fmt="o", ms=3, capsize=2, label=f"L={size}")
    ax.set(xlabel="Target occupancy", ylabel="Spanning fraction", ylim=(-0.04, 1.04), title="Electrode-to-electrode spanning; 95% Wilson intervals")
    ax.legend()
    fig.supxlabel(note, fontsize=8)
    save(fig, "01_spanning")

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    r0 = config["analysis"]["branch_resistance"]
    for size, group in sorted(groups.items()):
        ax.errorbar([number(r, "target_phi") for r in group], [r0*number(r, "conductance_mean") for r in group],
                    yerr=[r0*number(r, "conductance_se") for r in group], fmt="o", ms=3, capsize=2, label=f"L={size}")
    ax.set(xlabel="Target occupancy", ylabel=r"Mean $G_{eff}R_0$", title="Mean conductance; disconnected and missing-contact G=0; ±1 SE")
    ax.legend()
    fig.supxlabel(note + " Geometric size dependence is included.", fontsize=8)
    save(fig, "02_conductance")

    fig, (ax, sample_ax) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=(3, 1), layout="constrained")
    for size, group in sorted(groups.items()):
        x = [number(r, "target_phi") for r in group]
        y = np.array([number(r, "resistance_median") for r in group])
        low = np.array([number(r, "resistance_q25") for r in group])
        high = np.array([number(r, "resistance_q75") for r in group])
        ax.errorbar(x, y, yerr=np.array([y-low, high-y]), fmt="o", ms=3, capsize=2, label=f"L={size}")
        sample_ax.plot(x, [number(r, "n_connected") for r in group], "o", ms=3)
    ax.set(ylabel=r"$R_{eff}/R_0$", yscale="log", title="Connected resistance; median and interquartile range")
    sample_ax.set(xlabel="Target occupancy", ylabel="Connected n")
    ax.legend()
    fig.supxlabel(note_conditional + " IQR describes sample variation.", fontsize=8)
    save(fig, "03_resistance")

    fig, axes = plt.subplots(1, len(groups), figsize=(5*len(groups), 4.5), sharey=True, squeeze=False, layout="constrained")
    for ax, (size, group) in zip(axes[0], sorted(groups.items())):
        for exponent, marker in ((12, "o"), (10, "s"), (8, "^")):
            ax.errorbar([number(r, "target_phi") for r in group], [number(r, f"active_fraction_1e{exponent}_mean") for r in group],
                        yerr=[number(r, f"active_fraction_1e{exponent}_se") for r in group], fmt=marker, ms=3, capsize=2, label=f"threshold 1e-{exponent}")
        ax.set(xlabel="Target occupancy", title=f"L={size}", ylim=(-0.03, 1.03))
        ax.legend(fontsize=8)
    axes[0, 0].set_ylabel("Fraction of all bonds with current above threshold")
    fig.suptitle(r"Current-carrying bonds; $|I|/(\Delta V/R_0)$ thresholds, ±1 SE")
    fig.supxlabel(note_conditional, fontsize=8)
    save(fig, "04_current_carrying")

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for size, group in sorted(groups.items()):
        ax.errorbar([number(r, "target_phi") for r in group], [number(r, "power_participation_fraction_mean") for r in group],
                    yerr=[number(r, "power_participation_fraction_se") for r in group], fmt="o", ms=3, capsize=2, label=f"L={size}")
    ax.set(xlabel="Target occupancy", ylabel=r"$(\sum P_e)^2/(E\sum P_e^2)$", ylim=(-0.03, 1.03), title="Participation of bonds in power dissipation; ±1 SE")
    ax.legend()
    fig.supxlabel(note_conditional + " Larger values indicate a more distributed load.", fontsize=8)
    save(fig, "05_power_participation")

    # Preserve the count table alongside each set of plots for interpretation.
    (output / "figure_notes.txt").write_text(
        "01: spanning fraction among completed valid samples, including missing-contact samples as not spanning. 95% Wilson interval.\n"
        "02: mean G among all valid completed samples; disconnected and absent-contact realizations contribute G=0 by the chosen electrode convention. Raw missing-contact G stays undefined. Contact-conditioned mean is also exported. ±1 SE across samples.\n"
        "03: median and 25th/75th percentiles of connected R only, with connected counts below. IQR is not uncertainty of the median.\n"
        "04: connected-sample mean fraction of all bonds above three normalized current thresholds. ±1 SE. This is current support under the chosen electrodes.\n"
        "05: connected-sample power participation (sum P)^2/(E*sum P^2). ±1 SE. E includes all bonds.\n"
        "Markers and error bars only; no connecting lines or grid. Sizes share an axis but G includes sample geometry.\n"
        "Do not estimate a threshold or exponent from a smoke run. Consult resistance_summary.csv for counts and failures.\n",
        encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    plot(parser.parse_args().study.resolve())
