"""Analytic and CLI checks for the AB-SAW resistance postprocessor."""

import csv
import gzip
import json
import math
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "analyze_resistance.py"
sys.path.insert(0, str(SCRIPT.parent))
from analyze_resistance import DEFAULT_SETTINGS, solve_resistance


class ResistanceTest(unittest.TestCase):
    def test_analytic_networks(self):
        cases = [
            ("single", [0, 3], [(0, 1)], 1),
            ("series", [0, 1, 3], [(0, 1), (1, 2)], 2),
            ("parallel", [0, 3], [(0, 1), (0, 1)], 0.5),
            ("cycle", [0, 1, 3, 2], [(0, 1), (1, 2), (2, 3), (3, 0)], 1),
            ("dangling", [0, 1, 3, 2], [(0, 1), (1, 2), (1, 3)], 2),
            ("two spanning components", [0, 1, 3, 0, 1, 3], [(0, 1), (1, 2), (3, 4), (4, 5)], 1),
        ]
        for name, xs, endpoints, expected in cases:
            with self.subTest(name=name):
                cores = {i: (x, i % 4, i // 4) for i, x in enumerate(xs)}
                bonds = [(i, a, b) for i, (a, b) in enumerate(endpoints)]
                summary, nodes, branches = solve_resistance(4, cores, bonds, DEFAULT_SETTINGS)
                self.assertEqual(summary["status"], "connected")
                self.assertAlmostEqual(summary["resistance_ratio"], expected, places=12)
                self.assertAlmostEqual(summary["current_in"], 1 / expected, places=12)
                for field in ("kcl_residual_max", "current_balance_error", "power_balance_error"):
                    self.assertAlmostEqual(summary[field], 0, places=12)
                if name == "parallel":
                    self.assertEqual([row[3] for row in branches], [1, 1])
                if name == "dangling":
                    self.assertAlmostEqual(list(branches)[-1][3], 0, places=12)
                    self.assertAlmostEqual(list(nodes)[-1][2], 0.5, places=12)

    def test_one_sided_and_floating_components(self):
        cores = {7: (0, 0, 0), 9: (1, 0, 0), 20: (3, 0, 0), 30: (1, 1, 0), 32: (2, 1, 0)}
        bonds = [(10, 7, 9), (11, 30, 32)]
        summary, nodes, branches = solve_resistance(4, cores, bonds, DEFAULT_SETTINGS)
        self.assertEqual(summary["status"], "disconnected")
        self.assertEqual(summary["current_in"], 0)
        self.assertTrue(math.isinf(summary["resistance_ratio"]))
        self.assertEqual([row[2] for row in nodes], [1, 1, 0, "", ""])
        self.assertEqual([row[3] for row in branches], [0, 0])

        del cores[20]
        summary, _, _ = solve_resistance(4, cores, bonds, DEFAULT_SETTINGS)
        self.assertEqual(summary["status"], "missing_electrode_contact")
        self.assertIsNone(summary["resistance_ratio"])

    def test_cli_reads_gzip_and_keeps_outputs(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            sample, output = folder / "sample", folder / "analysis"
            sample.mkdir()
            (sample / "config.cfg").write_text("dim=3\nL=5\nseed=9\n", encoding="utf-8")
            (sample / "trajectories.csv").write_text(
                "molecule,arm,x,y,z\n3,-1,0,1,0\n5,-1,0,2,0\n8,-1,0,4,0\n", encoding="utf-8",
            )
            with gzip.open(sample / "bonds.csv.gz", "wt", encoding="utf-8") as stream:
                stream.write("bond_id,molecule_a,molecule_b\n1,3,5\n2,8,5\n")
            (sample / "run.json").write_text('{"status":"completed","run_id":"test"}', encoding="utf-8")
            settings = folder / "settings.json"
            settings.write_text(json.dumps({
                "axis": "y", "electrode_layers": 2, "branch_resistance": 2,
                "voltage_high": 5, "voltage_low": 1,
            }), encoding="utf-8")
            command = [sys.executable, str(SCRIPT), str(sample), "--output-dir", str(output), "--settings", str(settings)]
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with (output / "resistance-summary.csv").open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["status"], "connected")
            self.assertAlmostEqual(float(row["current_in"]), 1)
            self.assertAlmostEqual(float(row["effective_resistance"]), 4)
            self.assertAlmostEqual(float(row["resistance_ratio"]), 2)
            with (output / "bond-currents.csv").open(newline="", encoding="utf-8") as stream:
                currents = [float(row["current_a_to_b"]) for row in csv.DictReader(stream)]
            self.assertEqual(currents, [1, -1])
            metadata = json.loads((output / "analysis-run.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["sample_id"], "test")
            self.assertEqual(len(metadata["inputs"]["bonds.csv"]["sha256"]), 64)
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertTrue((output / "analysis-run.json").is_file())


if __name__ == "__main__":
    unittest.main()
