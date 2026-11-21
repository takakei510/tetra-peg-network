"""Regression checks for repeated-sample statistics and compressed graph inputs."""

import copy
import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
import resistance_study as study


class StudyTest(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / "configs/analysis/resistance-study-smoke.json").read_text())
        self.config.update(sizes=[32], phi_percent=[60], seeds=[3, 4, 5, 6])

    def test_unconditional_average_and_counts(self):
        common = {"L": 32, "target_phi": 0.6, "high_count": 2, "low_count": 2,
                  "active_fraction_1e12": 0.5, "active_fraction_1e10": 0.5,
                  "active_fraction_1e8": 0.5, "power_participation_fraction": 0.25,
                  "terminal_relative_error": 0.0}
        connected = dict(common, seed=3, status="connected", effective_conductance=2.0, resistance_ratio=0.5)
        disconnected = dict(common, seed=4, status="disconnected", effective_conductance=0.0, resistance_ratio="inf")
        absent = dict(common, seed=5, status="missing_electrode_contact", effective_conductance=None, resistance_ratio=None)
        failed = dict(common, seed=6, status="generation_incomplete")
        row, = study.summarize([connected, disconnected, absent, failed], self.config)
        self.assertEqual((row["n_valid"], row["n_connected"], row["n_contact_missing"], row["n_failed"]), (3, 1, 1, 1))
        self.assertAlmostEqual(row["conductance_mean"], 2 / 3)
        self.assertAlmostEqual(row["conductance_se"], 2 / 3)
        self.assertEqual(row["conductance_contact_mean"], 1)
        self.assertEqual(row["resistance_median"], 0.5)
        self.assertLess(row["spanning_ci_low"], 1 / 3)
        self.assertGreater(row["spanning_ci_high"], 1 / 3)

    def test_boundary_wilson_and_missing_conditions(self):
        low, high = study.wilson(0, 100)
        self.assertAlmostEqual(low, 0, places=14)
        self.assertGreater(high, 0)
        low, high = study.wilson(100, 100)
        self.assertLess(low, 1)
        self.assertAlmostEqual(high, 1, places=14)
        row, = study.summarize([], self.config)
        self.assertEqual(row["n_pending"], 4)
        self.assertIsNone(row["conductance_mean"])
        self.assertIsNone(row["resistance_median"])


if __name__ == "__main__":
    unittest.main()
