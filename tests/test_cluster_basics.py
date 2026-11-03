"""Verify sample-level weighting, missing Sfin, and input integrity."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

script = Path(__file__).with_name('plot_cluster_basics.py')
if not script.is_file():
    script = Path(__file__).resolve().parents[1]/'scripts/analysis/plot_cluster_basics.py'
spec = importlib.util.spec_from_file_location('cluster_basics', script)
basics = importlib.util.module_from_spec(spec)
spec.loader.exec_module(basics)


def row(N, N1, K, Sall, Sfin, seed):
    return dict(L='8', N=str(N), N1=str(N1), clusters=str(K), bonds='1',
                actual_phi=str(21*N/8**3), target_phi='0.3', seed=str(seed),
                case_id=str(seed), cluster_number_mean=str(N/K),
                cluster_weighted_all=str(Sall), cluster_weighted_excl_one=str(Sfin))


class BasicMetricsTests(unittest.TestCase):
    def test_equal_sample_weight_is_not_a_pooled_ratio(self):
        a = basics.parse_sample(row(4, 2, 3, 1.5, 1, 1))
        b = basics.parse_sample(row(8, 6, 2, 5, 2, 2))
        r = basics.summarize({(8, .3): [a, b]})[0]
        self.assertEqual(r['Pmax_mean'], .625)
        self.assertNotEqual(r['Pmax_mean'], (2+6)/(4+8))
        self.assertEqual(r['Sfin_mean'], 1.5)
        self.assertEqual(r['Sfin_se'], .5)

    def test_undefined_sfin_is_not_replaced_by_zero(self):
        a = basics.parse_sample(row(4, 4, 1, 4, '', 1))
        b = basics.parse_sample(row(4, 2, 3, 1.5, 1, 2))
        r = basics.summarize({(8, .3): [a, b]})[0]
        self.assertEqual(r['Sfin_n'], 1)
        self.assertEqual(r['Sfin_mean'], 1)
        self.assertIsNone(r['Sfin_se'])

    def test_inconsistent_saved_moment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Sfin'):
            basics.parse_sample(row(8, 6, 2, 5, 3, 1))

    def test_duplicate_samples_are_rejected(self):
        a = row(4, 2, 3, 1.5, 1, 1)
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory)/'samples.csv'
            with file.open('w', newline='') as handle:
                w = csv.DictWriter(handle, fieldnames=list(a))
                w.writeheader()
                w.writerows([a, a])
            with self.assertRaisesRegex(ValueError, '重複'):
                basics.load_groups(file, {'sites_per_molecule': 21, 'arms': 4})


if __name__ == '__main__':
    unittest.main()
