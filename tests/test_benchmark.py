"""性能数値を固定値で判定せず、再現性・段階分離・停止・復元を検証する。"""
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from benchmark import load_config, restore
from benchmark_resources import MIB, directory_bytes, disk_guard, measure


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


class BenchmarkTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.config = load_config(ROOT / 'configs/benchmarks/pc-pilot.json')

    def tearDown(self):
        self.temporary.cleanup()

    def launch(self, **changes):
        config = {**self.config, 'sizes': [16], 'target_phis': ['0.1'],
                  'seeds': [3], 'plot': False, **changes}
        cfg = self.folder / 'input.json'
        cfg.write_text(json.dumps(config))
        suite = self.folder / 'suite'
        process = subprocess.run([sys.executable, str(ROOT / 'scripts/benchmark.py'),
                                  '--config', str(cfg), '--output', str(suite)],
                                 cwd=ROOT, capture_output=True, text=True, timeout=60)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        return suite

    def test_reproducibility_stages_and_restore(self):
        suite = self.launch(repeats=2)
        cases = sorted((suite / 'runs').iterdir())
        first, second = [json.loads((p / 'performance.json').read_text()) for p in cases]
        self.assertEqual(first['status'], 'completed')
        self.assertEqual(int(first['placed_molecules']), 19)
        self.assertGreater(first['stages']['c']['peak_rss_bytes'], 0)
        self.assertGreater(first['stages']['postprocess']['peak_rss_bytes'], 0)
        self.assertEqual(first['stages']['plot']['stop_reason'], 'disabled')
        for a, b in zip(first['csv_files'], second['csv_files']):
            self.assertEqual(a['path'], b['path'])
            self.assertEqual(a['sha256'], b['sha256'])
            self.assertFalse((cases[0] / a['path']).exists())
            self.assertEqual(len(gzip.decompress((cases[0] / (a['path'] + '.gz')).read_bytes())), a['raw_bytes'])
        self.assertEqual(first['saved_bytes'], directory_bytes(cases[0]))
        self.assertEqual(json.loads((suite / 'benchmark.json').read_text())['saved_bytes'], directory_bytes(suite))
        summary = read_csv(suite / 'summary.csv')
        self.assertEqual(next(r for r in summary if r['stage'] == 'c')['completed'], '2')
        # 通常の実行入口と、物理CSVおよび派生CSVがバイト単位で一致する。
        baseline = subprocess.run([sys.executable, str(ROOT / 'scripts/run_demo.py'),
                                  '--mode', 'network', '--config', str(cases[0] / 'config.cfg'),
                                  '--no-plot'], cwd=ROOT, capture_output=True, text=True, check=True)
        baseline_run = Path(baseline.stdout.split('保存先: ', 1)[1].strip())
        try:
            for entry in first['csv_files']:
                compressed = cases[0] / (entry['path'] + '.gz')
                self.assertEqual(gzip.decompress(compressed.read_bytes()),
                                 (baseline_run / entry['path']).read_bytes())
        finally:
            shutil.rmtree(baseline_run)
        restore(cases[0])
        restore(cases[0])  # 同じハッシュなら既存を変更せず再照合。
        target = cases[0] / 'trajectories.csv'
        target.write_text('do not overwrite')
        with self.assertRaisesRegex(ValueError, '既存ファイル'):
            restore(cases[0])
        self.assertEqual(target.read_text(), 'do not overwrite')

    def test_incomplete_is_not_success(self):
        suite = self.launch(sizes=[2], target_phis=['0.625'], arm_length=1,
                            min_attempts=5, attempts_per_molecule=1)
        case = json.loads(next((suite / 'runs').glob('*/performance.json')).read_text())
        self.assertEqual(case['status'], 'incomplete')
        self.assertEqual(case['stages']['c']['status'], 'incomplete')
        self.assertEqual(case['stages']['postprocess']['stop_reason'], 'generation_not_completed')
        summary = next(r for r in read_csv(suite / 'summary.csv') if r['stage'] == 'c')
        self.assertEqual(summary['completed'], '0')
        self.assertEqual(summary['wall_seconds_median'], '')

    def test_plot_stage_has_separate_metrics(self):
        suite = self.launch(plot=True, plot_max_molecules=20)
        case = json.loads(next((suite / 'runs').glob('*/performance.json')).read_text())
        plot = case['stages']['plot']
        self.assertEqual(plot['status'], 'completed')
        self.assertGreater(plot['render_save_seconds'], 0)
        self.assertGreater(plot['read_validate_seconds'], 0)
        self.assertGreaterEqual(plot['wall_seconds'], plot['render_save_seconds'])
        self.assertEqual(len(list((suite / case['directory'] / 'plots').glob('*.png'))), 4)

    def test_timeout_records_partial(self):
        result = measure([sys.executable, '-c', 'import time;time.sleep(5)'],
                         self.folder / 'out', self.folder / 'err',
                         {**self.config, 'stage_timeout_seconds': .06}, self.folder, self.folder)
        self.assertEqual(result['status'], 'limited')
        self.assertEqual(result['stop_reason'], 'time_limit')
        self.assertLess(result['wall_seconds'], 2)
        self.assertLess(result['exit_code'], 0)

    def test_rss_limit(self):
        result = measure([sys.executable, '-c', 'import time;x=bytearray(96*1024**2);time.sleep(2)'],
                         self.folder / 'out', self.folder / 'err',
                         {**self.config, 'max_rss_mib': 48}, self.folder, self.folder)
        self.assertEqual(result['status'], 'limited')
        self.assertIn(result['stop_reason'], ('rss_limit', 'rss_limit_after_exit'))

    def test_disk_guard_before_and_during_execution(self):
        config = {**self.config, 'max_case_output_mib': 1}
        result = measure([sys.executable, '-c', 'import time;open("'+str(self.folder / 'large')+'","wb").write(b"x"*2*1024**2);time.sleep(2)'],
                         self.folder / 'out', self.folder / 'err', config, self.folder, self.folder)
        self.assertEqual(result['stop_reason'], 'case_output_limit')
        self.assertEqual(disk_guard(config, self.folder, self.folder), 'case_output_limit')

    def test_bad_configuration(self):
        path = self.folder / 'bad.json'
        path.write_text(json.dumps({**self.config, 'unknown': 1}))
        with self.assertRaises(ValueError):
            load_config(path)
        path.write_text('{"schema_version":1,"schema_version":1}')
        with self.assertRaisesRegex(ValueError, '重複'):
            load_config(path)


if __name__ == '__main__':
    unittest.main()
