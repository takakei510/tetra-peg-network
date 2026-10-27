"""Read-only plot input selection: avoid using the small capacity export by mistake."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('plot_h2',ROOT/'scripts/analysis/plot_h2_validation_results.py')
app=importlib.util.module_from_spec(spec);spec.loader.exec_module(app)


class PlotSelectionTests(unittest.TestCase):
    def make_result(self,root,name,n,status='completed'):
        p=root/name;p.mkdir()
        for f in app.REQUIRED:(p/f).write_text('placeholder')
        (p/'manifest.json').write_text(json.dumps(dict(status=status,completed_samples=n,requested_samples=n)))
        return p

    def test_full_export_wins_over_later_capacity_and_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            full=self.make_result(root,'results_full',18300)
            self.make_result(root,'results_capacity',18)
            incomplete=self.make_result(root,'results_incomplete',20000)
            (incomplete/'EXPORT_INCOMPLETE.json').write_text('{}')
            self.assertEqual(app.resolve_results(root),full)

    def test_partial_export_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.make_result(root,'results_partial',100,'partial')
            with self.assertRaises(ValueError):app.resolve_results(root)

if __name__=='__main__':unittest.main()
