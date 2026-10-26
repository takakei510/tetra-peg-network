"""Detail study planning and a real small C-to-CSV-to-analysis run."""
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_phi_detail as app


class DetailTests(unittest.TestCase):
    def test_full_grid_no_duplicates_or_old_endpoints(self):
        plan,groups=app.load_plan(app.DEFAULT_CONFIG)
        self.assertEqual([app.sample_count(cfg) for _,cfg in groups],[20,45,25,45])
        conditions={(L,float(phi),seed) for _,cfg in groups for L in cfg['sizes'] for phi in cfg['target_phis'] for seed in cfg['seeds']}
        expected={(L,k/100,seed) for L in [32,64,128] for k in range(41,50) for seed in [3,17,42,101,202]}
        self.assertEqual(conditions,expected)
        self.assertEqual(sum(app.sample_count(cfg) for _,cfg in groups),len(conditions))

    def test_duplicate_conditions_rejected(self):
        plan=json.loads(app.DEFAULT_CONFIG.read_text())
        duplicate=dict(plan['groups'][0],name='duplicate')
        plan['groups'].append(duplicate)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.json';path.write_text(json.dumps(plan))
            with self.assertRaises(ValueError):app.load_plan(path)

    def test_plan_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'not_created'
            app.run('plan',app.DEFAULT_CONFIG,out)
            self.assertFalse(out.exists())

    def test_partial_batch_not_silently_rerun(self):
        _,groups=app.load_plan(app.DEFAULT_CONFIG)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/groups[0][0]/'simulation').mkdir(parents=True)
            with patch.object(app.subprocess,'run') as invoke:
                with self.assertRaises(ValueError):app.run_group(root,*groups[0])
                invoke.assert_not_called()

    def test_real_pilot_then_full_reuses_completed_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);output=root/'study';config=root/'plan.json'
            plan=json.loads(app.DEFAULT_CONFIG.read_text())
            plan['benchmark']['seeds']=[3]
            plan['benchmark']['min_free_disk_mib']=64
            plan['groups']=[dict(name='pilot',sizes=[8],target_phis=['0.42']),dict(name='remaining',sizes=[8],target_phis=['0.46'])]
            config.write_text(json.dumps(plan))
            first=app.run('pilot',config,output)
            self.assertTrue(first.is_file())
            marker=output/'pilot/simulation/benchmark.json'
            before=marker.read_bytes();mtime=marker.stat().st_mtime_ns
            full=app.run('full',config,output)
            self.assertEqual(marker.read_bytes(),before)
            self.assertEqual(marker.stat().st_mtime_ns,mtime)
            self.assertTrue(full.is_file())
            with (full.parent/'h2_samples.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),2)
            self.assertEqual({r['case_id'].split('/')[0] for r in rows},{'pilot','remaining'})
            self.assertEqual({r['target_phi'] for r in rows},{'0.42','0.46'})
            manifest=json.loads((full.parent/'run_manifest.json').read_text())
            self.assertEqual(manifest['samples'],2)
            self.assertFalse((output/'running.lock').exists())


if __name__=='__main__':unittest.main()
