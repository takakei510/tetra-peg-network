"""Scientific definitions and the complete saved-run command, using synthetic runs."""
import csv
import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/analysis'))
import analyze_h2_clusters as app
from h2_cluster.raw import cluster_metrics, count_bonds
from h2_cluster.metrics import h2_metrics
from h2_cluster.stats import summarize


class H2ClusterTests(unittest.TestCase):
    def test_exact_two_excludes_three_and_four(self):
        for k in (2,3,4):
            bonds=[dict(bond_id=i,molecule_a=0,molecule_b=1,endpoint_a=i,endpoint_b=4+i) for i in range(k)]
            h=count_bonds(bonds,2,k)
            self.assertEqual(h[1],int(k==2))
        metrics=h2_metrics(1,6,5,512)
        self.assertEqual(metrics['H2_endpoint_fraction_bonded'],.4)
        self.assertEqual(metrics['H2_endpoint_fraction_all'],1/6)

    def test_duplicate_endpoint_rejected(self):
        rows=[dict(bond_id=i,molecule_a=0,molecule_b=1,endpoint_a=0,endpoint_b=4+i) for i in range(2)]
        with self.assertRaises(ValueError):count_bonds(rows,2,2)

    def test_weighted_ties_and_empty_residual(self):
        rows=[dict(cluster_id=i,molecule_count=s,bond_count=b) for i,s,b in [(9,3,2),(2,3,2),(7,1,0)]]
        m,_=cluster_metrics(rows,7,4,3,3)
        self.assertEqual(m['selected_largest_id'],2)
        self.assertEqual(m['cluster_weighted_excl_one'],2.5)
        self.assertEqual(m['cluster_weighted_excl_all_tied'],1)
        m,_=cluster_metrics([dict(cluster_id=0,molecule_count=2,bond_count=1)],2,1,1,2)
        self.assertIsNone(m['cluster_weighted_excl_one'])
        self.assertIsNone(h2_metrics(0,2,0,512)['H2_endpoint_fraction_bonded'])

    def test_mean_of_ratios_and_valid_sample_count(self):
        rows=[dict(L=8,target_phi=.2,x=.2),dict(L=8,target_phi=.2,x=.8),dict(L=8,target_phi=.2,x=None)]
        r=summarize(rows)[0]
        self.assertEqual(r['x_mean'],.5)
        self.assertEqual(r['x_n'],2)
        self.assertEqual(r['n_samples'],3)

    def test_default_plan(self):
        cfg,cases=app.load_config(app.DEFAULT_CONFIG)
        self.assertEqual(len(cases),105)
        self.assertEqual(len({(c['L'],c['target_phi'],c['seed']) for c in cases}),105)
        self.assertEqual(cases[0]['case_id'],'case00001_L32_seed3_r1')
        self.assertEqual(cases[-1]['case_id'],'case00105_L128_seed202_r1')

    def make_suite(self,root,compressed=False):
        config=dict(schema_version=1,sizes=[8],target_phis=['0.2'],seeds=[3],arm_length=5,a_probability='0.5')
        cfg=root/'plan.json';cfg.write_text(json.dumps(config))
        run=root/'suite/runs/case00001_L8_seed3_r1';run.mkdir(parents=True)
        (run/'config.cfg').write_text('dim=3\nL=8\nseed=3\narm_length=5\na_probability=0.5\ntarget_phi=0.2\n')
        (run/'run.json').write_text(json.dumps(dict(status='completed')))
        tables={'network-summary.csv':[dict(molecules=4,bonds=2,clusters=3,largest=2,largest_fraction=.5)],
            'placement-summary.csv':[dict(completed=1,placed_molecules=4,target_molecules=4,actual_phi=84/512)],
            'cluster-summary.csv':[dict(cluster_id=i,molecule_count=s,bond_count=e) for i,s,e in [(0,2,2),(1,1,0),(2,1,0)]],
            'bonds.csv':[dict(bond_id=i,molecule_a=0,molecule_b=1,endpoint_a=i,endpoint_b=4+i,arm_a=i,arm_b=i,type_a='A',type_b='B',x_a=0,y_a=i,z_a=0,x_b=1,y_b=i,z_b=0,manhattan_distance=1) for i in range(2)]}
        for name,rows in tables.items():
            path=run/(name+'.gz' if compressed else name)
            opener=gzip.open if compressed else open
            with opener(path,'wt',encoding='utf-8',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        return cfg,run

    def test_cli_raw_and_compressed(self):
        for compressed in (False,True):
            with self.subTest(compressed=compressed),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);cfg,run=self.make_suite(root,compressed)
                out=root/'output'
                args=[sys.executable,str(ROOT/'scripts/analysis/analyze_h2_clusters.py'),'--config',str(cfg),'--suite',str(root/'suite'),'--output',str(out),'--no-plots']
                result=subprocess.run(args,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                with (out/'h2_samples.csv').open(encoding='utf-8-sig') as f:row=next(csv.DictReader(f))
                self.assertEqual(int(row['H2']),1)
                self.assertEqual(float(row['cluster_weighted_excl_one']),1)
                self.assertEqual(float(row['H2_endpoint_fraction_bonded']),1)
                self.assertTrue((out/'upload_h2_cluster_results.zip').is_file())
                self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0)

    def test_bad_status_or_config_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg,folder=self.make_suite(root)
            config,cases=app.load_config(cfg)
            (folder/'run.json').write_text(json.dumps(dict(status='incomplete')))
            with self.assertRaises(ValueError):app.analyze_case(folder,cases[0],config)
            (folder/'run.json').write_text(json.dumps(dict(status='completed')))
            config['arm_length']=4
            with self.assertRaises(ValueError):app.analyze_case(folder,cases[0],config)

    def test_missing_and_ambiguous_suites_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg,folder=self.make_suite(root)
            _,cases=app.load_config(cfg)
            with self.assertRaises(ValueError):app.resolve_runs(root,cases)
            first=root/'data/a'/folder.name;first.mkdir(parents=True)
            second=root/'data/b'/folder.name;second.mkdir(parents=True)
            with self.assertRaises(ValueError):app.resolve_runs(root,cases)
            self.assertEqual(app.resolve_runs(root,cases,first.parent),first.parent)


if __name__=='__main__':unittest.main()
