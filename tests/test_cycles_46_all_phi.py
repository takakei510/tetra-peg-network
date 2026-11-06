import csv
import gzip
import importlib.util
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
import analyze_cycles_46_all_phi as m

class Tests(unittest.TestCase):
    def test_mean_of_trial_ratios_not_pooled_and_zero_included(self):
        c=dict(L=128,phi_values=[.43],samples_per_phi=2)
        rows=[]
        for n,c4 in [(4,1),(12,0)]:
            r=dict(target_phi=.43,actual_phi=.429,N=n,H2=0,C4_mol=c4,C6_mol=0,C4_per_100_molecules=100*c4/n,C6_per_100_molecules=0)
            rows.append(r)
        s=m.summarize(rows,c)[0]
        self.assertEqual(s['C4_mol_mean'],.5)
        self.assertEqual(s['C4_per_100_molecules_mean'],12.5)
        self.assertEqual(s['C4_per_100_molecules_se'],12.5)
        self.assertNotEqual(s['C4_per_100_molecules_mean'],100/16)
        with self.assertRaises(ValueError):m.summarize(rows[:1],c)

    def test_engine_hash_guard_and_bad_counts(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.py';p.write_text('')
            with self.assertRaises(ValueError):m.load_engine(p,m.PINNED_ENGINE)
        bad={k:0 for k in m.COUNT_KEYS};bad['C4_mol']=1
        with self.assertRaises(ValueError):m.validate_totals(bad)

    def test_multi_phi_resume_legacy_cache_and_changed_raw(self):
        engine_path=Path(m.__file__).with_name('analyze_cycles_46.py')
        engine=m.load_engine(engine_path,m.PINNED_ENGINE)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);samples=[]
            # Two disconnected squares: C4=2, C6=0. No synthetic result is delivered as research evidence.
            for phi in [.10,.43,.70]:
                for trial in range(2):
                    case=f'L128_phi{phi:.2f}_i{trial:06d}';folder=root/'cases'/case/'runs'/case;folder.mkdir(parents=True)
                    edges=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4)]
                    used={};bb=[]
                    for bid,(a,b) in enumerate(edges):
                        ea=4*a+used.get(a,0);eb=4*b+used.get(b,0);used[a]=used.get(a,0)+1;used[b]=used.get(b,0)+1
                        bb.append(dict(bond_id=bid,molecule_a=a,molecule_b=b,endpoint_a=ea,endpoint_b=eb))
                    clusters,owner,degree=engine.count_cycles(8,bb)
                    members=[dict(molecule=k,cluster_id=owner[k],cluster_size=clusters[owner[k]]['size'],bond_degree=degree[k]) for k in range(8)]
                    for name,rr in [('bonds',bb),('cluster-members',members)]:
                        with gzip.open(folder/(name+'.csv.gz'),'wt',newline='') as f:
                            w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
                    samples.append(dict(case_id=case,L=128,target_phi=phi,seed=trial,N=8,N1=4,bonds=8,clusters=2,H2=0))
            source=root/'samples.csv';m.write_csv(source,samples)
            c=dict(schema_version=1,L=128,phi_values=[.10,.43,.70],samples_per_phi=2,pilot_phi_values=[.10,.43,.70],engine_sha256=m.PINNED_ENGINE,scope='all')
            cfg=root/'config.json';cfg.write_text(json.dumps(c))
            _,groups=m.select_samples(root,c,source)
            for sample in groups[.43]:
                rec=engine.inspect_case(root,sample)
                m.atomic_json(root/'cycle46_cache'/m.PINNED_ENGINE[:16]/(sample['case_id']+'.json'),rec)
            cmd=[sys.executable,m.__file__,'--study',str(root),'--samples-csv',str(source),'--config',str(cfg)]
            pilot=subprocess.run(cmd+['--mode','pilot'],check=True,capture_output=True,text=True)
            self.assertIn('PILOT ONLY',pilot.stdout)
            for index in range(2):
                out=root/f'out{index}'
                subprocess.run(cmd+['--mode','run','--output',str(out)],check=True,capture_output=True,text=True)
                meta=json.loads((out/'run_metadata.json').read_text())
                self.assertEqual(meta['samples'],6)
                if index==1:self.assertEqual(meta['cache_sources'],{'compact-cache':6})
                with (out/'phi_summary.csv').open(encoding='utf-8-sig') as f:
                    for r in csv.DictReader(f):
                        self.assertEqual(float(r['C4_mol_mean']),2)
                        self.assertEqual(float(r['C6_mol_mean']),0)
                        self.assertEqual(float(r['C4_per_100_molecules_mean']),25)
                from zipfile import ZipFile
                with ZipFile(out/'cycle46_all_phi_share.zip') as z:
                    self.assertIsNone(z.testzip());self.assertEqual(len(z.namelist()),8)
            sample=groups[.10][0];folder=root/'cases'/sample['case_id']/'runs'/sample['case_id'];bp=folder/'bonds.csv.gz'
            raw=gzip.decompress(bp.read_bytes());bp.write_bytes(gzip.compress(raw,mtime=1234))
            with patch.object(engine,'inspect_case',wraps=engine.inspect_case) as inspected:
                _,_,origin=m.process_case(engine,root,sample,m.PINNED_ENGINE)
                self.assertEqual(origin,'new-dfs');self.assertEqual(inspected.call_count,1)
            # A repeated seed in the same phi must be rejected.
            broken=list(samples);broken[1]={**broken[1],'seed':0};m.write_csv(source,broken)
            with self.assertRaises(ValueError):m.select_samples(root,c,source)

if __name__=='__main__':unittest.main()
