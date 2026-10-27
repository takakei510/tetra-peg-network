import json
from pathlib import Path
import sys
import tempfile
import unittest
import subprocess
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'scripts/analysis'))
import run_h2_validation as app
from cluster_h2_relation import graph_rows, bin_rows


def bonds(pairs):
    return [dict(molecule_a=a,molecule_b=b) for a,b in pairs]


class RelationTests(unittest.TestCase):
    def test_known_components_h2_h3_isolates(self):
        rows,owner=graph_rows(7,bonds([(0,1),(0,1),(1,2),(3,4),(3,4),(3,4)]))
        self.assertEqual([(r['size'],r['bonds'],r['h2']) for r in rows],[(3,3,1),(2,3,0),(1,0,0),(1,0,0)])
        self.assertEqual(owner,[0,0,0,1,1,2,3])
        self.assertEqual(rows[0]['endpoint_fraction_all'],1/3)
        self.assertIsNone(rows[2]['endpoint_fraction_bonded'])

    def test_zero_h2_clusters_included_and_tie_rule(self):
        rows,_=graph_rows(4,bonds([(0,1),(0,1),(2,3)]))
        bs=bin_rows(rows)
        all_bin=next(b for b in bs if b['scope']=='all')
        finite=next(b for b in bs if b['scope']!='all')
        self.assertEqual(all_bin['fraction_mean'],.25)
        self.assertEqual(finite['fraction_mean'],0)

    def test_se_not_sd(self):
        self.assertAlmostEqual(app.stats([1,3])['se'],1)
        self.assertIsNone(app.stats([1])['se'])

    def test_stable_extension_unique_seeds(self):
        config=app.read(app.DEFAULT)
        a=app.jobs(config,5);b=app.jobs(config,10)
        lookup={r['case_id']:r['seed'] for r in b}
        self.assertTrue(all(lookup[r['case_id']]==r['seed'] for r in a))
        self.assertEqual(len({r['seed'] for r in b}),len(b))

    def test_real_parallel_resume_extension_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);cfg=app.read(app.DEFAULT)
            cfg['sizes']=[8];cfg['target_phis']=['0.30','0.40']
            cfg['benchmark'].update(max_rss_mib=128,min_free_disk_mib=32,max_case_output_mib=16)
            path=tmp/'cfg.json';path.write_text(json.dumps(cfg));out=tmp/'study'
            command=[sys.executable,str(ROOT/'scripts/run_h2_validation.py'),'run','--config',str(path),'--output',str(out),'--workers','2','--max-storage-gib','1','--no-plots']
            subprocess.run(command+['--samples','2'],check=True)
            first=next((out/'cases').glob('*/result.json'));before=first.stat().st_mtime_ns
            subprocess.run(command+['--samples','3'],check=True)
            self.assertEqual(first.stat().st_mtime_ns,before)
            self.assertEqual(len(list((out/'cases').glob('*/result.json'))),6)
            result=app.read(first)
            cluster=first.parent/'clusters.csv.gz';cluster.write_bytes(cluster.read_bytes()+b'bad')
            failed=subprocess.run(command+['--samples','3'],capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            self.assertIn('hash mismatch',failed.stderr)

if __name__=='__main__':unittest.main()
