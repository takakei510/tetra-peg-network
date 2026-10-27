"""Storage limits, scientific equivalence and budgeted export regression tests."""
import csv
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from zipfile import ZipFile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'scripts/analysis'))
import benchmark
import run_h2_validation as app
from storage_budget import WriteBudget, StorageLimit
from benchmark_resources import directory_bytes


class StorageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-B','build/demo_network','build/check_network_config'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)

    def test_budget_counts_across_files_and_blocks_before_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);budget=WriteBudget(10)
            with budget.open(tmp/'a') as f:f.write(b'123456')
            with budget.open(tmp/'b') as f:
                with self.assertRaises(StorageLimit):f.write(b'12345')
            self.assertEqual(directory_bytes(tmp),6)

    def test_scheduler_reserves_export_and_preserves_empty_study(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);config=app.read(app.DEFAULT)
            config['sizes']=[8];config['target_phis']=['0.44']
            config['benchmark'].update(max_rss_mib=128,min_free_disk_mib=16,max_case_output_mib=16)
            path=tmp/'config.json';path.write_text(json.dumps(config))
            command=[sys.executable,str(ROOT/'scripts/run_h2_validation.py'),'run','--config',str(path),
                     '--output',str(tmp/'study'),'--samples','1','--workers','1','--no-plots']
            result=subprocess.run(command+['--max-storage-gib','0.02'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Storage budget reached',result.stderr)
            self.assertFalse(list((tmp/'study').glob('cases/*/result.json')))
            subprocess.run(command+['--max-storage-gib','1'],check=True,stdout=subprocess.DEVNULL)
            self.assertEqual(len(list((tmp/'study').glob('cases/*/result.json'))),1)

    def test_zip_writes_cannot_bypass_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.zip';budget=WriteBudget(512)
            with self.assertRaises(StorageLimit):
                with budget.open(path) as target:
                    with ZipFile(target,'w') as z:
                        with z.open('payload','w') as member:
                            member.write(random.Random(4).randbytes(4096))
            self.assertLessEqual(path.stat().st_size,512)

    def test_gzip_csv_valid_and_bounded(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'clusters.csv.gz'
            rows=[dict(size=2,h2=1),dict(size=1,h2=0)]
            WriteBudget(1024).csv(p,rows,compressed=True)
            with gzip.open(p,'rt',encoding='utf-8-sig') as f:
                self.assertEqual(list(csv.DictReader(f)),[{'size':'2','h2':'1'},{'size':'1','h2':'0'}])

    def test_full_compact_same_bonds_and_measurements(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cfg=app.read(app.DEFAULT)
            cfg['sizes']=[8];cfg['target_phis']=['0.44']
            cfg['benchmark'].update(max_rss_mib=128,min_free_disk_mib=16,max_case_output_mib=16)
            case=app.jobs(cfg,1)[0]
            bcfg=dict(cfg['benchmark'],sizes=[8],target_phis=['0.44'],seeds=[1])
            results=[]
            for policy in ('compact','full'):
                folder=root/policy;folder.mkdir();(folder/'study.json').write_text(json.dumps({'storage_policy':policy}))
                path=Path(app.worker(str(folder),bcfg,case,'test',policy));results.append(app.read(path))
            for field in results[0]['sample']:
                if field not in ('wall_seconds','retained_bytes_before_manifest'):
                    self.assertEqual(results[0]['sample'][field],results[1]['sample'][field],field)
            self.assertEqual(results[0]['bins'],results[1]['bins'])
            for name in ('bonds.csv.gz','cluster-members.csv.gz','cluster-summary.csv.gz'):
                a=root/'compact/cases'/case['case_id']/'runs'/case['case_id']/name
                b=root/'full/cases'/case['case_id']/'runs'/case['case_id']/name
                self.assertEqual(gzip.decompress(a.read_bytes()),gzip.decompress(b.read_bytes()))
            compact=root/'compact/cases'/case['case_id']
            self.assertFalse(list(compact.rglob('trajectories.csv*')))
            self.assertFalse(list(compact.rglob('candidates.csv*')))
            self.assertTrue(list(compact.rglob('bonds.csv.gz')))
            self.assertFalse((compact/'clusters.csv').exists())
            out=app.export(root/'compact',[case],1,'test',True,maximum=8*1024**2,min_free=0)
            self.assertFalse((out/'cluster_samples.csv').exists())
            self.assertFalse((out/'EXPORT_INCOMPLETE.json').exists())
            with ZipFile(out/'upload_h2_validation.zip') as z:
                self.assertIsNone(z.testzip());self.assertIn('storage_report.json',z.namelist())
            self.assertLess(directory_bytes(root/'compact'),8*1024**2)
            # No new samples: force export itself to exhaust the available allocation.
            # An incompressible existing result payload defeats ZIP compression.
            result_path=compact/'result.json';value=app.read(result_path)
            value['sample']['padding']=random.Random(0).randbytes(2*1024**2).hex()
            result_path.write_text(json.dumps(value))
            maximum=directory_bytes(root/'compact')+1024*1024
            with self.assertRaises(StorageLimit):app.export(root/'compact',[case],1,'test',False,maximum=maximum,min_free=0)
            self.assertLessEqual(directory_bytes(root/'compact'),maximum)
            self.assertTrue(result_path.exists())

if __name__=='__main__':unittest.main()
