#!/usr/bin/env python3
"""Run the phi=0.41..0.49 study in resumable completed groups (Linux/WSL)."""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import benchmark
from benchmark_resources import atomic_json, digest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts/analysis'))
from h2_cluster.stats import summarize, write_csv

DEFAULT_CONFIG = ROOT/'configs/experiments/phi040_050.json'
DEFAULT_OUTPUT = ROOT/'data/refinement/phi040_050'


def load_plan(path):
    plan = json.loads(path.read_text(encoding='utf-8-sig'), object_pairs_hook=benchmark.no_duplicate_keys)
    if set(plan) != {'schema_version','benchmark','groups'} or plan['schema_version'] != 1:
        raise ValueError('Invalid study config keys or version')
    if not plan['groups']:
        raise ValueError('Study groups must not be empty')
    names, conditions = set(), set()
    groups = []
    for group in plan['groups']:
        if set(group) != {'name','sizes','target_phis'} or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',group['name']):
            raise ValueError('Invalid group name or settings')
        name = group['name']
        if name in names:
            raise ValueError('Duplicate group name')
        names.add(name)
        cfg = dict(plan['benchmark'], sizes=group['sizes'], target_phis=group['target_phis'])
        with tempfile.TemporaryDirectory() as tmp:
            file = Path(tmp)/'config.json'
            file.write_text(json.dumps(cfg))
            benchmark.load_config(file)
        if cfg['repeats'] != 1 or cfg['postprocess'] or cfg['plot']:
            raise ValueError('Study requires repeats=1, postprocess=false, plot=false; analysis follows separately')
        for size in cfg['sizes']:
            for phi in cfg['target_phis']:
                for seed in cfg['seeds']:
                    key = size,Decimal(str(phi)),seed
                    if key in conditions:
                        raise ValueError(f'Duplicate simulation condition: {key}')
                    conditions.add(key)
        groups.append((name,cfg))
    return plan,groups


def sample_count(cfg):
    return len(cfg['sizes'])*len(cfg['target_phis'])*len(cfg['seeds'])


def analysis_config(cfg):
    return {key:cfg[key] for key in ['schema_version','sizes','target_phis','seeds','arm_length','a_probability']}


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def run_group(root, name, cfg):
    group = root/name
    group.mkdir(exist_ok=True)
    bench_cfg = group/'benchmark-config.json'
    analysis_cfg = group/'analysis-config.json'
    for path,value in [(bench_cfg,cfg),(analysis_cfg,analysis_config(cfg))]:
        if path.exists() and read_json(path) != value:
            raise ValueError(f'{path}: settings changed; choose a new --output study directory')
        if not path.exists():atomic_json(path,value)
    simulation = group/'simulation'
    if simulation.exists():
        marker = simulation/'benchmark.json'
        if not marker.is_file() or read_json(marker).get('status')!='completed':
            raise ValueError(f'{simulation}: incomplete previous batch. It is preserved and will not be silently rerun; send the error and benchmark.json.')
        if read_json(simulation/'config.json') != cfg:
            raise ValueError(f'{simulation}: completed batch has different settings')
        print(f'Reuse completed simulation group: {name}',flush=True)
    else:
        subprocess.run([sys.executable,str(ROOT/'scripts/benchmark.py'),'--config',str(bench_cfg),
                        '--output',str(simulation)],cwd=ROOT,check=True)
    analysis = group/'analysis'
    if analysis.exists():
        marker = analysis/'run_manifest.json'
        if not marker.is_file() or read_json(marker).get('status')!='completed':
            raise ValueError(f'{analysis}: incomplete analysis; existing files preserved. Send the error message.')
        if read_json(analysis/'analysis_config.json') != analysis_config(cfg):
            raise ValueError(f'{analysis}: analysis config mismatch')
        if read_json(marker)['samples'] != sample_count(cfg):
            raise ValueError(f'{analysis}: incomplete sample count')
    else:
        subprocess.run([sys.executable,str(ROOT/'scripts/analysis/analyze_h2_clusters.py'),
                        '--config',str(analysis_cfg),'--suite',str(simulation),
                        '--output',str(analysis),'--no-plots'],cwd=ROOT,check=True)
    return analysis


def export_results(root, selected, mode):
    samples,histograms,sources = [],[],[]
    binary_hashes = set()
    conditions = set()
    for name,cfg in selected:
        folder = root/name/'analysis'
        with (folder/'h2_samples.csv').open(encoding='utf-8-sig',newline='') as stream:
            raw_samples = list(csv.DictReader(stream))
        if len(raw_samples) != sample_count(cfg):
            raise ValueError(f'{name}: sample count mismatch')
        expected = {(L,float(phi),seed) for L in cfg['sizes'] for phi in cfg['target_phis'] for seed in cfg['seeds']}
        observed = set()
        for raw in raw_samples:
            row = {k:(v if k=='case_id' else None if not v else float(v) if any(c in v.lower() for c in '.e') else int(v)) for k,v in raw.items()}
            key = row['L'],row['target_phi'],row['seed']
            if key not in expected or key in conditions:
                raise ValueError(f'Duplicate/unexpected result condition: {key}')
            observed.add(key);conditions.add(key)
            row['case_id'] = name+'/'+row['case_id']
            samples.append(row)
        if observed != expected:
            raise ValueError(f'{name}: missing condition')
        with (folder/'cluster_histograms.csv').open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                row['case_id'] = name+'/'+row['case_id']
                histograms.append(row)
        simulation_manifest=read_json(root/name/'simulation/benchmark.json')
        binary_hashes.add(simulation_manifest['binary_sha256'])
        sources.append(dict(group=name,path=str(folder),simulation_binary_sha256=simulation_manifest['binary_sha256'],
                            manifest=read_json(folder/'run_manifest.json'),
                            sample_sha256=digest(folder/'h2_samples.csv'),histogram_sha256=digest(folder/'cluster_histograms.csv')))
    if len(binary_hashes) != 1:
        raise ValueError('Simulation binaries differ across groups; inspect source/compiler changes before combining results')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output = root/('results_'+mode+'_'+stamp)
    output.mkdir()
    write_csv(output/'h2_samples.csv',samples)
    summary = summarize(samples)
    write_csv(output/'h2_summary.csv',summary)
    write_csv(output/'cluster_histograms.csv',histograms)
    manifest = dict(schema_version=1,status='completed',mode=mode,samples=len(samples),conditions=len(summary),
                    source='new detailed simulations; original phi=0.40/0.50 results are reused separately',
                    sources=sources,study_plan=read_json(root/'study.json')['plan'],
                    runner_sha256=digest(Path(__file__)),case_id='group/original case_id; model condition is (L,target_phi,seed)',
                    definitions=read_json(root/selected[0][0]/'analysis/run_manifest.json')['definitions'])
    atomic_json(output/'run_manifest.json',manifest)
    archive = output/('upload_phi_detail_'+mode+'.zip')
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for path in sorted(output.iterdir()):
            if path != archive:z.write(path,path.name)
    print(f'\nSUCCESS: {len(samples)} samples, {len(summary)} conditions\nUPLOAD THIS FILE:\n{archive}',flush=True)
    return archive


def run(mode, config, output):
    plan,groups = load_plan(config)
    selected = groups[:1] if mode=='pilot' else groups
    print('Detailed concentration study (old phi=0.40/0.50 samples are not regenerated)',flush=True)
    for name,cfg in groups:
        print(f"  {name}: L={cfg['sizes']}, phi={cfg['target_phis']}, {sample_count(cfg)} samples",flush=True)
    if mode=='plan':return None
    output = output.resolve()
    state_path = output/'study.json'
    if output.exists():
        if not state_path.is_file() or read_json(state_path).get('plan')!=plan:
            raise ValueError(f'{output}: unrecognized study directory or changed plan; choose a new --output')
    else:
        output.mkdir(parents=True)
        atomic_json(state_path,dict(schema_version=1,plan=plan,status='created'))
    # A crashed process leaves the lock for explicit inspection; never run concurrently in one study.
    lock = output/'running.lock'
    fd = os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        for name,cfg in selected:run_group(output,name,cfg)
        archive = export_results(output,selected,mode)
        atomic_json(state_path,dict(schema_version=1,plan=plan,status=mode+'_completed',latest_archive=str(archive)))
        return archive
    finally:
        lock.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['plan','pilot','full'],nargs='?',default='pilot')
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT,help='Study directory; completed groups are reused')
    args = parser.parse_args()
    try:run(args.mode,args.config,args.output)
    except (OSError,ValueError,KeyError,TypeError,InvalidOperation,subprocess.SubprocessError) as exc:
        parser.exit(1,f'ERROR: {exc}\nExisting experiment files were preserved. Send the complete error message.\n')


if __name__=='__main__':main()
