#!/usr/bin/env python3
"""Resumable, process-parallel independent-seed H2/cluster validation (Linux/WSL)."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
import uuid
from zipfile import ZipFile, ZIP_DEFLATED
import benchmark
from benchmark_resources import atomic_json, digest, available_memory, directory_bytes, MIB, environment
from run_demo import snapshot_sources
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts/analysis'))
from analyze_h2_clusters import analyze_case
from cluster_h2_relation import inspect, bin_rows
from h2_cluster.stats import write_csv
DEFAULT = ROOT/'configs/experiments/h2-validation.json'


def read(path):
    return json.loads(path.read_text())


def jobs(config, count):
    result = []
    for L in config['sizes']:
        for phi in config['target_phis']:
            for index in range(count):
                # Stable across extensions and scheduling; different stream per condition.
                seed = int.from_bytes(hashlib.sha256(f"{config['seed_namespace']}:{L}:{phi}:{index}".encode()).digest()[:8], 'big')
                result.append(dict(case_id=f'L{L}_phi{phi}_i{index:06d}', L=L,
                                   target_phi=float(phi), seed=seed, index=index))
    if len({c['seed'] for c in result}) != len(result):
        raise ValueError('Seed collision')
    return result


def worker(root_text, cfg, case, identity):
    root = Path(root_text); suite = root/'cases'/case['case_id']
    suite.mkdir(parents=True, exist_ok=True)
    destination = suite/'result.json'
    if destination.exists():
        value = read(destination)
        if value['identity'] != identity or value['case'] != case:
            raise ValueError('Cached sample identity mismatch')
        for path, expected in value['hashes'].items():
            if digest(suite/path) != expected:
                raise ValueError('Cached input/output hash mismatch: '+path)
        return str(destination)
    if (suite/'runs'/case['case_id']).exists():
        raise ValueError(f'Incomplete case preserved at {suite}; use --retry-failed to move aside and retry same seed')
    start = time.perf_counter()
    result = benchmark.run_case(cfg, suite, case['case_id'], case['L'], case['target_phi'], case['seed'], 1)
    if result['status'] != 'completed':
        raise ValueError(f"{case['case_id']}: {result['status']} {result.get('stop_reason')}")
    folder = suite/'runs'/case['case_id']
    acfg = dict(arm_length=cfg['arm_length'], a_probability=cfg['a_probability'])
    sample, histogram, source = analyze_case(folder, {k:v for k,v in case.items() if k!='index'}, acfg)
    clusters = inspect(folder, sample)
    write_csv(suite/'clusters.csv', clusters)
    bins = bin_rows(clusters)
    sample['K_per_volume'] = sample['clusters']/case['L']**3
    sample['wall_seconds'] = time.perf_counter()-start
    hashes = {str(Path(f['path']).relative_to(suite)):f['sha256'] for f in source['files']}
    # Membership was additionally checked by BFS.
    from h2_cluster.raw import locate
    member = locate(folder, 'cluster-members.csv')
    hashes[str(member.relative_to(suite))] = digest(member)
    hashes['clusters.csv'] = digest(suite/'clusters.csv')
    atomic_json(destination, dict(identity=identity, case=case, sample=sample, bins=bins, hashes=hashes))
    return str(destination)


def stats(values):
    n = len(values)
    sd = statistics.stdev(values) if n>1 else None
    return dict(n=n, mean=statistics.mean(values) if n else None, sd=sd,
                se=sd/math.sqrt(n) if sd is not None else None)


def export(root, cases, requested, identity, plots):
    out = root/('results_'+uuid.uuid4().hex[:10]); out.mkdir()
    records = [read(root/'cases'/c['case_id']/'result.json') for c in cases
               if (root/'cases'/c['case_id']/'result.json').exists()]
    if not records:
        raise ValueError('No completed samples to export')
    samples = [r['sample'] for r in records]; write_csv(out/'samples.csv', samples)
    grouped = {}; bins = {}; convergence = []
    for r in records:
        key = r['case']['L'], r['case']['target_phi']
        grouped.setdefault(key, []).append(r)
        for b in r['bins']:
            bk = (*key, b['scope'], b['size_low'], b['size_high'])
            bins.setdefault(bk, []).append(b)
    summary = []
    metrics = ['H2','H2_per_volume','H2_endpoint_fraction_all','H2_endpoint_fraction_bonded',
               'clusters','K_per_volume','Pmax','S1','p','cluster_weighted_excl_one','wall_seconds']
    for (L, phi), group in sorted(grouped.items()):
        row = dict(L=L, target_phi=phi, actual_phi=group[0]['sample']['actual_phi'], n_samples=len(group))
        for metric in metrics:
            values = [r['sample'][metric] for r in group if r['sample'][metric] is not None]
            row.update({metric+'_'+k:v for k,v in stats(values).items()})
        summary.append(row)
        ordered = sorted(group, key=lambda r:r['case']['index'])
        checkpoints = sorted(set([5,10,20,50,100,200,500,1000,len(group)]))
        for n in checkpoints:
            if n>len(group): continue
            for metric in metrics[:-1]:
                values = [r['sample'][metric] for r in ordered[:n] if r['sample'][metric] is not None]
                convergence.append(dict(L=L,target_phi=phi,first_n=n,metric=metric,**stats(values)))
    write_csv(out/'summary.csv', summary); write_csv(out/'convergence.csv', convergence)
    conditional = []
    for (L, phi, scope, low, high), bs in sorted(bins.items()):
        conditional.append(dict(L=L,target_phi=phi,scope=scope,size_low=low,size_high=high,
            total_samples=len(grouped[(L,phi)]), samples_with_bin=len(bs),
            clusters=sum(b['cluster_count'] for b in bs),
            **{'fraction_'+k:v for k,v in stats([b['fraction_mean'] for b in bs]).items()}))
    write_csv(out/'size_relation.csv', conditional)
    with (out/'cluster_samples.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer = None
        for r in records:
            with (root/'cases'/r['case']['case_id']/'clusters.csv').open(encoding='utf-8-sig',newline='') as source:
                for c in csv.DictReader(source):
                    row = {k:r['case'][k] for k in ('case_id','L','target_phi','seed')}; row.update(c)
                    if writer is None:
                        writer=csv.DictWriter(stream,fieldnames=list(row)); writer.writeheader()
                    writer.writerow(row)
    plot_status = 'disabled'
    if plots:
        try:
            make_plots(out, summary, conditional); plot_status='completed'
        except ImportError:
            plot_status='matplotlib_unavailable'
    atomic_json(out/'manifest.json', dict(identity=identity,requested_samples=requested,
        completed_samples=len(records),status='completed' if len(records)==requested else 'partial',
        plot_status=plot_status,study=read(root/'study.json'),
        definitions={'h2':'exactly two bonds per molecular pair; H3/H4 excluded',
        'fraction':'h2_i/s_i = 4h2_i/(4s_i)',
        'conditional':'mean of per-sample cluster-equal mean in fixed size bin; empty bins omitted, not zero',
        'se':'sample SD / sqrt(number of independent seeds contributing)',
        'limitations':'association only; denominator size can induce trends; no causality; no equilibrium SAW guarantee'}))
    archive=out/'upload_h2_validation.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for p in out.iterdir():
            if p != archive:z.write(p,p.name)
    print(f'Export: {len(records)}/{requested} samples\nUPLOAD THIS FILE:\n{archive}',flush=True)
    return out


def make_plots(out, summary, conditional):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Points + SE only. No interpolation, fit, or connecting lines.
    for metric in ['H2_per_volume','H2_endpoint_fraction_all','clusters','K_per_volume','cluster_weighted_excl_one','S1']:
        fig,ax=plt.subplots(layout='constrained')
        for L in sorted({r['L'] for r in summary}):
            rows=[r for r in summary if r['L']==L]
            ax.errorbar([r['actual_phi'] for r in rows],[r[metric+'_mean'] for r in rows],
                yerr=[r[metric+'_se'] or 0 for r in rows],fmt='o',linestyle='none',capsize=3,label=f'L={L}')
        ax.set(xlabel='Actual occupied fraction',ylabel=metric,title='Mean +/- SE (n=1: SE unavailable)')
        ax.legend();fig.savefig(out/(metric+'.png'),dpi=160);plt.close(fig)
    for L,phi in sorted({(r['L'],r['target_phi']) for r in conditional}):
        fig,ax=plt.subplots(layout='constrained')
        for scope in ['all','exclude_one_largest']:
            rows=[r for r in conditional if (r['L'],r['target_phi'],r['scope'])==(L,phi,scope)]
            ax.errorbar([math.sqrt(r['size_low']*r['size_high']) for r in rows],
                [r['fraction_mean'] for r in rows],yerr=[r['fraction_se'] or 0 for r in rows],
                fmt='o',linestyle='none',capsize=3,label=scope)
        ax.set(xscale='log',xlabel='Cluster size bin (geometric midpoint)',ylabel='Mean h2_i/s_i',
               title=f'L={L}, target phi={phi}: seed mean +/- SE')
        ax.legend();fig.savefig(out/f'relation_L{L}_phi{phi}.png',dpi=160);plt.close(fig)


def run(args):
    config=read(args.config)
    if set(config)!={'schema_version','seed_namespace','sizes','target_phis','benchmark'} or config['schema_version']!=1:
        raise ValueError('Invalid study config')
    if args.samples<1 or args.workers<1 or args.max_storage_gib<=0:
        raise ValueError('samples/workers/storage must be positive')
    # Reuse existing strict benchmark configuration validation.
    import tempfile
    cfg=dict(config['benchmark'],sizes=config['sizes'],target_phis=config['target_phis'],seeds=[1])
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'cfg.json';path.write_text(json.dumps(cfg)); benchmark.load_config(path)
    if cfg['repeats']!=1 or cfg['plot'] or cfg['postprocess'] or cfg['retention']!='gzip':
        raise ValueError('Require repeats=1, plot/postprocess=false, retention=gzip')
    cases=jobs(config,args.samples)
    if args.mode=='pilot':
        cases=[c for c in cases if c['L']==max(config['sizes']) and c['target_phi']==min((float(x) for x in config['target_phis']),key=lambda x:abs(x-0.44))][:args.samples]
    print(f'{len(cases)} planned samples, up to {args.workers} workers; raw CSVs retained as gzip',flush=True)
    if args.mode=='plan':return
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=True)
    import fcntl
    with (root/'running.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        subprocess.run(['make','-B','build/demo_network','build/check_network_config'],cwd=ROOT,check=True)
        source_paths=sorted((ROOT/'src').glob('*.c'))+sorted((ROOT/'include').glob('*.h'))+sorted((ROOT/'scripts').rglob('*.py'))+[ROOT/'Makefile']
        state=dict(config=config,binary=digest(ROOT/'build/demo_network'),
                   source_hashes={str(p.relative_to(ROOT)):digest(p) for p in source_paths})
        identity=hashlib.sha256(json.dumps(state,sort_keys=True).encode()).hexdigest()
        if (root/'study.json').exists() and read(root/'study.json')!=state:
            raise ValueError('Changed source/binary/config; choose a new --output')
        atomic_json(root/'study.json',state)
        atomic_json(root/'environment.json',environment())
        if not (root/'source.zip').exists():snapshot_sources(root/'source.zip')
        pending=[]
        for c in cases:
            suite=root/'cases'/c['case_id']
            if (suite/'result.json').exists():
                worker(str(root),cfg,c,identity) # validate hashes even on resume/export
            elif args.mode!='analyze':
                if (suite/'runs'/c['case_id']).exists() and args.retry_failed:
                    suite.rename(suite.with_name(suite.name+'_failed_'+uuid.uuid4().hex[:8]))
                pending.append(c)
        failure=None
        if pending:
            # Whole worker includes Python CSV analysis, so reserve memory per worker.
            free=available_memory()
            workers=min(args.workers,os.cpu_count() or 1)
            if free is not None:workers=min(workers,max(1,free//(int(cfg['max_rss_mib'])*MIB)))
            print(f'Using {workers} worker processes; Ctrl-C preserves completed samples.',flush=True)
            with ProcessPoolExecutor(max_workers=workers) as pool:
                iterator=iter(pending); active={}
                def submit():
                    nonlocal failure
                    # Reserve worst-case active raw output. Soft guard; not OS quota.
                    reserve=(len(active)+1)*cfg['max_case_output_mib']*MIB*2
                    if (directory_bytes(root)+reserve>args.max_storage_gib*1024**3 or
                        shutil.disk_usage(root).free<cfg['min_free_disk_mib']*MIB+reserve):
                        failure='Storage budget reached; increase --max-storage-gib or free space and rerun';return False
                    c=next(iterator,None)
                    if c is None:return False
                    active[pool.submit(worker,str(root),cfg,c,identity)]=c;return True
                for _ in range(workers):
                    if not submit():break
                while active:
                    done,_=wait(active,timeout=5,return_when=FIRST_COMPLETED)
                    for future in done:
                        c=active.pop(future)
                        try:
                            future.result(); print(f"Done {c['case_id']}",flush=True)
                        except Exception as exc:
                            failure=str(exc);print('ERROR: '+failure,flush=True)
                    if not failure:
                        while len(active)<workers and submit():pass
        if failure and not any((root/'cases'/c['case_id']/'result.json').exists() for c in cases):
            raise ValueError(failure)
        out=export(root,cases,len(cases),identity,not args.no_plots)
        if failure:raise ValueError(failure+'; partial export at '+str(out))
        if read(out/'manifest.json')['status']!='completed' and args.mode!='analyze':
            raise ValueError('Incomplete study; see manifest')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['plan','pilot','run','analyze'])
    p.add_argument('--config',type=Path,default=DEFAULT)
    p.add_argument('--output',type=Path,default=ROOT/'data/validation/h2_relation')
    p.add_argument('--samples',type=int,default=100,help='Independent seeds per condition; extend to 500/1000 with same output')
    p.add_argument('--workers',type=int,default=max(1,min(4,(os.cpu_count() or 2)//2)))
    p.add_argument('--max-storage-gib',type=float,default=20)
    p.add_argument('--retry-failed',action='store_true')
    p.add_argument('--no-plots',action='store_true')
    args=p.parse_args()
    try:run(args)
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as exc:p.exit(1,f'ERROR: {exc}\n')

if __name__=='__main__':main()
