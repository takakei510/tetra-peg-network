#!/usr/bin/env python3
"""Analyze existing benchmark-style runs. Never starts a simulation."""
import argparse
import json
import math
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

from h2_cluster.metrics import h2_metrics
from h2_cluster.raw import locate, iter_csv, cluster_metrics, count_bonds, filehash
from h2_cluster.stats import summarize, write_csv

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / 'configs/analysis/h2-cluster-report.json'


def load_config(path):
    config = json.loads(path.read_text(encoding='utf-8-sig'))
    expected = {'schema_version','sizes','target_phis','seeds','arm_length','a_probability'}
    if set(config) != expected or config['schema_version'] != 1:
        raise ValueError('Unsupported analysis config keys or version')
    for key in ['sizes','seeds']:
        values = config[key]
        if not isinstance(values,list) or not values or any(type(v) is not int or v < (1 if key=='sizes' else 0) for v in values) or len(set(values))!=len(values):
            raise ValueError(f'Invalid or duplicate {key}')
    phis = [Decimal(str(v)) for v in config['target_phis']]
    if not phis or any(not v.is_finite() or not 0 < v <= 1 for v in phis) or len(set(phis))!=len(phis):
        raise ValueError('Invalid or duplicate target_phis')
    if type(config['arm_length']) is not int or config['arm_length']<1:
        raise ValueError('arm_length must be a positive integer')
    prob = Decimal(str(config['a_probability']))
    if not prob.is_finite() or not 0 <= prob <= 1:
        raise ValueError('Invalid a_probability')
    cases = []
    for size in config['sizes']:
        for phi in phis:
            for seed in config['seeds']:
                cases.append(dict(case_id=f'case{len(cases)+1:05d}_L{size}_seed{seed}_r1',L=size,target_phi=float(phi),seed=seed))
    return config, cases


def resolve_runs(root, cases, explicit=None):
    if explicit is not None:
        candidates = [explicit/'runs', explicit]
    else:
        search = root/'data'
        candidates = [p.parent for p in search.rglob(cases[0]['case_id']) if p.is_dir()] if search.is_dir() else []
    valid = {p.resolve() for p in candidates if all((p/c['case_id']).is_dir() for c in cases)}
    if len(valid)!=1:
        raise ValueError(f'Expected one directory with all {len(cases)} runs; found {len(valid)}. Use --suite PATH to select the original simulation/runs folder.')
    return valid.pop()


def read_cfg(path):
    config = {}
    for raw in path.read_text(encoding='utf-8-sig').splitlines():
        line = raw.split('#',1)[0].strip()
        if not line:continue
        key,value = (v.strip() for v in line.split('=',1))
        if not key or key in config:raise ValueError(f'{path}: duplicate or empty setting')
        config[key] = value
    return config


def one_row(path):
    rows = list(iter_csv(path))
    if len(rows)!=1:raise ValueError(f'{path}: expected exactly one summary row')
    return rows[0]


def analyze_case(folder, case, config):
    used = [folder/'config.cfg',folder/'run.json']
    cfg = read_cfg(used[0])
    required = {'dim':3,'L':case['L'],'seed':case['seed'],'arm_length':config['arm_length']}
    if any(int(cfg[k])!=v for k,v in required.items()):
        raise ValueError(f'{folder}: model settings differ from the analysis config')
    for key,value in [('target_phi',case['target_phi']),('a_probability',config['a_probability'])]:
        if Decimal(cfg[key]) != Decimal(str(value)):
            raise ValueError(f'{folder}: {key} differs from the analysis config')
    if json.loads(used[1].read_text())['status']!='completed':
        raise ValueError(f'{folder}: run is not completed')
    np = locate(folder,'network-summary.csv'); used.append(np)
    net = one_row(np)
    n,e,k,n1 = [int(net[key]) for key in ['molecules','bonds','clusters','largest']]
    if not (n>0 and 0<=e<=2*n and 1<=k<=n and 1<=n1<=n):
        raise ValueError(f'{folder}: invalid network counts')
    if not math.isclose(float(net['largest_fraction']),n1/n,rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError(f'{folder}: largest_fraction mismatch')
    pp = locate(folder,'placement-summary.csv'); used.append(pp)
    placement = one_row(pp)
    if placement['completed']!='1':raise ValueError(f'{folder}: placement incomplete')
    if int(placement['placed_molecules'])!=n or int(placement['target_molecules'])!=n:
        raise ValueError(f'{folder}: placement molecule count mismatch')
    phi = (1+4*config['arm_length'])*n/case['L']**3
    if not math.isclose(float(placement['actual_phi']),phi,rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError(f'{folder}: actual_phi mismatch')
    cp=locate(folder,'cluster-summary.csv');bp=locate(folder,'bonds.csv');used.extend([cp,bp])
    cm,hist = cluster_metrics(iter_csv(cp),n,e,k,n1)
    h = count_bonds(iter_csv(bp),n,e)
    if sum((i+1)*v for i,v in enumerate(h))!=e:raise ValueError('Sum(k Hk) != E')
    row = dict(**case,actual_phi=phi,N=n,bonds=e,clusters=k,N1=n1,Pmax=n1/n,
               S1=(1+4*config['arm_length'])*n1/case['L']**3,p=e/(2*n),
               H2=h[1],**h2_metrics(h[1],n,e,case['L']**3),**cm)
    if not math.isclose(row['H2_endpoint_fraction_all'],row['p']*(row['H2_endpoint_fraction_bonded'] or 0),abs_tol=1e-12):
        raise ValueError('Endpoint fraction identity failed')
    source = dict(case_id=case['case_id'],files=[dict(path=str(p),sha256=filehash(p)) for p in used],
                  validation_Hk={f'H{i+1}':v for i,v in enumerate(h)})
    histogram=[dict(case_id=case['case_id'],L=case['L'],target_phi=case['target_phi'],seed=case['seed'],size=s,count=count) for s,count in sorted(hist.items())]
    return row,histogram,source


def git_revision():
    try:
        return subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError):return None


def run(config_path, suite=None, output=None, no_plots=False):
    config,cases=load_config(config_path)
    runs=resolve_runs(ROOT,cases,suite)
    out=output or ROOT/'data/analysis'/('h2_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8])
    if out.exists():raise ValueError('Output exists; select a new --output. Existing results are never overwritten.')
    rows=[];hist=[];sources=[]
    print(f'Found {len(cases)} runs: {runs}',flush=True)
    for i,case in enumerate(cases,1):
        print(f"[{i}/{len(cases)}] {case['case_id']} phi={case['target_phi']}",flush=True)
        try:row,histogram,source=analyze_case(runs/case['case_id'],case,config)
        except (OSError,ValueError,KeyError) as exc:
            raise ValueError(f"{case['case_id']}: {exc}") from exc
        rows.append(row);hist.extend(histogram);sources.append(source)
    summary=summarize(rows)
    out.mkdir(parents=True,exist_ok=False)
    write_csv(out/'h2_samples.csv',rows);write_csv(out/'h2_summary.csv',summary)
    write_csv(out/'cluster_histograms.csv',hist)
    (out/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    plot_status='disabled'
    if not no_plots:
        try:
            from h2_cluster.plots import make_plots
            make_plots(summary,out)
            plot_status='completed'
        except ImportError:
            plot_status='matplotlib_unavailable'
            print('Matplotlib unavailable: CSV results and ZIP are still produced.',flush=True)
    manifest = dict(schema_version=1,status='completed',samples=len(rows),conditions=len(summary),
        source='original run CSVs; no aggregate input tables required',sources=sources,
        analysis_git_commit=git_revision(),
        analysis_source_sha256={str(p.relative_to(ROOT)):filehash(p) for p in [Path(__file__).resolve(), *sorted((Path(__file__).parent/'h2_cluster').glob('*.py'))]},
        python=sys.version,plot_status=plot_status,
        config_sha256=filehash(config_path),
        definitions={'H2':'unordered molecular pairs with exactly two bonds; H3/H4 pairs excluded',
        'H2_per_volume':'H2/L^3','H2_endpoint_fraction_bonded':'2*H2/E',
        'H2_endpoint_fraction_all':'H2/N','actual_phi':'(1+4*l)*N/L^3',
        'cluster_weighted_excl_one':'sum(s^2)/sum(s) excluding one largest cluster per sample',
        'ties':'exclude smallest cluster_id; also export all-tied exclusion sensitivity',
        'undefined':'zero denominator is null/blank, never zero','statistics':'mean of per-sample statistics; sample SD ddof=1'})
    (out/'run_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    zipped=out/'upload_h2_cluster_results.zip'
    with ZipFile(zipped,'w',ZIP_DEFLATED) as z:
        for p in sorted(out.iterdir()):
            if p.is_file() and p!=zipped:z.write(p,p.name)
    print(f'\nSUCCESS: {len(rows)} samples, {len(summary)} conditions\nUPLOAD THIS FILE:\n{zipped.resolve()}',flush=True)
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--suite',type=Path,help='Existing simulation or runs directory; default: discover under repository data/')
    parser.add_argument('--output',type=Path,help='New output directory (default: data/analysis/h2_TIMESTAMP_UUID)')
    parser.add_argument('--no-plots',action='store_true',help='Only CSV/JSON/ZIP; no Matplotlib needed')
    args=parser.parse_args()
    try:run(args.config,args.suite,args.output,args.no_plots)
    except (OSError,ValueError,KeyError,TypeError,InvalidOperation) as exc:
        parser.exit(1,f'ERROR: {exc}\nNo new simulations were run. Send the complete error message.\n')


if __name__=='__main__':main()
