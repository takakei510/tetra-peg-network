#!/usr/bin/env python3
"""L128 all-phi molecular C4/C6 counts. Reuse validated DFS and hash-bound caches.
One sample is processed at a time. No new network generation; no partial plots.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import time
from zipfile import ZipFile, ZIP_DEFLATED

REPO=Path(__file__).resolve().parents[2]
VERSION='cycle46-all-phi-v1'
PINNED_ENGINE='d16eede57d87d37ae3578454847ea6655b4124bb6fdd19be35804f7dbbdd9faa'
COUNT_KEYS=('C4_mol','C6_mol','C4_edge','C6_edge','B4','B6')
METRICS=('N','H2','C4_mol','C6_mol','C4_per_100_molecules','C6_per_100_molecules')

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def atomic_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')
    temp.replace(path)

def load_engine(path,expected):
    if expected!=PINNED_ENGINE or digest(path)!=expected:
        raise ValueError('DFS engine hash changed; reviewed cache compatibility is required before running.')
    spec=importlib.util.spec_from_file_location('cycle46_engine',path)
    engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    return engine

def check_config(c):
    if c['schema_version']!=1 or c['scope']!='all' or c['samples_per_phi']<2 or c['L']<1:
        raise ValueError('Invalid configuration')
    phis=c['phi_values']
    if not phis or phis!=sorted(set(phis)) or any(not 0<float(p)<1 for p in phis):
        raise ValueError('phi_values must be increasing and unique, between 0 and 1')
    if not set(c['pilot_phi_values']).issubset(phis):raise ValueError('Pilot phi must belong to phi_values')

def select_samples(study,c,explicit=None):
    if explicit:source=Path(explicit)
    else:
        options=[]
        for f in study.glob('results_*'):
            if (f/'manifest.json').is_file() and not (f/'EXPORT_INCOMPLETE.json').exists():
                m=json.loads((f/'manifest.json').read_text())
                if m.get('status')=='completed' and m.get('completed_samples')==m.get('requested_samples') and (f/'samples.csv').is_file():
                    options.append((int(m['completed_samples']),f.stat().st_mtime,f/'samples.csv'))
        if not options:raise ValueError('No completed export; specify --samples-csv')
        source=max(options)[2]
    phis={round(float(p),8) for p in c['phi_values']};groups=defaultdict(list)
    with source.open(encoding='utf-8-sig',newline='') as f:
        for r in csv.DictReader(f):
            phi=round(float(r['target_phi']),8)
            if int(r['L'])!=c['L'] or phi not in phis:continue
            if abs(float(r['target_phi'])-phi)>1e-10:raise ValueError('Ambiguous phi precision')
            case=r['case_id']
            if Path(case).name!=case or case in ('','.','..'):raise ValueError('Unsafe case ID')
            groups[phi].append(r)
    if set(groups)!=phis:raise ValueError('Missing target phi conditions')
    ids=[]
    for phi,rr in groups.items():
        if len(rr)!=c['samples_per_phi'] or len({r['seed'] for r in rr})!=len(rr):
            raise ValueError(f'phi={phi}: expected {c["samples_per_phi"]} distinct trials/seeds')
        rr.sort(key=lambda r:r['case_id']);ids.extend(r['case_id'] for r in rr)
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate case IDs')
    return source,groups

def raw_hashes(engine,study,sample):
    folder=study/'cases'/sample['case_id']/'runs'/sample['case_id']
    return dict(bonds=digest(engine.raw_path(folder,'bonds.csv')),members=digest(engine.raw_path(folder,'cluster-members.csv')))

def compact_totals(rec,sample):
    rows=rec['clusters']
    if not rows:raise ValueError('Empty cluster cache')
    measured=(len(rows),sum(r['size'] for r in rows),max(r['size'] for r in rows),sum(r['H2'] for r in rows),sum(r['bonds'] for r in rows))
    if measured!=tuple(int(sample[k]) for k in ('clusters','N','N1','H2','bonds')):
        raise ValueError('Cluster cache totals disagree with sample')
    out={k:sum(r[k] for r in rows) for k in COUNT_KEYS}
    validate_totals(out)
    return out

def validate_totals(out):
    if set(out)!=set(COUNT_KEYS) or any(not isinstance(out[k],int) or isinstance(out[k],bool) or out[k]<0 for k in COUNT_KEYS):
        raise ValueError('Invalid compact cycle counts')
    if out['C4_edge']<out['C4_mol'] or out['C6_edge']<out['C6_mol']:raise ValueError('Multiplicity counts inconsistent')

def process_case(engine,study,sample,engine_hash):
    # Rehash raw files before accepting either cache. Existing 0.43 full caches stay untouched.
    hashes=raw_hashes(engine,study,sample)
    cp=study/'cycle46_all_phi_cache'/engine_hash[:16]/(sample['case_id']+'.json')
    binding=dict(version=VERSION,engine_sha256=engine_hash,sample=sample,input_sha256=hashes)
    if cp.is_file():
        old=json.loads(cp.read_text())
        if all(old.get(k)==v for k,v in binding.items()):
            validate_totals(old['totals']);return old['totals'],hashes,'compact-cache'
    legacy=study/'cycle46_cache'/engine_hash[:16]/(sample['case_id']+'.json')
    if legacy.is_file():
        old=json.loads(legacy.read_text())
        if old.get('version')==engine.VERSION and old.get('sample')==sample and old.get('input_sha256')==hashes:
            totals=compact_totals(old,sample);origin='existing-cycle46-cache'
        else:totals=None
    else:totals=None
    if totals is None:
        rec=engine.inspect_case(study,sample)
        if rec['input_sha256']!=hashes:raise ValueError('Raw files changed while processing')
        totals=compact_totals(rec,sample);origin='new-dfs'
    atomic_json(cp,{**binding,'totals':totals})
    return totals,hashes,origin

def total_row(sample,totals):
    row={k:sample[k] for k in ('case_id','L','target_phi','seed','N','N1','bonds','clusters','H2')}
    row['actual_phi']=float(sample.get('actual_phi',21*int(sample['N'])/int(sample['L'])**3))
    row.update(totals)
    for k in (4,6):row[f'C{k}_per_100_molecules']=100*totals[f'C{k}_mol']/int(sample['N'])
    return row

def summarize(rows,config):
    groups=defaultdict(list)
    for r in rows:groups[round(float(r['target_phi']),8)].append(r)
    if set(groups)!={round(p,8) for p in config['phi_values']}:raise ValueError('Incomplete phi coverage')
    out=[]
    for phi,rr in sorted(groups.items()):
        if len(rr)!=config['samples_per_phi']:raise ValueError('Incomplete trial count')
        item=dict(L=config['L'],target_phi=phi,n_samples=len(rr),actual_phi_mean=statistics.mean(float(r['actual_phi']) for r in rr))
        for metric in METRICS:
            v=[float(r[metric]) for r in rr];sd=statistics.stdev(v)
            item.update({metric+'_mean':statistics.mean(v),metric+'_sd':sd,metric+'_se':sd/math.sqrt(len(v)),metric+'_min':min(v),metric+'_max':max(v)})
        out.append(item)
    return out

def write_csv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def draw(summary,output,config):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    specs=[('01_cycle_counts_vs_phi',('C4_mol','C6_mol'),'Cycles per trial [count]'),('02_cycle_density_vs_phi',('C4_per_100_molecules','C6_per_100_molecules'),'Cycles per 100 molecules [count / 100 molecules]')]
    for name,metrics,label in specs:
        fig,ax=plt.subplots(figsize=(9.5,6))
        for metric,color,legend in zip(metrics,('#BD681B','#218477'),('4-molecule cycles','6-molecule cycles')):
            ax.errorbar([r['actual_phi_mean'] for r in summary],[r[metric+'_mean'] for r in summary],yerr=[r[metric+'_se'] for r in summary],fmt='o',markersize=3.5,capsize=2,color=color,label=legend)
        ax.set_xlabel('Actual occupied fraction phi');ax.set_ylabel(label);ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2);ax.legend()
        ax.set_title(f'L={config["L"]}; all clusters including largest; {config["samples_per_phi"]} trials / condition')
        fig.text(.1,.015,'Mean +/- SE across independent trials. Molecular routes counted once under rotation/reversal.\nChords allowed; parallel edge choices collapsed. Not independent cycle rank. No causal inference.',fontsize=9)
        fig.tight_layout(rect=(0,.065,1,1))
        for ext in ('png','svg'):fig.savefig(output/(name+'.'+ext),dpi=180)
        plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=REPO/'configs/analysis/cycle-46-all-phi.json')
    ap.add_argument('--study',type=Path);ap.add_argument('--samples-csv',type=Path)
    ap.add_argument('--engine',type=Path,default=Path(__file__).with_name('analyze_cycles_46.py'))
    ap.add_argument('--mode',choices=('pilot','run'),default='pilot');ap.add_argument('--output',type=Path)
    args=ap.parse_args();config=json.loads(args.config.read_text());check_config(config)
    engine=load_engine(args.engine,config['engine_sha256'])
    study=(args.study or REPO/config['study']).resolve()
    source,groups=select_samples(study,config,args.samples_csv)
    source_hash=digest(source)
    if args.mode=='run':
        import matplotlib  # Fail before DFS if plotting dependency is missing.
    chosen=[groups[round(p,8)][0] for p in config['pilot_phi_values']] if args.mode=='pilot' else [r for p in sorted(groups) for r in groups[p]]
    # Check required raw paths before starting DFS. Does not load graph data.
    for sample in chosen:
        folder=study/'cases'/sample['case_id']/'runs'/sample['case_id']
        engine.raw_path(folder,'bonds.csv');engine.raw_path(folder,'cluster-members.csv')
    started=time.perf_counter();rows=[];inputs=[];origins=defaultdict(int)
    for i,sample in enumerate(chosen,1):
        begin=time.perf_counter();totals,hashes,origin=process_case(engine,study,sample,config['engine_sha256']);origins[origin]+=1
        rows.append(total_row(sample,totals));inputs.append(dict(case_id=sample['case_id'],sha256=hashes))
        if args.mode=='pilot' or i%config['samples_per_phi']==0:
            print(f'{i}/{len(chosen)} {sample["case_id"]}: C4={totals["C4_mol"]}, C6={totals["C6_mol"]}, last_case={time.perf_counter()-begin:.2f}s, elapsed={time.perf_counter()-started:.1f}s, sources={dict(origins)}',flush=True)
    if args.mode=='pilot':
        print('PILOT ONLY: one trial per selected phi; no 100-trial figures. Runtime over the full sweep is not guaranteed.\nNext: rerun with --mode run. Checkpoints are reused after raw-file hashing.')
        return
    if digest(source)!=source_hash:raise ValueError('Input samples.csv changed during analysis')
    summary=summarize(rows,config)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=args.output or study/'cycle46_all_phi_outputs'/stamp;output.mkdir(parents=True,exist_ok=False)
    write_csv(output/'sample_totals.csv',rows);write_csv(output/'phi_summary.csv',summary)
    meta=dict(status='completed',created_at_utc=datetime.now(timezone.utc).isoformat(),config=config,mode=args.mode,samples=len(rows),conditions=len(summary),wall_seconds=time.perf_counter()-started,cache_sources=dict(origins),script_sha256=digest(__file__),engine_sha256=config['engine_sha256'],samples_sha256=source_hash,raw_inputs=inputs,
              definitions=dict(scope='all clusters including one largest and isolates',counts='distinct-molecule simple C4/C6 routes, modulo rotation/reversal; chords allowed; parallel pair choices collapsed',normalization='100 * trial cycle count / trial N; then equal trial mean',uncertainty='sample SD and SE=SD/sqrt(n), zeros included',warning='not endpoint percentage, not independent cycle rank, no causal growth conclusion'))
    atomic_json(output/'run_metadata.json',meta)
    (output/'figure_guide.md').write_text('''# L=128 全占有率の4・6分子閉路
目的：占有率によって短い閉路の存在量がどう変わるかを確認する。
01：全クラスタ（最大クラスタを含む）の1試料あたり閉路個数。
02：試料ごとに100*Ck/Nを計算後に平均。単位は個/100分子で、百分率・末端割合ではない。
各点は独立試料平均±SE。0個の試料も含む。横軸は実占有率、CSVには目標占有率も保存。
閉路の定義とDFSは既存phi=0.43解析と同じ。弦ありも含み、開始点・向きの重複を除き、並列辺の選び分けは主個数では区別しない。
異なるphiは別の試料であり、同じクラスタの時間発展ではない。閉路増加から成長抑制の因果は言えない。全クラスタの個数と、以前の最大クラスタ除外・サイズ別平均を混同しない。
処理は各試料ごとに再開可能。既存の0.43キャッシュと、新しい小さい総数キャッシュを、サンプル定義・元ファイルSHA256・DFSコードSHA256が合う場合に再利用する。今回の入力CSVとキャッシュ由来はmetadataに保存。
''',encoding='utf-8')
    draw(summary,output,config)
    archive=output/'cycle46_all_phi_share.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for f in sorted(output.iterdir()):
            if f!=archive:z.write(f,f.name)
    print('OUTPUT:',output,flush=True);print('SHARE:',archive,flush=True)

if __name__=='__main__':main()
