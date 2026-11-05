#!/usr/bin/env python3
"""Count H2-positive clusters from validated cycle46-v1 cache; no new DFS.
H2 is a molecular pair with exactly two bonds, not all short cycles.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
from datetime import datetime, timezone
from zipfile import ZipFile, ZIP_DEFLATED

CACHE_HASH = 'd16eede57d87d37a'

def bounds(size):
    if size == 1:
        return 1, 1
    upper = 1 << (size - 1).bit_length()
    return (2 if upper == 2 else upper // 2 + 1), upper

def stats(values):
    values = list(values)
    n = len(values)
    return (statistics.mean(values) if n else None,
            statistics.stdev(values) / math.sqrt(n) if n > 1 else None, n)

def load_cache(cache, L=128, phi=.43, expected=100):
    records = []
    for path in sorted(Path(cache).glob('*.json')):
        rec = json.loads(path.read_text())
        sample = rec.get('sample', {})
        if int(sample.get('L', -1)) != L or not math.isclose(float(sample.get('target_phi', -1)), phi, abs_tol=1e-9):
            continue
        if rec.get('version') != 'cycle46-v1':
            raise ValueError(f'{path}: unsupported cache version')
        rows = rec['clusters']
        ids = [int(r['cluster_id']) for r in rows]
        if len(set(ids)) != len(ids) or not rows:
            raise ValueError(f'{path}: duplicate cluster IDs or empty cache')
        for r in rows:
            if any(not isinstance(r[k], int) or isinstance(r[k], bool) for k in ('size', 'H2', 'bonds', 'cluster_id')):
                raise ValueError(f'{path}: counts must be integers')
            if r['size'] < 1 or r['H2'] < 0 or r['bonds'] < 0 or r['cluster_id'] < 0 or 2*r['H2'] > r['bonds'] or r['bonds'] > 2*r['size']:
                raise ValueError(f'{path}: invalid counts')
        measured = (len(rows), sum(r['size'] for r in rows), max(r['size'] for r in rows), sum(r['H2'] for r in rows), sum(r['bonds'] for r in rows))
        saved = tuple(int(sample[k]) for k in ('clusters', 'N', 'N1', 'H2', 'bonds'))
        if measured != saved:
            raise ValueError(f'{path}: cluster data disagree with saved sample totals')
        records.append((path, rec))
    if len(records) != expected:
        raise ValueError(f'Expected {expected} matching trials, found {len(records)}; use the completed run cache, not pilot only.')
    cases = [r['sample']['case_id'] for _, r in records]
    seeds = [int(r['sample']['seed']) for _, r in records]
    if len(set(cases)) != expected or len(set(seeds)) != expected:
        raise ValueError('Duplicate cases or seeds')
    return records

def aggregate(records):
    # Every trial contributes to counts, including zeros. Fractions require a present bin.
    from collections import defaultdict
    by_trial = []
    keys = set()
    for _, rec in records:
        rows = rec['clusters']
        largest = min(rows, key=lambda r: (-r['size'], r['cluster_id']))['cluster_id']
        for scope in ('all', 'exclude_one_largest'):
            counts = defaultdict(lambda: [0, 0])
            for r in rows:
                if scope != 'all' and r['cluster_id'] == largest:
                    continue
                present = int(r['H2'] > 0)
                counts[('TOTAL', 0, 0)][present] += 1
                lo, hi = bounds(r['size'])
                counts[('BIN', lo, hi)][present] += 1
            by_trial.append((rec['sample']['case_id'], scope, dict(counts)))
            keys.update(counts)
    detail, summary = [], []
    for case, scope, counts in by_trial:
        for kind, lo, hi in sorted(keys):
            no, yes = counts.get((kind, lo, hi), (0, 0))
            detail.append(dict(case_id=case, scope=scope, kind=kind, size_low=lo, size_high=hi,
                               without_H2=no, with_H2=yes, cluster_count=no+yes,
                               fraction_with_H2=yes/(no+yes) if no+yes else None))
    for scope in ('all', 'exclude_one_largest'):
        for kind, lo, hi in sorted(keys):
            rows = [r for r in detail if r['scope']==scope and (r['kind'],r['size_low'],r['size_high'])==(kind,lo,hi)]
            out = dict(scope=scope, kind=kind, size_low=lo, size_high=hi)
            for metric in ('without_H2','with_H2','cluster_count','fraction_with_H2'):
                mean, se, n = stats(r[metric] for r in rows if r[metric] is not None)
                out.update({metric+'_mean':mean, metric+'_se':se, metric+'_n':n})
            assert math.isclose(out['without_H2_mean']+out['with_H2_mean'], out['cluster_count_mean'])
            summary.append(out)
    return detail, summary

def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def figures(summary, output, L, phi, n):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    colors = ['#8c99a5', '#258578']
    totals=[next(r for r in summary if r['scope']==s and r['kind']=='TOTAL') for s in ('all','exclude_one_largest')]
    fig, ax=plt.subplots(figsize=(8,5.6))
    x=np.arange(2)
    no=np.array([r['without_H2_mean'] for r in totals]);yes=np.array([r['with_H2_mean'] for r in totals])
    ax.bar(x,no,label='No H2',color=colors[0]);ax.bar(x,yes,bottom=no,label='At least one H2',color=colors[1])
    for i,r in enumerate(totals):
        ax.text(i,no[i]/2,f"{no[i]:,.2f}",ha='center',va='center',color='white')
        ax.text(i,no[i]+yes[i]/2,f"{yes[i]:,.2f}",ha='center',va='center',color='white')
    ax.set_xticks(x,['All clusters','Exclude one largest']);ax.set_ylabel('Mean clusters per trial')
    ax.set_title(f'H2 presence: L={L}, target phi={phi:g}, n={n} trials')
    ax.legend(loc='upper left',bbox_to_anchor=(1,1));ax.set_ylim(0,(no+yes).max()*1.12)
    fig.text(.08,.015,'Means of trial counts, including zeros. H2 = exactly-double molecular pair. SE in CSV.',fontsize=9)
    fig.tight_layout(rect=(0,.045,1,1))
    for ext in ('png','svg'):fig.savefig(output/f'01_H2_presence_totals.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    rows=[r for r in summary if r['scope']=='exclude_one_largest' and r['kind']=='BIN' and 2<=r['size_low'] and r['size_high']<=512]
    rows.sort(key=lambda r:r['size_low'])
    fig, axes=plt.subplots(2,1,figsize=(11,8),sharex=True)
    x=np.arange(len(rows));no=np.array([r['without_H2_mean'] for r in rows]);yes=np.array([r['with_H2_mean'] for r in rows])
    axes[0].bar(x,no,label='No H2',color=colors[0]);axes[0].bar(x,yes,bottom=no,label='At least one H2',color=colors[1])
    axes[0].set_ylabel('Mean clusters per trial');axes[0].legend()
    for i,r in enumerate(rows):
        if r['fraction_with_H2_n']>=50:
            axes[1].errorbar(i,100*r['fraction_with_H2_mean'],yerr=100*r['fraction_with_H2_se'] if r['fraction_with_H2_se'] is not None else None,fmt='o',color=colors[1],capsize=3)
        axes[1].text(i,103,str(r['fraction_with_H2_n']),ha='center',fontsize=9)
    axes[1].set_ylim(-3,111);axes[1].set_ylabel('Clusters with H2 [%]');axes[1].set_xlabel('Cluster size [molecules]')
    axes[1].set_xticks(x,[str(r['size_low']) if r['size_low']==r['size_high'] else f"{r['size_low']}-{r['size_high']}" for r in rows])
    fig.suptitle(f'H2 presence by size: L={L}, target phi={phi:g}; exclude one largest')
    fig.text(.08,.015,f'Counts: all {n} trials, absent bin = 0. Fractions: equal trial means among present bins.\nNumbers above points = contributing trials; fractions shown only for n>=50, mean +/- SE. Bins 2-512 shown; all bins in CSV.',fontsize=9)
    fig.tight_layout(rect=(0,.075,1,.97))
    for ext in ('png','svg'):fig.savefig(output/f'02_H2_presence_by_size.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,default=Path('data/validation/h2_relation_fine_compact/cycle46_cache')/CACHE_HASH)
    ap.add_argument('--L',type=int,default=128);ap.add_argument('--phi',type=float,default=.43)
    ap.add_argument('--samples',type=int,default=100)
    ap.add_argument('--output',type=Path,default=Path('data/validation/h2_relation_fine_compact/h2_presence_report'))
    args=ap.parse_args()
    if args.samples<2:raise ValueError('At least two trials required')
    records=load_cache(args.cache,args.L,args.phi,args.samples)
    detail, summary=aggregate(records)
    # Refuse overwriting an earlier result.
    output=args.output/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');output.mkdir(parents=True,exist_ok=False)
    write_csv(output/'sample_counts.csv',detail);write_csv(output/'presence_summary.csv',summary)
    figures(summary,output,args.L,args.phi,len(records))
    metadata=dict(status='completed',L=args.L,target_phi=args.phi,samples=len(records),
                  created_at_utc=datetime.now(timezone.utc).isoformat(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  inputs=[dict(case_id=r['sample']['case_id'],cache_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),raw_sha256=r.get('input_sha256')) for p,r in records],
                  definitions={'H2':'exactly two parallel bonds between a molecule pair',
                               'counts':'per trial counts; zeros included, all trials equal weight',
                               'fraction':'per trial positive cluster count / total bin cluster count, then equal contributing-trial mean',
                               'exclude':'remove one largest cluster; tied largest: lowest cluster_id',
                               'limits':'cached results; raw files not revalidated; H2 absent does not imply all cycles absent; no causal inference'})
    (output/'run_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    (output/'README_JA.md').write_text('目的：H2を1個以上含むクラスタと、H2を含まないクラスタがそれぞれ何個あるかを調べる。\nH2はちょうど2本の分子ペア。H2なしは閉路なしを意味しない。\n図1：全クラスタ／最大1個除外の1試料平均個数。\n図2上：サイズ別の平均個数。区間なしを0として全試料で平均。\n図2下：H2ありクラスタ割合。区間ありの試料のみで試料割合を等重み平均。50は表示上の目安。\n孤立分子は全体個数とCSVに含む。図2は2〜512分子のみ。全サイズはCSV。\n全体と区間別の結果を混同しない。全試料を合算した割合ではない。キャッシュの総数整合性を確認するが元の結合ファイルは再検査しない。\n',encoding='utf-8')
    zip_path=output/'h2_presence_share.zip'
    with ZipFile(zip_path,'w',ZIP_DEFLATED) as z:
        for p in output.iterdir():
            if p!=zip_path:z.write(p,p.name)
    for r in summary:
        if r['kind']=='TOTAL':print(r['scope'], 'H2あり:',round(r['with_H2_mean'],3),'H2なし:',round(r['without_H2_mean'],3))
    print('OUTPUT:',output.resolve());print('SHARE:',zip_path.resolve())

if __name__=='__main__':main()
