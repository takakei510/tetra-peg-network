#!/usr/bin/env python3
"""Bounded DFS of distinct-molecule 4/6 cycles in saved AB bond graphs.

No network regeneration. Counts all simple cycles, including cycles with chords.
Counts molecular routes once under rotation/reversal and separately counts
edge-distinct choices. Endpoint fractions use UNION of participating bonds.
"""
import argparse
from collections import Counter, defaultdict, deque
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import time
from zipfile import ZipFile, ZIP_DEFLATED

REPO = Path(__file__).resolve().parents[2]
VERSION = 'cycle46-v1'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def csv_rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        yield from csv.DictReader(f)


def raw_path(folder, name):
    paths = [p for p in (folder/name, folder/(name+'.gz')) if p.is_file()]
    if len(paths) != 1:
        raise ValueError(f'{folder}: {name}(.gz) must exist in exactly one copy')
    return paths[0]


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    tmp.replace(path)


def write_csv(path, rows):
    if not rows: raise ValueError(f'Empty table: {path}')
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def graph(n, bonds):
    """Validate physical endpoint ownership and reconstruct components by BFS."""
    if n < 1: raise ValueError('N must be positive')
    pairs = Counter(); endpoints = set(); bids = set(); types = {}
    for r in bonds:
        a, b = int(r['molecule_a']), int(r['molecule_b'])
        ea, eb, bid = int(r['endpoint_a']), int(r['endpoint_b']), int(r['bond_id'])
        if not (0 <= a < n and 0 <= b < n and a != b and bid >= 0):
            raise ValueError('Invalid molecule/bond ID or self-loop')
        if ea//4 != a or eb//4 != b or min(ea,eb) < 0 or ea in endpoints or eb in endpoints or ea == eb or bid in bids:
            raise ValueError('Endpoint owner/reuse or duplicate bond ID')
        for m, key in ((a, 'type_a'), (b, 'type_b')):
            if key in r:
                t = r[key]
                if t not in ('A','B') or (m in types and types[m] != t): raise ValueError('Inconsistent AB type')
                types[m] = t
        if a in types and b in types and types[a] == types[b]: raise ValueError('Non-AB bond')
        endpoints.update((ea, eb)); bids.add(bid)
        pairs[tuple(sorted((a,b)))] += 1
    adjacent = [set() for _ in range(n)]; degrees = [0]*n
    for (a,b), k in pairs.items():
        adjacent[a].add(b); adjacent[b].add(a); degrees[a]+=k; degrees[b]+=k
    if max(degrees) > 4: raise ValueError('More than four bonded arms')
    owner = [-1]*n; color = [-1]*n; clusters=[]
    for start in range(n):
        if owner[start] >= 0: continue
        cid=len(clusters); queue=deque([start]); owner[start]=cid; color[start]=0; size=0
        while queue:
            a=queue.popleft(); size+=1
            for b in adjacent[a]:
                if owner[b] < 0:
                    owner[b]=cid; color[b]=1-color[a]; queue.append(b)
                elif color[b] == color[a]: raise ValueError('Graph is not bipartite')
        clusters.append(dict(cluster_id=cid,size=size,bonds=0,H2=0,C4_mol=0,C6_mol=0,C4_edge=0,C6_edge=0,B4=0,B6=0))
    for (a,b), k in pairs.items():
        row=clusters[owner[a]]; row['bonds']+=k; row['H2']+=k==2
    return adjacent,pairs,owner,clusters,degrees


def count_cycles(n, bonds):
    adj,pairs,owner,rows,degrees=graph(n,bonds)
    # Leaf removal on the SIMPLE projection preserves all simple cycles.
    live=set(range(n)); degree=[len(a) for a in adj]
    queue=deque(a for a in range(n) if degree[a]<2)
    while queue:
        a=queue.popleft()
        if a not in live: continue
        live.remove(a)
        for b in adj[a]:
            if b in live:
                degree[b]-=1
                if degree[b]==1:queue.append(b)
    core={a:tuple(sorted(adj[a]&live)) for a in live}
    participating={4:set(),6:set()}
    for start in sorted(live):
        path=[start]; visited={start}
        def dfs(a):
            depth=len(path)
            for b in core[a]:
                if b==start:
                    if depth in (4,6) and path[1]<path[-1]:
                        # start is the smallest vertex; comparison removes reversal.
                        edges=[tuple(sorted((path[i],path[(i+1)%depth]))) for i in range(depth)]
                        row=rows[owner[start]]
                        row[f'C{depth}_mol']+=1
                        row[f'C{depth}_edge']+=math.prod(pairs[e] for e in edges)
                        participating[depth].update(edges)
                elif depth<6 and b>start and b not in visited:
                    visited.add(b);path.append(b);dfs(b);path.pop();visited.remove(b)
        dfs(start)
    for k, edges in participating.items():
        for a,b in edges: rows[owner[a]][f'B{k}']+=pairs[(a,b)]
    largest=min(rows,key=lambda r:(-r['size'],r['cluster_id']))['cluster_id']
    for r in rows:
        s=r['size'];r['is_largest']=int(r['cluster_id']==largest)
        r['qH2']=r['H2']/s
        for k in (4,6):
            r[f'q{k}']=r[f'B{k}']/(2*s)
            r[f'C{k}_per_molecule']=r[f'C{k}_mol']/s
            if not 0<=r[f'q{k}']<=1:raise ValueError('Invalid union endpoint fraction')
    return rows,owner,degrees


def size_bin(s):
    hi=1<<(s-1).bit_length()
    return (1 if hi==1 else hi//2+1),hi


def summarize_case(case, rows):
    result=[]
    for scope in ('all','exclude_one_largest'):
        bins=defaultdict(list)
        for r in rows:
            if scope=='exclude_one_largest' and r['is_largest']:continue
            bins[size_bin(r['size'])].append(r)
        for (lo,hi),rr in sorted(bins.items()):
            item=dict(case_id=case['case_id'],L=int(case['L']),target_phi=float(case['target_phi']),scope=scope,size_low=lo,size_high=hi,cluster_count=len(rr))
            for m in ('qH2','q4','q6','C4_per_molecule','C6_per_molecule'):
                item[m]=statistics.mean(r[m] for r in rr)
            result.append(item)
    return result


def inspect_case(study, sample):
    case_id=sample['case_id']
    if Path(case_id).name!=case_id:raise ValueError('Invalid case ID')
    suite=study/'cases'/case_id
    folder=suite/'runs'/case_id
    bp=raw_path(folder,'bonds.csv');mp=raw_path(folder,'cluster-members.csv')
    hashes={'bonds':digest(bp),'members':digest(mp)}
    t=time.perf_counter();bonds=list(csv_rows(bp))
    n=int(sample['N']); rows,owner,degrees=count_cycles(n,bonds)
    if (len(bonds),len(rows),max(r['size'] for r in rows),sum(r['H2'] for r in rows)) != (int(sample['bonds']),int(sample['clusters']),int(sample['N1']),int(sample['H2'])):
        raise ValueError(f'{case_id}: graph totals disagree with samples.csv')
    seen=set()
    for r in csv_rows(mp):
        m=int(r['molecule'])
        if m in seen or not 0<=m<n:raise ValueError('Duplicate/invalid membership')
        seen.add(m)
        if (int(r['cluster_id']),int(r['cluster_size']),int(r['bond_degree'])) != (owner[m],rows[owner[m]]['size'],degrees[m]):raise ValueError('BFS differs from saved cluster membership')
    if len(seen)!=n:raise ValueError('Missing molecules in membership')
    return dict(version=VERSION,sample=sample,input_sha256=hashes,elapsed_seconds=time.perf_counter()-t,clusters=rows,bins=summarize_case(sample,rows))


def aggregate(records, expected):
    groups=defaultdict(list)
    for rec in records:
        for r in rec['bins']:groups[(r['scope'],r['size_low'],r['size_high'])].append(r)
    out=[]
    for (scope,lo,hi), rr in sorted(groups.items()):
        n=len(rr); item=dict(scope=scope,size_low=lo,size_high=hi,n_contributing=n,n_total=expected,
                           cluster_count_mean=sum(r['cluster_count'] for r in rr)/expected)
        for m in ('qH2','q4','q6','C4_per_molecule','C6_per_molecule'):
            vv=[r[m] for r in rr];item[m+'_mean']=statistics.mean(vv)
            item[m+'_se']=statistics.stdev(vv)/math.sqrt(n) if n>1 else None
        out.append(item)
    return out


def draw(output, summary, config):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[r for r in summary if r['scope']==config['scope']]
    labels=[str(r['size_low']) if r['size_low']==r['size_high'] else f"{r['size_low']}-{r['size_high']}" for r in rows]
    x=list(range(len(rows)));colors=['#2873B8','#BD681B','#218477']
    definitions=[('01_endpoint_participation',('qH2','q4','q6'),('H2: exactly two parallel bonds','4-molecule cycles','6-molecule cycles'),'Participating endpoints / all endpoints [%]',100),
                 ('02_cycle_counts',('C4_per_molecule','C6_per_molecule'),('4-molecule cycles','6-molecule cycles'),'Molecular cycles per 100 molecules',100)]
    for name, metrics, legends, ylabel, factor in definitions:
        fig,ax=plt.subplots(figsize=(12,6))
        for j,(m,label) in enumerate(zip(metrics,legends)):
            xx=[v+(j-(len(metrics)-1)/2)*.12 for v in x]
            yy=[factor*r[m+'_mean'] for r in rows];err=[factor*(r[m+'_se'] or 0) for r in rows]
            ax.errorbar(xx,yy,yerr=err,fmt='none',ecolor=colors[j],capsize=2)
            for adequate in (True,False):
                mask=[i for i,r in enumerate(rows) if (r['n_contributing']>=config['display_min_samples'])==adequate]
                ax.scatter([xx[i] for i in mask],[yy[i] for i in mask],color=colors[j],facecolors=colors[j] if adequate else 'none',label=label if adequate else None,s=30)
        ax.set_xticks(x,labels,rotation=45,ha='right');ax.set_xlabel('Cluster size bin [molecules]');ax.set_ylabel(ylabel);ax.set_ylim(bottom=0)
        ax.set_title(f"L={config['L']}, target phi={config['target_phi']:.2f}; one largest cluster excluded")
        ax.grid(axis='y',alpha=.2);ax.legend(fontsize=9)
        fig.text(.1,.015,'Mean +/- SE. Hollow: fewer than '+str(config['display_min_samples'])+' contributing samples. Shared endpoints counted once per length.',fontsize=9)
        fig.tight_layout(rect=(0,.04,1,1))
        for ext in ('png','svg'):fig.savefig(output/f'{name}.{ext}',dpi=180)
        plt.close(fig)
    fig,(a,b)=plt.subplots(2,1,figsize=(12,6),sharex=True)
    a.bar(x,[r['n_contributing'] for r in rows],color='#2873B8');a.axhline(config['display_min_samples'],color='gray',linestyle='--');a.set_ylabel('Contributing samples');a.set_ylim(0,config['samples']*1.05)
    b.bar(x,[r['cluster_count_mean'] for r in rows],color='#218477');b.set_yscale('log');b.set_ylabel('Clusters per sample');b.set_xticks(x,labels,rotation=45,ha='right');b.set_xlabel('Cluster size bin [molecules]')
    fig.suptitle('Support for size-conditioned cycle statistics');fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(output/f'03_support.{ext}',dpi=180)
    plt.close(fig)


def choose_samples(study, config, explicit):
    if explicit:
        source=explicit
    else:
        options=[]
        for f in study.glob('results_*'):
            if not (f/'manifest.json').is_file() or (f/'EXPORT_INCOMPLETE.json').exists():continue
            m=json.loads((f/'manifest.json').read_text())
            if m.get('status')=='completed' and m.get('completed_samples')==m.get('requested_samples'):
                options.append((m['completed_samples'],f.stat().st_mtime,f/'samples.csv'))
        if not options:raise ValueError('No completed results_* export; specify --samples-csv')
        source=max(options)[2]
    selected=[r for r in csv_rows(source) if int(r['L'])==config['L'] and math.isclose(float(r['target_phi']),config['target_phi'],abs_tol=1e-9)]
    selected.sort(key=lambda r:r['case_id'])
    if len(selected)!=config['samples'] or len({r['seed'] for r in selected})!=len(selected) or len({r['case_id'] for r in selected})!=len(selected):raise ValueError('Expected exactly 100 distinct cases/seeds for selected condition')
    # Preflight all required raw files before any expensive DFS.
    for r in selected:
        folder=study/'cases'/r['case_id']/'runs'/r['case_id']
        raw_path(folder,'bonds.csv');raw_path(folder,'cluster-members.csv')
    return source,selected


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=REPO/'configs/analysis/cycle-46.json')
    ap.add_argument('--study',type=Path)
    ap.add_argument('--samples-csv',type=Path)
    ap.add_argument('--mode',choices=['pilot','run'],default='pilot')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args();config=json.loads(args.config.read_text())
    if config['schema_version']!=1 or config['lengths']!=[4,6] or config['scope']!='exclude_one_largest' or config['samples']!=100:raise ValueError('Unsupported configuration')
    study=(args.study or REPO/config['study']).resolve()
    source,cases=choose_samples(study,config,args.samples_csv)
    code_hash=digest(Path(__file__))
    checkpoint=study/'cycle46_cache'/code_hash[:16]
    records=[];start=time.perf_counter()
    chosen=cases[:1] if args.mode=='pilot' else cases
    for i,sample in enumerate(chosen,1):
        cp=checkpoint/(sample['case_id']+'.json');rec=None
        if cp.is_file():
            old=json.loads(cp.read_text())
            folder=study/'cases'/sample['case_id']/'runs'/sample['case_id']
            hashes={'bonds':digest(raw_path(folder,'bonds.csv')),'members':digest(raw_path(folder,'cluster-members.csv'))}
            if old.get('sample')==sample and old.get('input_sha256')==hashes and old.get('version')==VERSION:rec=old
        if rec is None:
            rec=inspect_case(study,sample);write_json(cp,rec)
        records.append(rec)
        totals={f'C{k}':sum(r[f'C{k}_mol'] for r in rec['clusters']) for k in (4,6)}
        print(f'{i}/{len(chosen)} {sample["case_id"]}: {totals}, {rec["elapsed_seconds"]:.2f}s',flush=True)
    if args.mode=='pilot':
        print(f'PILOT ONLY. No 100-sample conclusions. Raw graph verified. Cached at {checkpoint}')
        print('Next: run the same command with --mode run. Pilot timing is not a guaranteed runtime for all cases.')
        return
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output=args.output or study/'cycle46_outputs'/stamp;output.mkdir(parents=True,exist_ok=False)
    bins=[r for rec in records for r in rec['bins']];summary=aggregate(records,config['samples'])
    write_csv(output/'sample_bins.csv',bins);write_csv(output/'size_summary.csv',summary)
    totals=[]
    for rec in records:
        item={k:rec['sample'][k] for k in ('case_id','L','target_phi','seed','N','bonds','clusters','H2')}
        for k in (4,6):
            for m in (f'C{k}_mol',f'C{k}_edge',f'B{k}'):item[m]=sum(r[m] for r in rec['clusters'])
        item['elapsed_seconds']=rec['elapsed_seconds'];totals.append(item)
    write_csv(output/'sample_totals.csv',totals)
    write_json(output/'run_metadata.json',dict(status='completed',created_at_utc=datetime.now(timezone.utc).isoformat(),config=config,mode=args.mode,version=VERSION,script_sha256=code_hash,samples_sha256=digest(source),samples=len(records),wall_seconds=time.perf_counter()-start,raw_inputs=[dict(case_id=r['sample']['case_id'],sha256=r['input_sha256']) for r in records],definitions={'Ck_mol':'simple distinct-molecule cycles modulo rotation/reversal; chords allowed; parallel pairs collapsed','Ck_edge':'sum of products of edge multiplicities over molecular cycles','qk':'2 * union of original participating bonds / (4 * cluster size)','qH2':'exactly-double pair count / cluster size','aggregation':'within-sample equal cluster weights, then equal contributing sample weights; missing bin excluded, existing zero included','warning':'q4 and q6 can overlap with each other and qH2; do not add. Counts are not independent cycle rank.'}))
    guide='''# 4分子・6分子閉路の図の目的

目的：H2の2分子間閉路に加え、4・6個の異なる分子を経由する閉路が、どのサイズのクラスタに現れるかを調べる。L128、目標φ0.43の同じ100試料を使う。濃度による一般化や成長抑制の因果は、この1条件では検証しない。

01_endpoint_participation：閉路に参加する末端が全末端の何%か。4/6分子閉路で使う成立結合の和集合を取り、共有結合・末端を同じ長さ内で重複計上しない。各割合の分母は4s。H2はちょうど2本のペアのみ。長さ間には重複があるので足し合わせない。4C4/(2s)や6C6/(2s)は重複を含むため、この割合には使わない。

02_cycle_counts：クラスタの分子100個あたり、異なる分子経路として何個の閉路があるか。これは%ではなく個数密度。多重辺の選択は主図で区別せず、辺の選択まで区別したC4_edge/C6_edgeはsample_totals.csvに別出力する。独立閉路数E-N+Kとも異なる。

03_support：どのサイズ区間を何試料で比較できるか。上段は区間が存在した試料数、下段は全100試料（区間なしは0個）で平均したクラスタ数。n>=50は表示上の目安で、精度の保証ではない。

全図は最大クラスタを試料ごとに1個除く。同率最大はcluster_id最小を除く。割合と個数密度は試料内クラスタ等重み平均を求めた後、その区間がある試料間で等重み平均。存在するが閉路0のクラスタも含む。区間自体がない試料は割合平均から除く。誤差棒はSE=標本SD/√n。異なる区間は寄与試料が異なり得る。

4分子閉路はs<4、6分子閉路はs<6で構造的に0。この制約だけによる違いを物理的な発見として扱わない。5〜8区間には6分子閉路が不可能な5分子クラスタが含まれる。グラフは100試料の静的な比較で、同一クラスタの時間発展ではない。

DFSは長さ6で打ち切り。最小分子IDを開始点とし、逆向きの重複を除く。弦を持つ単純閉路も含み、最小環だけを数える定義ではない。小さい検証グラフの結果と実試料を混同しない。
'''
    (output/'figure_guide.md').write_text(guide,encoding='utf-8')
    try:draw(output,summary,config)
    except ImportError:
        print('CSVs saved, but matplotlib is missing. Install matplotlib and rerun (cached DFS will be reused).',file=sys.stderr)
        return
    archive=output/'cycle46_share.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for f in sorted(output.iterdir()):
            if f.is_file() and f!=archive:z.write(f,f.name)
    with ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Corrupt ZIP')
    print(f'COMPLETE: {archive}',flush=True)


if __name__=='__main__':
    main()
