#!/usr/bin/env python3
"""Plot completed H2 CSV exports; never runs or modifies a simulation."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import uuid
from zipfile import ZipFile, ZIP_DEFLATED

REPO = Path(__file__).resolve().parents[2]
REQUIRED = ('summary.csv','size_relation.csv','convergence.csv','manifest.json')


def resolve_results(source):
    if (source/'manifest.json').is_file():
        candidates=[source]
    else:
        candidates=list(source.glob('results_*'))
    valid=[]
    for path in candidates:
        if (path/'EXPORT_INCOMPLETE.json').exists() or not all((path/name).is_file() for name in REQUIRED):
            continue
        m=json.loads((path/'manifest.json').read_text(encoding='utf-8'))
        if m.get('status')=='completed' and m.get('completed_samples')==m.get('requested_samples') and m.get('completed_samples',0)>0:
            valid.append((m['completed_samples'],path.stat().st_mtime_ns,path))
    if not valid:
        raise ValueError('No completed result export found. Set --results to the completed results_* directory.')
    # Full study wins over earlier capacity/pilot, then newest complete export.
    return max(valid)[2]


def read(name):
    with (ROOT/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)


def concentration(ax,key,label,factor=1,xlim=None):
    for L in colors:
        g=sorted([r for r in s if int(r['L'])==L],key=lambda r:float(r['actual_phi']))
        ax.errorbar([float(r['actual_phi']) for r in g],[factor*float(r[key+'_mean']) for r in g],
            yerr=[factor*float(r[key+'_se']) if r[key+'_se'] else 0 for r in g],fmt='o',markersize=3,capsize=2,linestyle='none',color=colors[L],label=f'L={L}')
    ax.set(xlabel='Actual occupied fraction Φ',ylabel=label)
    if xlim:ax.set_xlim(*xlim)
    ax.grid(alpha=.15);ax.legend()


def main():
    global ROOT,OUT,s,rel,conv,colors,plt
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=REPO/'data/validation/h2_relation_fine_compact',help='Completed results_* folder or study folder (select most complete export)')
    p.add_argument('--output',type=Path,help='New output directory; existing directories are not overwritten')
    p.add_argument('--phi',type=float,default=.44,help='Target phi for cluster-size relation')
    p.add_argument('--size',type=int,default=128,help='L for the relation + seed count figure')
    args=p.parse_args()
    ROOT=resolve_results(args.results.resolve())
    manifest=json.loads((ROOT/'manifest.json').read_text())
    s,rel,conv=read('summary.csv'),read('size_relation.csv'),read('convergence.csv')
    if not s or not rel or not conv:raise ValueError('Empty CSV export')
    if len({(int(r['L']),float(r['target_phi'])) for r in s})!=len(s):raise ValueError('Duplicate summary condition')
    if sum(int(r['n_samples']) for r in s)!=manifest['completed_samples']:raise ValueError('Summary sample count disagrees with manifest')
    sizes=sorted({int(r['L']) for r in s})
    if len(sizes)>3:raise ValueError('This comparison layout supports at most three L values')
    phi=args.phi;selected_size=args.size
    if not math.isfinite(phi) or selected_size not in sizes:raise ValueError('Invalid --phi or --size')
    for L in sizes:
        if not any(int(r['L'])==L and abs(float(r['target_phi'])-phi)<1e-8 and r['scope']=='all' for r in rel):
            raise ValueError(f'No size relation for L={L}, target phi={phi}; choose --phi from your experiment')
    counts=sorted({int(r['n_samples']) for r in s})
    n_label=str(counts[0]) if len(counts)==1 else f'{counts[0]}..{counts[-1]}'
    max_n=max(counts);threshold=max(1,math.ceil(max_n/2))
    peaks=[]
    for L in sizes:
        peak=max((r for r in s if int(r['L'])==L),key=lambda r:float(r['cluster_weighted_excl_one_mean']))
        peaks.append((L,float(peak['target_phi'])))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ValueError('Matplotlib is missing. Ubuntu/WSL: sudo apt install python3-matplotlib') from exc
    OUT=args.output or REPO/'data/figures'/('h2_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:6])
    OUT=OUT.resolve()
    OUT.mkdir(parents=True,exist_ok=False)
    colors={L:c for L,c in zip(sizes,['#306aa0','#dd8630','#25846f'])}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(15,4.4),layout='constrained')
    for ax,key,label,factor in zip(axs,['H2_per_volume','H2_endpoint_fraction_all','H2_endpoint_fraction_bonded'],['H₂ per 1,000 lattice sites','All endpoints in H₂ pairs (%)','Bonded endpoints in H₂ pairs (%)'],[1000,100,100]):concentration(ax,key,label,factor)
    fig.suptitle(f'H₂ measures: n={n_label} independent seeds per condition; mean ± SE; no fit/interpolation')
    save(fig,'01_h2_measures')
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    concentration(axs[0],'cluster_weighted_excl_one','Weighted size Sfin (molecules)');concentration(axs[1],'cluster_weighted_excl_one','Weighted size Sfin (molecules)',xlim=(.37,.53))
    fig.suptitle('One largest cluster excluded per sample; mean ± SE (not SD)')
    save(fig,'02_sfin')
    fig,axs=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    for ax,key,label in zip(axs,['clusters','K_per_volume','S1'],['Cluster count K (including isolated molecules)','Cluster count per lattice site K/L³','Largest cluster occupancy 21N₁/L³']):concentration(ax,key,label)
    fig.suptitle('Cluster count and largest-cluster occupancy: mean ± SE')
    save(fig,'03_cluster_counts_and_largest')
    fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for ax,scope in zip(axs,['all','exclude_one_largest']):
     for L in colors:
      g=[r for r in rel if int(r['L'])==L and abs(float(r['target_phi'])-phi)<1e-8 and r['scope']==scope]
      for valid in [True,False]:
       rows=[r for r in g if (int(r['samples_with_bin'])>=threshold)==valid]
       ax.errorbar([math.sqrt(int(r['size_low'])*int(r['size_high'])) for r in rows],
        [100*float(r['fraction_mean']) for r in rows],yerr=[100*float(r['fraction_se']) if r['fraction_se'] else 0 for r in rows],
        fmt='o',markersize=4,capsize=2,linestyle='none',color=colors[L],mfc=colors[L] if valid else 'white',label=f'L={L}'+('' if valid else f', n<{threshold}'))
     ax.set(xscale='log',xlabel='Cluster size bin (geometric midpoint)',ylabel='Mean h₂,i/sᵢ (%)',title=scope.replace('_',' '));ax.grid(alpha=.15);ax.legend(fontsize=8)
    fig.suptitle(f'Target Φ={phi:.2f}: cluster-equal mean within each seed, then seed-equal mean ± SE')
    save(fig,'04_size_relation')
    fig,axs=plt.subplots(2,1,figsize=(9,7),layout='constrained',sharex=True,gridspec_kw={'height_ratios':[2,1]})
    g=[r for r in rel if int(r['L'])==selected_size and abs(float(r['target_phi'])-phi)<1e-8 and r['scope']=='all']
    x=[math.sqrt(int(r['size_low'])*int(r['size_high'])) for r in g]
    axs[0].errorbar(x,[100*float(r['fraction_mean']) for r in g],yerr=[100*float(r['fraction_se']) if r['fraction_se'] else 0 for r in g],fmt='o',capsize=3,color=colors[selected_size],linestyle='none')
    axs[0].set(ylabel='Mean H₂ endpoint fraction (%)',title=f'L={selected_size}, target Φ={phi:.2f}: all clusters, mean ± SE')
    axs[1].scatter(x,[int(r['samples_with_bin']) for r in g],color=colors[selected_size]);axs[1].set(ylim=(0,max_n*1.05),ylabel=f'Contributing seeds / {max_n}',xlabel='Cluster size bin (geometric midpoint)',xscale='log')
    for ax in axs:ax.grid(alpha=.15)
    save(fig,'05_relation_with_sample_counts')
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for L,peak_phi in peaks:
     for ax,metric,label in zip(axs,['H2_endpoint_fraction_all','cluster_weighted_excl_one'],['All-endpoint fraction in H₂ pairs','Sfin (molecules)']):
      g=sorted([r for r in conv if int(r['L'])==L and abs(float(r['target_phi'])-peak_phi)<1e-8 and r['metric']==metric],key=lambda r:int(r['first_n']))
      ax.errorbar([int(r['first_n']) for r in g],[float(r['mean']) for r in g],yerr=[float(r['se']) if r['se'] else 0 for r in g],fmt='o',capsize=3,linestyle='none',color=colors[L],label=f'L={L}, Φ={peak_phi:.2f}')
      ax.set(xlabel='First n seeds (nested subsets)',ylabel=label);ax.grid(alpha=.15);ax.legend(fontsize=9)
    fig.suptitle('Mean ± SE versus sample count; nested subsets are not independent')
    save(fig,'06_convergence')

    input_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in REQUIRED}
    metadata=dict(input_directory=str(ROOT),input_sha256=input_hashes,samples=manifest['completed_samples'],
        conditions=len(s),relation_target_phi=phi,relation_size=selected_size,
        error_bars='SE; per-seed means, not SD or confidence intervals',
        rendering='points only; no fitted/interpolated/connecting lines',
        convergence='nested seed subsets; selected per-L observed Sfin maximum, not true critical point',
        matplotlib=matplotlib.__version__,plot_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'plot_manifest.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (OUT/'READ_ME_JA.md').write_text("""# PCで作成したH₂検証図

PNGは閲覧・スライド用、SVGは拡大・編集用。誤差棒は標準誤差SEです。有効試料が1個でSEが未定義なら、点だけ表示します。
線・fit・補間は使っていません。背景の格子線は読み取り補助です。

- 01：H₂密度、全末端割合H₂/N、結合済み末端割合2H₂/E。横軸は実占有率。
- 02：最大クラスタを1個除いた重み付き平均Sfin。ピークは転移点の確定ではありません。
- 03：生のクラスタ数K、K/L³、最大クラスタ占有率21N1/L³。
- 04：指定ΦでのクラスタサイズとH₂末端割合。左は全クラスタ、右は最大1個除外。白抜きは有効seed数が全試料の半数未満。
- 05：指定L・ΦのクラスタサイズとH₂末端割合。下段はサイズ区間に寄与したseed数。
- 06：先頭n試料での平均とSE。入れ子の部分集合で、独立の実験ではありません。条件は各Lの観測Sfin最大点。

04/05の横軸は固定サイズ区間の幾何中心、対数表示。縦軸は100h2_i/s_i (%)。
まず各試料内でクラスタを等重み平均し、次に寄与seedを等重み平均。
H₂=0のクラスタも含めます。空区間は0に置き換えません。

サイズと閉路割合の関係だけから、閉路が成長を妨げたという因果までは言えません。
元のCSV・結合データは変更していません。入力ハッシュ・条件はplot_manifest.jsonを参照。
次の解析やスライド構成に進む前に、目的・定義・計算・観測・限界の5点を確認しましょう。
""",encoding='utf-8')
    archive=OUT/'upload_h2_plots.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for path in sorted(OUT.iterdir()):
            if path!=archive:z.write(path,path.name)
    with ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Plot ZIP CRC mismatch')
    print(f'Created 6 figures (PNG + SVG). No simulations were started.\nInput: {ROOT}\nOutput: {OUT}\nUPLOAD THIS FILE:\n{archive}')


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError) as exc:
        sys.exit(f'ERROR: {exc}')
