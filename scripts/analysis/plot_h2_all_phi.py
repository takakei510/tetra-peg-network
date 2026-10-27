#!/usr/bin/env python3
"""Read completed size_relation.csv and compare ALL measured concentrations (no simulation)."""
import argparse
import csv
import json
import hashlib
import math
from pathlib import Path
import sys
from datetime import datetime,timezone
from zipfile import ZipFile,ZIP_DEFLATED


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path(__file__).resolve().parents[2]/'data/validation/h2_relation_fine_compact')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    from plot_h2_validation_results import resolve_results
    source=resolve_results(args.results.resolve())
    with (source/'size_relation.csv').open(encoding='utf-8-sig',newline='') as f:raw=list(csv.DictReader(f))
    data={}
    for r in raw:
        key=int(r['L']),round(float(r['target_phi'])*100),r['scope'],int(r['size_low']),int(r['size_high'])
        if key in data:raise ValueError('Duplicate size-relation row')
        n,total=int(r['samples_with_bin']),int(r['total_samples'])
        if not 1<=n<=total:raise ValueError('Invalid contributing seed count')
        data[key]=dict(mean=100*float(r['fraction_mean']),se=100*float(r['fraction_se']) if r['fraction_se'] else None,n=n,total=total,clusters=int(r['clusters']))
    sizes=sorted({k[0] for k in data});phis=sorted({k[1] for k in data});bins=sorted({(k[3],k[4]) for k in data})
    if sizes!=[32,64,128] or phis!=list(range(10,71)):raise ValueError('This layout expects L=32/64/128, target phi=0.10..0.70 at 0.01 steps')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    out=(args.output or source.parent.parent.parent/'figures'/('h2_all_phi_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))).resolve()
    out.mkdir(parents=True,exist_ok=False)
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(out/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    colors=['#306aa0','#dd8630','#25846f'];scopes=['all','exclude_one_largest']
    # Same color limits and size bins across all six panels. Nearest-neighbor cells, no smoothing.
    vmax=math.ceil(max(v['mean'] for v in data.values() if v['n']>=math.ceil(v['total']/2)))
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#dddddd')
    fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained',sharex=True,sharey=True)
    for row,scope in enumerate(scopes):
        for col,L in enumerate(sizes):
            ax=axes[row,col];values=np.full((len(bins),len(phis)),np.nan);rare=[]
            for i,(lo,hi) in enumerate(bins):
                for j,p in enumerate(phis):
                    v=data.get((L,p,scope,lo,hi))
                    if v and v['n']>=math.ceil(v['total']/2):values[i,j]=v['mean']
                    elif v:rare.append((p/100,i))
            im=ax.imshow(values,origin='lower',aspect='auto',interpolation='nearest',extent=(.095,.705,-.5,len(bins)-.5),cmap=cmap,vmin=0,vmax=vmax)
            if rare:ax.scatter(*zip(*rare),marker='.',s=4,c='#555555')
            ax.set_title(f'L={L}; '+scope.replace('_',' '));ax.set_xticks([.1,.2,.3,.4,.5,.6,.7])
            ax.set_yticks(range(len(bins)),[str(a) if a==b else f'{a}–{b}' for a,b in bins],fontsize=8)
            if row==1:ax.set_xlabel('Target occupied fraction Φ')
            if col==0:ax.set_ylabel('Cluster size bin (molecules)')
    fig.colorbar(im,ax=axes,label='Mean H₂ endpoint fraction (%)',shrink=.8)
    fig.suptitle('All 61 concentrations: gray = absent or <50% of seeds; dot on gray = present but rare\nSame size bins/color scale; cluster mean within seed, then seed mean; no causal inference',fontsize=13)
    save(fig,'07_all_phi_fraction')
    fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained',sharex=True,sharey=True)
    for row,scope in enumerate(scopes):
        for col,L in enumerate(sizes):
            ax=axes[row,col];counts=np.zeros((len(bins),len(phis)))
            for i,(lo,hi) in enumerate(bins):
                for j,p in enumerate(phis):
                    v=data.get((L,p,scope,lo,hi))
                    if v:counts[i,j]=100*v['n']/v['total']
            im=ax.imshow(counts,origin='lower',aspect='auto',interpolation='nearest',extent=(.095,.705,-.5,len(bins)-.5),cmap='cividis',vmin=0,vmax=100)
            ax.set_title(f'L={L}; '+scope.replace('_',' '));ax.set_xticks([.1,.2,.3,.4,.5,.6,.7])
            ax.set_yticks(range(len(bins)),[str(a) if a==b else f'{a}–{b}' for a,b in bins],fontsize=8)
            if row==1:ax.set_xlabel('Target occupied fraction Φ')
            if col==0:ax.set_ylabel('Cluster size bin (molecules)')
    fig.colorbar(im,ax=axes,label='Seeds with >=1 cluster in bin (%)',shrink=.8)
    fig.suptitle('Size-bin availability across all concentrations (NOT a cluster-count histogram)',fontsize=13)
    save(fig,'08_all_phi_availability')
    # These concentrations are selected by regime, not by maximizing an observed effect.
    representatives=[20,44,70]
    fig,axes=plt.subplots(2,3,figsize=(16,8),layout='constrained',sharex=True,sharey='row')
    for col,p in enumerate(representatives):
        for row,scope in enumerate(scopes):
            ax=axes[row,col]
            for L,color in zip(sizes,colors):
                for frequent in [True]:
                    found=[(lo,hi,v) for (l,q,sc,lo,hi),v in data.items() if (l,q,sc)==(L,p,scope) and (v['n']>=math.ceil(v['total']/2))==frequent]
                    ax.errorbar([math.sqrt(lo*hi) for lo,hi,v in found],[v['mean'] for lo,hi,v in found],
                        yerr=[v['se'] or 0 for lo,hi,v in found],fmt='o',markersize=3.5,capsize=2,linestyle='none',color=color,mfc=color if frequent else 'white',label=f'L={L}'+('' if frequent else ', rare'))
            ax.set(xscale='log',xlim=(.85,8e4),ylim=(-.5,15),title=f'Φ={p/100:.2f}; '+scope.replace('_',' '));ax.grid(alpha=.15)
            if row==1:ax.set_xlabel('Cluster size bin (geometric midpoint)')
            if col==0:ax.set_ylabel('Mean H₂ endpoint fraction (%)')
            if col==2:ax.legend(fontsize=7)
    fig.suptitle('Low density (0.20), previous reference (0.44), high density (0.70): mean ± SE\nOnly bins with >=50 contributing seeds; rare/absent bins omitted; zero-H₂ clusters included',fontsize=13)
    save(fig,'09_density_comparison')
    fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout='constrained',sharey=True)
    bin_colors=['#306aa0','#dd8630','#25846f','#9d4e98']
    selected_bins=[(2,2),(5,8),(17,32),(65,128)]
    for ax,L in zip(axes,sizes):
        for (lo,hi),color in zip(selected_bins,bin_colors):
            rows=[(p,data[L,p,'all',lo,hi]) for p in phis if (L,p,'all',lo,hi) in data and data[L,p,'all',lo,hi]['n']>=50]
            ax.errorbar([p/100 for p,v in rows],[v['mean'] for p,v in rows],yerr=[v['se'] or 0 for p,v in rows],fmt='o',markersize=2.5,capsize=1,linestyle='none',color=color,label=f'{lo}' if lo==hi else f'{lo}–{hi}')
        ax.set(title=f'L={L}',xlabel='Target Φ');ax.grid(alpha=.15);ax.legend(title='Cluster size',fontsize=8)
    axes[0].set_ylabel('Mean H₂ endpoint fraction (%)')
    fig.suptitle('Fixed cluster-size bins across Φ; all clusters, only bins with >=50 contributing seeds')
    save(fig,'10_fixed_sizes_across_phi')
    table=[];diagnostics=[]
    for key,v in sorted(data.items()):
        L,p,scope,lo,hi=key
        table.append(dict(L=L,target_phi=p/100,scope=scope,size_low=lo,size_high=hi,mean_percent=v['mean'],se_percentage_points=v['se'],contributing_seeds=v['n'],cluster_count=v['clusters'],frequent=v['n']>=math.ceil(v['total']/2)))
    for L in sizes:
        for scope in scopes:
            eligible=greater=0;opposite=[]
            for p in phis:
                a=data.get((L,p,scope,2,2));b=data.get((L,p,scope,5,8))
                if a and b and min(a['n'],b['n'])>=50:
                    eligible+=1
                    if a['mean']>b['mean']:greater+=1
                    else:opposite.append(p/100)
            diagnostics.append(dict(L=L,scope=scope,eligible_phi_count=eligible,mean_size2_greater_than_size5to8=greater,opposite_phi=opposite,note='Descriptive comparison of means only; no hypothesis test or causality.'))
    with (out/'all_phi_table.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    metadata=dict(source=str(source),size_relation_sha256=hashlib.sha256((source/'size_relation.csv').read_bytes()).hexdigest(),rows=len(data),
        representative_phis=[p/100 for p in representatives],rare_threshold=50,diagnostics=diagnostics,
        limits='Means at different bins may condition on different sets of seeds; no individual clusters in this export; ratios alone do not identify causal growth effects.')
    (out/'all_phi_manifest.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    archive=out/'upload_h2_all_phi.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for p in sorted(out.iterdir()):
            if p!=archive:z.write(p,p.name)
    print(json.dumps(diagnostics,ensure_ascii=False,indent=2))
    print('OUTPUT:',out,'\nUPLOAD THIS FILE:\n'+str(archive))

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError,ImportError) as exc:sys.exit('ERROR: '+str(exc))
