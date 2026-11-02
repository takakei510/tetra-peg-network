"""Report figures: Python3 + matplotlib. Input is the compact size_relation.csv.
python3 plot_report.py --input size_relation.csv --output figures
"""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

def japanese_font(explicit=None):
    if explicit:
        font_manager.fontManager.addfont(explicit)
        name=font_manager.FontProperties(fname=explicit).get_name()
    else:
        names={f.name for f in font_manager.fontManager.ttflist}
        name=next((n for n in ['Noto Sans CJK JP','Noto Sans JP','IPAexGothic','Yu Gothic','Meiryo'] if n in names),None)
        if name is None:
            candidates=[Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'),Path('/mnt/c/Windows/Fonts/meiryo.ttc')]
            f=next((f for f in candidates if f.is_file()),None)
            if f:
                font_manager.fontManager.addfont(str(f));name=font_manager.FontProperties(fname=str(f)).get_name()
        if name is None:raise RuntimeError('日本語フォントが必要です。Ubuntu/WSL: sudo apt install fonts-noto-cjk または --font フォントファイル')
    plt.rcParams.update({'font.family':name,'font.size':11,'axes.unicode_minus':False,'svg.fonttype':'path','axes.spines.top':False,'axes.spines.right':False})

def plot(args):
    japanese_font(args.font)
    rows=list(csv.DictReader(Path(args.input).open(encoding='utf-8-sig')))
    phis=args.phi;L=args.L
    data={(round(float(r['target_phi']),2),r['scope'],int(r['size_low'])):r for r in rows if int(r['L'])==L}
    allbins=sorted({(int(r['size_low']),int(r['size_high'])) for r in rows})
    for p in phis:
        if not any(k[0]==round(p,2) for k in data):raise ValueError(f'L={L}, phi={p}: データなし')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    colors=['#2864a0','#bd601c','#25816c']
    for kind in ['count_log','count_linear','h2_fraction']:
        count=kind.startswith('count');scope='all' if count else 'exclude_one_largest'
        used=[i for i,(lo,hi) in enumerate(allbins) if any(data.get((round(p,2),scope,lo)) for p in phis)]
        bins=allbins[:max(used)+1]
        if not count:bins=[b for b in bins if b[0]>1]
        labels=[str(lo) if lo==hi else f'{lo:,}〜{hi:,}' for lo,hi in bins]
        n=len(phis)
        fig=plt.figure(figsize=(12,3.0*n+1.4))
        gs=fig.add_gridspec(n*2,1,height_ratios=[1,.13]*n,hspace=.25)
        axes=[];maxcount=max(int(r['clusters'])/int(r['total_samples']) for k,r in data.items() if k[0] in phis and k[1]=='all')
        for j,p in enumerate(phis):
            ax=fig.add_subplot(gs[2*j,0],sharex=axes[0] if axes else None);axes.append(ax)
            axn=fig.add_subplot(gs[2*j+1,0],sharex=ax)
            color=colors[j%len(colors)]
            for x,(lo,hi) in enumerate(bins):
                r=data.get((round(p,2),scope,lo))
                if count:
                    y=int(r['clusters'])/int(r['total_samples']) if r else 0
                    if y>0 or kind=='count_linear':ax.scatter(x,y,color=color,s=38,zorder=3)
                    axn.text(x,.5,'0' if y==0 else '',ha='center',va='center',color='#7d8790',fontsize=10)
                else:
                    if r:
                        y=100*float(r['fraction_mean']);ns=int(r['fraction_n'])
                        if ns>=50:ax.errorbar(x,y,yerr=100*float(r['fraction_se']),fmt='o',color=color,ms=6,capsize=3)
                        else:ax.scatter(x,y,facecolors='none',edgecolors='#999999',s=35)
                        s=str(ns)
                    else:s='—'
                    axn.text(x,.5,s,ha='center',va='center',color='#556575',fontsize=10)
            if kind=='count_log':ax.set_yscale('log');ax.set_ylim(.005,10**(__import__('math').ceil(__import__('math').log10(maxcount))))
            elif kind=='count_linear':ax.set_ylim(-.015*maxcount,1.1*maxcount)
            else:ax.set_ylim(0,11)
            ax.set_ylabel('平均クラスタ数\n［個／試料］' if count else '$H_2$末端割合［％］')
            ax.text(.99,.85,f'φ = {p:.2f}',transform=ax.transAxes,ha='right',fontsize=15,color=color)
            ax.grid(axis='y',alpha=.2);ax.set_xlim(-.65,len(bins)-.35)
            ax.tick_params(axis='x',labelbottom=False,bottom=False)
            axn.set_ylim(0,1);axn.set_yticks([])
            for sp in axn.spines.values():sp.set_visible(False)
            axn.set_ylabel('0個の区間' if count else '寄与試料数',rotation=0,ha='right',va='center',labelpad=12,fontsize=10)
            if j==n-1:
                axn.set_xticks(range(len(bins)),labels,rotation=40,ha='right',fontsize=10)
                axn.set_xlabel('クラスタサイズ区間［分子］',labelpad=9)
            else:axn.tick_params(axis='x',bottom=False,labelbottom=False)
        title='クラスタサイズ区間ごとの平均個数' if count else '最大クラスタを除いた$H_2$末端割合のサイズ依存'
        scope_text='全クラスタ' if count else '各試料の最大クラスタを1個除外'
        fig.suptitle(f'{title}\nL = {L}　｜　{scope_text}　｜　各φ：100試料',fontsize=15,y=.98)
        if count:
            note=('縦軸：対数目盛。' if kind=='count_log' else '縦軸：通常目盛。')+'区間内のクラスタ数の合計 ÷ 100試料（該当なしは0として集計）。\n区間幅では割っていない。個数のSEは未算出。下段の0は、その区間のクラスタが100試料で0個。'
        else:
            note='色付き点：平均 ± SE（寄与試料数50以上）。白抜き点：50未満。—：該当クラスタなし。\n試料内でクラスタを等重み平均し、次に寄与試料を等重み平均。$H_2$=0も含む。1分子は常に0のため省略。'
        note+='\n横軸は区間を等間隔に配置したカテゴリ軸。代表3条件の記述であり、全φの傾向・因果はこの図だけでは判断しない。'
        fig.text(.14,.025,note,fontsize=9,color='#45576a',va='bottom')
        fig.subplots_adjust(left=.14,right=.98,top=.88,bottom=.26)
        for ext in ['png','svg']:fig.savefig(out/f'{kind}.{ext}',dpi=200,bbox_inches='tight',pad_inches=.2)
        plt.close(fig)
    print(f'出力先: {out.resolve()}')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',default='size_relation.csv');p.add_argument('--output',default='figures')
    p.add_argument('--L',type=int,default=128);p.add_argument('--phi',nargs='+',type=float,default=[.41,.43,.60]);p.add_argument('--font')
    plot(p.parse_args())
