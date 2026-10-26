"""Optional scientific plots (Matplotlib is imported only on request)."""
from pathlib import Path

def make_plots(summary,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    styles={32:('#3068a0','o'),64:('#db8635','s'),128:('#258471','^')}
    configs=[('H2_per_volume',r'$H_2/L^3$',1,'h2_density'),
             ('H2_endpoint_fraction_bonded','Bonded endpoints in exactly-two-bond pairs (%)',100,'h2_endpoint_bonded'),
             ('H2_endpoint_fraction_all','All endpoints in exactly-two-bond pairs (%)',100,'h2_endpoint_all'),
             ('S1',r'$S_1=(1+4l)N_1/L^3$',1,'largest_lattice_fraction'),
             ('cluster_number_mean',r'Number-average cluster size $N/K$',1,'cluster_number_mean'),
             ('cluster_weighted_excl_one',r'Weighted cluster size (one largest excluded)',1,'cluster_weighted')]
    plt.rcParams.update({'font.size':13,'svg.fonttype':'none'})
    for key,label,factor,name in configs:
        if not any(key+'_mean' in r for r in summary):continue
        fig,ax=plt.subplots(figsize=(8.6,5.1),layout='constrained')
        for L in sorted({r['L'] for r in summary}):
            g=[r for r in summary if r['L']==L and key+'_mean' in r]
            if not g:continue
            color,marker=styles.get(L,('black','o'))
            ax.errorbar([r['actual_phi_mean'] for r in g], [factor*r[key+'_mean'] for r in g],
                        yerr=[factor*(r.get(key+'_sd') or 0) for r in g],color=color,marker=marker,capsize=3,label=f'L = {L}')
        ax.set(xlabel=r'Actual occupied fraction $\phi$',ylabel=label,ylim=(0,None))
        ax.grid(alpha=.22);ax.legend();
        for ext in ['png','svg']:fig.savefig(Path(out)/(name+'.'+ext),dpi=190)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    available_phis=sorted({r['target_phi'] for r in summary})
    selected_phis=[v for v in available_phis if any(abs(v-x)<1e-9 for x in [.3,.4,.5,.7])] or available_phis[:3]
    sizes=sorted({r['L'] for r in summary})
    for phi in selected_phis:
        g=[r for r in summary if abs(r['target_phi']-phi)<1e-9]
        for ax,key,label in zip(axes,['H2','H2_per_volume'],[r'$H_2$',r'$H_2/L^3$']):
            ax.errorbar([r['L']**3 for r in g],[r[key+'_mean'] for r in g],yerr=[r[key+'_sd'] or 0 for r in g],marker='o',capsize=3,label=f'phi = {phi}')
            ax.set(xscale='log',xlabel=r'Volume $L^3$',ylabel=label,ylim=(0,None))
            ax.set_xticks([L**3 for L in sizes],[f'{L}³' for L in sizes]);ax.grid(alpha=.2)
    axes[0].legend(fontsize=10)
    for ext in ['png','svg']:fig.savefig(Path(out)/('h2_size.'+ext),dpi=190)
    plt.close(fig)
