"""Render Figure 4 forest plot and exploratory spline panels."""
from lib.plotting import *
le=read('longitudinal_estimates');lp=read('longitudinal_meta')
d=le[le.analysis=='Discrete-time primary'].set_index('cohort').loc[LONG].reset_index()
pp=lp.iloc[0].to_dict();pp['cohort']='Pooled';d=pd.concat([d,pd.DataFrame([pp])],ignore_index=True)
cur=read('longitudinal_splines')
def spline(ax,z,ylabel,ylim):
    ax.fill_between(z.score,z.lcl,z.ucl,color=MIDBLUE,alpha=.22,lw=0)
    ax.plot(z.score,z.ratio,color=NAVY,lw=1.3);ax.axhline(1,color=MIDGRAY,ls=(0,(4,4)),lw=.8)
    ax.set_yscale('log');ax.set_ylim(*ylim);ax.set_xlim(z.score.min(),z.score.max())
    ax.yaxis.set_major_locator(FixedLocator([.5,1,2,3,4,6,10]));ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'))
    ax.yaxis.set_minor_locator(NullLocator())
    # FixedLocator can expand limits to its ticks, so restore the intended common limits.
    ax.set_ylim(*ylim);ax.set_xlabel('Baseline symptom score',fontsize=8)
    if ylabel:ax.set_ylabel(ylabel)
    else:ax.tick_params(axis='y',left=False,labelleft=False);ax.spines['left'].set_visible(False)
fig=plt.figure(figsize=(7.1,5.25),label='Figure4');title(fig,.02,.978,'A  Adjusted prospective association')
ax=fig.add_axes([.095,.605,.88,.305]);forest(ax,d,(.94,1.72),[1,1.2,1.4],'Grouped hazard ratio per 1-SD higher symptom score',1.52)
ax.text(.95,.47,f'Random-effects I² = {lp.iloc[0].I2:.1f}%',fontsize=7.5,color=TEXT,va='bottom')
title(fig,.02,.502,'B  Dose-response associations')
ylim=(min(.4,cur.lcl.min()*.9),max(3.05,cur.ucl.max()*1.05))
for j,co in enumerate(LONG):
    ax=fig.add_axes([.095+j*.303,.12,.267,.30]);z=cur[cur.cohort==co]
    spline(ax,z,'Grouped hazard ratio' if j==0 else None,ylim)
    ax.set_title(name(co)+', '+('CES-D-8' if co=='HRS' else 'CES-D-10'),pad=8)
    tick(ax,[0,2,4,6,8] if co=='HRS' else [0,5,10,15,20,25,30],0)
fig.text(.095,.027,'Shaded areas show 95% confidence intervals. Spline shapes are exploratory.',fontsize=7,color=TEXT)
save(fig,'Figure4_prospective')
