"""Render Figure 5 incidence curves and numbers at risk."""
from lib.plotting import *
cur=read('descriptive_cif');risk=read('cif_at_risk')
fig=plt.figure(figsize=(7.1,4.25));title(fig,.02,.977,'Cumulative incidence of CMM by baseline symptom score')
palette=['#B3B3B3','#7F9DB9','#4F81BD',NAVY]
risk_labels=[]
max_cif=cur[cur.n_risk>=.1*cur.initial_n].CIF.max();ymax=np.ceil(max_cif/.05)*.05
for j,co in enumerate(LONG):
    ax=fig.add_axes([.10+j*.315,.43,.235,.38]);ar=fig.add_axes([.10+j*.315,.13,.235,.18])
    d=cur[cur.cohort==co];rr=risk[risk.cohort==co];groups=list(d.group.unique())
    cols=palette if len(groups)==4 else [palette[0],palette[2],palette[3]]
    xmax=float(d.time.max());end=np.ceil(xmax/2)*2;times=sorted(rr.time.unique())
    times=times[::2] if len(times)>7 else times
    def score_range(g):
        a,b=g.strip('[]()').split(',');lower=int(float(a))+(1 if g.startswith('(') else 0)
        upper=int(float(b))
        return str(lower) if lower==upper else f'{lower}-{upper}'
    for k,(g,col) in enumerate(zip(groups,cols)):
        z=d[(d.group==g)&(d.n_risk>=.1*d.initial_n)]
        ax.step(z.time,z.CIF,where='post',color=col,lw=1.15,label=score_range(g))
        for t in times:
            row=rr[(rr.group==g)&(rr.time==t)];assert len(row)==1
            risk_labels.append(ar.text(t,len(groups)-k,f'{int(row.iloc[0].n):,}',ha='center',va='center',fontsize=6.1,color=TEXT))
    ax.set_ylim(0,ymax);ax.set_xlim(0,end);ax.set_title(name(co),pad=10)
    ax.set_xticks(times);ax.set_yticks(np.arange(0,ymax+.01,.1))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.1f}'))
    ax.legend(title='Symptom score',title_fontsize=6.8,fontsize=6.7,frameon=False,ncol=2,loc='upper left',columnspacing=.7,handlelength=1.4)
    if j==0:ax.set_ylabel('Cumulative incidence of detected CMM',fontsize=8)
    else:ax.tick_params(axis='y',left=False,labelleft=False);ax.spines['left'].set_visible(False)
    ar.set_xlim(0,end);ar.set_ylim(.4,len(groups)+.6)
    ar.set_yticks(range(len(groups),0,-1),[score_range(g) for g in groups],fontsize=6.4)
    ar.set_xticks(times);ar.tick_params(axis='both',length=0,pad=3,labelsize=7)
    ar.tick_params(axis='y',pad=14)
    ar.spines['left'].set_visible(False);ar.set_xlabel('Follow-up (years)',fontsize=8)
    ar.set_title('Number at risk',loc='left',fontsize=7.5,pad=3)
    risk_labels.extend(ar.get_yticklabels())
fig.text(.10,.025,'Death is a competing event. Curves stop below 10% of initial numbers at risk.',fontsize=7,color=TEXT)
fig.canvas.draw()
for i,t in enumerate(risk_labels):
    b=t.get_window_extent()
    for u in risk_labels[i+1:]:
        assert not b.overlaps(u.get_window_extent()),('Risk-table text overlaps',t.get_text(),u.get_text())
save(fig,'Figure5_cumulative_incidence')
