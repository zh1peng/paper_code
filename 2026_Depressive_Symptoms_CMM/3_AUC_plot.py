"""Render Figure 3, including optimism-corrected AUC and delta AUC."""
from lib.plotting import *
d=read('discrimination_bootstrap').set_index('cohort').loc[ORDER].reset_index()
flow=read('cross_sectional_flow').set_index('cohort')
assert d.replicates.nunique()==1
replicates=int(d.replicates.iloc[0])
assert all(int(r.n)==int(flow.loc[r.cohort,'complete']) for _,r in d.iterrows())
fig=plt.figure(figsize=(7.1,3.65))
fig.text(.02,.97,'A  Cross-sectional discrimination',fontsize=10.8,weight='bold',va='top')
fig.text(.535,.97,'B  Incremental value',fontsize=10.8,weight='bold',va='top')
ax=fig.add_axes([.09,.16,.355,.66]);ax.set_xlim(.697,.824);ax.set_ylim(.4,4.7)
for y,(_,r) in zip([4,3,2,1],d.iterrows()):
    ax.axhline(y,color=PALE,lw=.65)
    ax.plot([r.AUC_core,r.AUC_added],[y,y],color=MIDGRAY,lw=1.5)
    ax.plot(r.AUC_core,y,'o',mfc='white',mec=GRAY,ms=5,mew=1)
    ax.plot(r.AUC_added,y,'o',color=NAVY,ms=4.5)
    ax.text(r.AUC_core-.003,y,f'{r.AUC_core:.3f}',ha='right',va='center',fontsize=7,color=GRAY)
    ax.text(r.AUC_added+.003,y,f'{r.AUC_added:.3f}',va='center',fontsize=7,color=NAVY)
    ax.plot(r.corrected_AUC_added,y-.27,'D',ms=3.8,color=MIDBLUE)
    ax.text(r.corrected_AUC_added+.003,y-.27,f'{r.corrected_AUC_added:.3f}',va='center',fontsize=7,color=BLUE)
ax.set_yticks([4,3,2,1],ORDER[:-1]+['KLoSA']);ax.tick_params(axis='y',length=0)
ax.spines['left'].set_visible(False);tick(ax,[.70,.72,.74,.76,.78,.80,.82],2);ax.set_xlabel('Area under the ROC curve')
fig.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',mec=GRAY,label='Base model'),
    Line2D([],[],marker='o',ls='',color=NAVY,label='+ depressive symptoms')],
    loc='upper left',bbox_to_anchor=(.06,.90),ncol=2,fontsize=7,frameon=False,handletextpad=.3,columnspacing=1,markerscale=.7)
fig.legend(handles=[Line2D([],[],marker='D',ls='',color=MIDBLUE,label='Optimism-corrected augmented model')],
    loc='upper left',bbox_to_anchor=(.06,.84),fontsize=7,frameon=False,handletextpad=.3,markerscale=.7)
ax=fig.add_axes([.555,.16,.425,.66]);ax.set_xlim(-.003,.101);ax.set_ylim(.4,4.7)
table_text=[]
for y,(_,r) in zip([4,3,2,1],d.iterrows()):
    ax.axhline(y,color=PALE,lw=.65)
    ax.hlines(y,0,r.delta,color=MIDBLUE,lw=1.4)
    ax.errorbar(r.delta,y,xerr=[[r.delta-r.delta_lcl],[r.delta_ucl-r.delta]],fmt='o',color=NAVY,ms=4.5,capsize=2,lw=.85)
    table_text.append(ax.text(.035,y,f'{r.delta:.3f} ({r.delta_lcl:.3f}-{r.delta_ucl:.3f})',va='center',fontsize=6.7,color=TEXT))
    table_text.append(ax.text(.079,y,f'{int(r.n):,} / {int(flow.loc[r.cohort,"events"]):,}',va='center',fontsize=6.7,color=TEXT))
ax.axvline(0,color=MIDGRAY,lw=.8);ax.set_yticks([4,3,2,1],ORDER[:-1]+['KLoSA']);ax.tick_params(axis='y',length=0)
ax.spines['left'].set_visible(False);tick(ax,[0,.01,.02,.03],2);ax.set_xlabel('Change in AUC')
table_text.append(ax.text(.035,4.52,'Change (95% CI)',weight='bold',fontsize=7.1))
table_text.append(ax.text(.079,4.52,'N / events',weight='bold',fontsize=7.1))
fig.text(.09,.035,f'Paired refitted bootstrap, {replicates:,} resamples per cohort. The AUC axis in panel A is truncated.',fontsize=6.8,color=TEXT)
fig.canvas.draw()
for i,t in enumerate(table_text):
    for u in table_text[i+1:]:assert not t.get_window_extent().overlaps(u.get_window_extent()),(t.get_text(),u.get_text())
save(fig,'Figure3_discrimination')
