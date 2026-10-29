"""Render Figure 6 joint-category prevalence ratios and cell counts."""
from lib.plotting import *
d=read('biomarker_groups');d=d[d.analysis=='Primary']
cmap=LinearSegmentedColormap.from_list('submitted',[PALE,'#C4D2E5',MIDBLUE,NAVY]);norm=Normalize(1,3.1,clip=True)
fig=plt.figure(figsize=(7.1,2.70))
for j,co in enumerate(ORDER[:3]):
    ax=fig.add_axes([.105+j*.33,.245,.218,.61]);z=d[d.cohort==co]
    for _,r in z.iterrows():
        x=int(r.high_inflammation);y=1-int(r.depressive_symptoms)
        ax.add_patch(Rectangle((x-.5,y-.5),1,1,facecolor=cmap(norm(r.PR)),edgecolor='white',lw=2))
        ink='white' if norm(r.PR)>.82 else 'black'
        label='Reference' if r.depressive_symptoms==r.high_inflammation==0 else f'{r.PR:.2f}\n({r.lcl:.2f}-{r.ucl:.2f})'
        ax.text(x,y+.08,label,ha='center',va='center',fontsize=7.1,color=ink,linespacing=1.1)
        ax.text(x,y-.29,f'{int(r.cases):,} / {int(r.n):,}',ha='center',va='center',fontsize=6.5,color=ink)
    ax.set_xlim(-.5,1.5);ax.set_ylim(-.5,1.5);ax.set_xticks([0,1],['Low\ninflammation','High\ninflammation'],fontsize=7.3)
    ax.set_yticks([0,1],['Elevated\nsymptoms','Below\nthreshold'],fontsize=7.3);ax.tick_params(length=0,pad=5)
    for s in ax.spines.values():s.set_visible(False)
    ax.set_title(co,pad=12)
fig.text(.105,.045,'Cells show prevalence ratio (95% CI) and CMM cases / participants.',fontsize=7,color=TEXT)
save(fig,'Figure6_joint_categories')
