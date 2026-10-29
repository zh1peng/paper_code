"""Render Figure 2 from newly computed main-analysis tables."""
from lib.plotting import *
est=read('cross_sectional_estimates');pool=read('cross_sectional_meta')
d=est[est.analysis=='Primary binary'].set_index('cohort').loc[ORDER].reset_index()
pp=pool[pool.analysis=='Primary binary'].iloc[0].to_dict();pp['cohort']='Pooled'
d=pd.concat([d,pd.DataFrame([pp])],ignore_index=True)
fig=plt.figure(figsize=(7.1,4.6),label='Figure2')
title(fig,.02,.973,'A  Prevalent CMM');title(fig,.535,.973,'B  Individual conditions')
ax=fig.add_axes([.09,.14,.40,.72]);forest(ax,d,(.85,10),[1,1.5,2,3,4],'Adjusted odds ratio',4.5,True)
cp=read('component_estimates');ax=fig.add_axes([.615,.14,.365,.72])
ax.set_xlim(.9,3.65);ax.set_ylim(.2,15.5);ax.axvline(1,color=MIDGRAY,ls=(0,(4,4)),lw=.8)
ticks=[];labels=[]
for k,oc in enumerate(['diabetes','heart_disease','stroke']):
    top=14.3-k*5
    ax.text(-.19,top+.65,oc.replace('_',' ').capitalize(),transform=ax.get_yaxis_transform(),fontsize=7.7,weight='bold')
    for j,co in enumerate(ORDER):
        rr=cp[(cp.cohort==co)&(cp.analysis==oc)].iloc[0];y=top-j
        ax.axhline(y,color=PALE,lw=.5)
        ax.errorbar(rr.ratio,y,xerr=[[rr.ratio-rr.lcl],[rr.ucl-rr.ratio]],fmt='o',color=GRAY,ms=3.3,lw=.9,capsize=1.7)
        ax.text(2.02,y,ci(rr),va='center',fontsize=7.0,color=TEXT,bbox={'fc':'white','ec':'none','pad':.3})
        ticks.append(y);labels.append(name(co))
ax.set_yticks(ticks,labels,fontsize=7);ax.tick_params(axis='y',length=0,pad=2);ax.spines['left'].set_visible(False)
tick(ax,[1,1.4,1.8]);ax.set_xlabel('Odds ratio per 1-SD higher score',fontsize=8)
ax.text(2.02,15.25,'OR (95% CI)',fontsize=7.5,weight='bold')
save(fig,'Figure2_cross_sectional')
