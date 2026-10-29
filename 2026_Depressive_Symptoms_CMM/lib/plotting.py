from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from matplotlib.lines import Line2D

OUTPUT = Path(sys.argv[1] if len(sys.argv)>1 else 'output')
TAB = OUTPUT / 'tables'
OUT = OUTPUT / 'figures'
OUT.mkdir(parents=True,exist_ok=True)
NAVY='#17365D'; BLUE='#2F5597'; MIDBLUE='#6F91BC'
GRAY='#666666'; MIDGRAY='#A6A6A6'; LIGHT='#D9D9D9'; PALE='#F2F2F2'; TEXT='#404040'
ORDER=['NHANES','CHARLS','HRS','KLOSA']; LONG=ORDER[1:]
NAMES={'KLOSA':'KLoSA'}
def read(name):
    return pd.read_csv(TAB/(name+'.csv'))

plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8.5,
 'axes.titlesize':9,'axes.titleweight':'bold','axes.linewidth':.65,
 'xtick.labelsize':8,'ytick.labelsize':8,'xtick.major.width':.6,'ytick.major.width':.6,
 'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white',
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False})

def name(co):return NAMES.get(co,co)
def ci(r):return f'{r.ratio:.2f} ({r.lcl:.2f}-{r.ucl:.2f})'
def title(fig,x,y,text):fig.text(x,y,text,fontsize=10.8,fontweight='bold',va='top')
def tick(ax,values,decimals=1):
    ax.xaxis.set_major_locator(FixedLocator(values));ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.{decimals}f}'))
    ax.xaxis.set_minor_locator(NullLocator())
def save(fig,stem):
    fig.canvas.draw()
    # Catch text extending outside the physical canvas before export.
    renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height()
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():continue
        box=t.get_window_extent(renderer)
        if box.width and box.height:
            assert box.x0>=-2 and box.y0>=-2 and box.x1<=w+2 and box.y1<=h+2,(stem,t.get_text(),box.bounds,(w,h))
    for ext in ['png','pdf','tif']:
        kw={'dpi':300 if ext=='png' else 600}
        if ext=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(OUT/(stem+'.'+ext),**kw)
    plt.close(fig)

def forest(ax,d,xlim,xticks,xlabel,ci_x,log=False):
    ys=np.arange(len(d),0,-1)
    ax.set_xlim(*xlim);ax.set_ylim(.4,len(d)+.7)
    if log:ax.set_xscale('log')
    for y,(_,r) in zip(ys,d.iterrows()):
        ax.axhline(y,color=PALE,lw=.65,zorder=0)
        ax.errorbar(r.ratio,y,xerr=[[r.ratio-r.lcl],[r.ucl-r.ratio]],fmt='D' if r.cohort=='Pooled' else 'o',
                    color=NAVY if r.cohort=='Pooled' else GRAY,ecolor=GRAY,lw=1.05,capsize=2.2,ms=4.5 if r.cohort=='Pooled' else 4)
        ax.text(ci_x,y,ci(r),va='center',fontsize=7.7,color=TEXT,bbox={'fc':'white','ec':'none','pad':.5})
    ax.axvline(1,color=MIDGRAY,ls=(0,(4,4)),lw=.8,zorder=0)
    ax.set_yticks(ys,[name(c) for c in d.cohort]);ax.tick_params(axis='y',length=0,pad=3)
    ax.spines['left'].set_visible(False);tick(ax,xticks);ax.set_xlabel(xlabel,labelpad=4)
    ax.text(ci_x,len(d)+.48,'HR (95% CI)' if 'hazard' in xlabel.lower() else 'OR (95% CI)',weight='bold',fontsize=8)
