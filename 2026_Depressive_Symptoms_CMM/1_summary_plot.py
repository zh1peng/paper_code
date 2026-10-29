"""Render Table 1 (HTML and CSV) and Figure 1 from calculated summaries."""
from lib.plotting import *
from html import escape

table = read('table1').fillna('')
flow = read('cross_sectional_flow').set_index('cohort')
table_order = ['CHARLS','NHANES','HRS','KLOSA']
labels = {'age':'Age, years','gender':'Sex','education_group':'Education',
          'marital_status':'Marital status','smoking':'Current smoking',
          'drinking':'Alcohol consumption','bmi':'BMI, kg/m²',
          'depression_score':'Symptom score','elevated_symptoms':'Elevated symptoms',
          'component_count':'CMM components','hypertension':'Diagnosed hypertension'}
level_labels = {'Less than High School':'< high school','More than High School':'> high school',
                'High School':'High school'}
keys = [(v,level) for v in labels for level in dict.fromkeys(table.loc[table.variable==v,'level'])]
rows = []
for variable,level in keys:
    label = labels[variable] + (': '+level_labels.get(level,level) if level else '')
    cells = [label]
    for co in table_order:
        z = table[(table.cohort==co)&(table.variable==variable)&(table.level==level)]
        if len(z):
            r = z.iloc[0]
            cells += [r.without_CMM,r.with_CMM,'—' if r.SMD=='' else f'{float(r.SMD):.2f}']
        else:
            assert level=='Missing'
            cells += ['0','0','—']
    rows.append(cells)
headers = ['Characteristic'] + [name(co)+' '+x for co in table_order for x in ['Without CMM','With CMM','SMD']]
pd.DataFrame(rows,columns=headers).to_csv(TAB/'Table1_display.csv',index=False)
footnote = ('Values are unweighted mean (SD) or n (%). CMM, cardiometabolic multimorbidity; '
 'BMI, body mass index; SD, standard deviation; SMD, standardized mean difference '
 '(with CMM minus without CMM). SMDs describe the magnitude of differences rather than '
 'hypothesis-test P values. Category-specific SMDs use binary indicators. Percentages '
 'use nonmissing denominators; missing counts are shown separately. SMDs are omitted '
 'for CMM component counts because these counts define CMM status. Symptom scales '
 'and alcohol definitions differ by cohort; see Supplementary Methods.')
html = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Table 1</title>',
 '<style>body{font:12px Arial;margin:25px;color:#222}table{border-collapse:collapse;width:100%}'
 'th,td{padding:6px;border:1px solid #ddd;text-align:center;white-space:nowrap}'
 'th{background:#edf1f5}td:first-child{text-align:left}p{max-width:1300px;line-height:1.5}'
 '@page{size:A3 landscape;margin:12mm}</style>',
 '<h2>Table 1. Characteristics of the cross-sectional samples</h2><table><thead><tr>'
 '<th rowspan="2">Characteristic</th>']
html += ['<th colspan="3">'+name(co)+'</th>' for co in table_order]
html.append('</tr><tr>')
for co in table_order:
    r = flow.loc[co]
    html.append(f'<th>Without CMM<br>(n = {int(r.complete-r.events):,})</th>'
                f'<th>With CMM<br>(n = {int(r.events):,})</th><th>SMD</th>')
html.append('</tr></thead><tbody>')
for row in rows:
    html.append('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in row)+'</tr>')
html.append('</tbody></table><p>'+escape(footnote)+'</p></html>')
(TAB/'Table1.html').write_text('\n'.join(html),encoding='utf-8')
d = read('sample_flow_stages')

fig=plt.figure(figsize=(7.1,6.6))
ax=fig.add_axes([0,0,1,1]);ax.set_axis_off()
box_text=[]

def panel(letter,title,cohorts,ys,heading_y,title_y,labels,xs,width):
    data=d[d.panel==letter]
    final=data[data.stage==len(labels)]
    expected=(int(final.n.sum()),int(final.final_cases_or_events.sum()))
    fig.text(.035,title_y,letter+'  '+title,fontsize=10.3,weight='bold',va='top')
    outcome='detected CMM' if letter=='B' else 'prevalent CMM'
    fig.text(.975,title_y-.003,f'Total N = {expected[0]:,}  |  {outcome} = {expected[1]:,}',
             ha='right',va='top',fontsize=7.1,color=TEXT)
    for x,label in zip(xs,labels):
        ax.text(x+width/2,heading_y,label,ha='center',va='bottom',fontsize=7.3,
                weight='bold',color=TEXT,linespacing=1.1)
    for co,y in zip(cohorts,ys):
        rows=data[data.cohort==co].sort_values('stage')
        ax.text(.143,y,'KLoSA' if co=='KLOSA' else co,ha='right',va='center',weight='bold',fontsize=7.8)
        for j,(_,r) in enumerate(rows.iterrows()):
            last=j==len(labels)-1
            rect=Rectangle((xs[j],y-.020),width,.040,facecolor=NAVY if last else 'white',
                           edgecolor=NAVY if last else LIGHT,lw=1.1)
            ax.add_patch(rect)
            text=f'N = {r.n:,}'
            if last:text+=f'\n{"Events" if letter=="B" else "Cases"} = {int(r.final_cases_or_events):,}'
            t=ax.text(xs[j]+width/2,y,text,ha='center',va='center',fontsize=7.1,
                      color='white' if last else TEXT,weight='bold' if last else 'normal',linespacing=1.05)
            box_text.append((rect,t))
            if j:
                left=xs[j-1]+width;right=xs[j]
                ax.annotate('',xy=(right,y),xytext=(left,y),arrowprops={
                    'arrowstyle':'-|>','color':MIDGRAY,'lw':1.1,'shrinkA':.5,'shrinkB':.5})
                if r.excluded_from_previous:
                    ax.text((left+right)/2,y+.023,f'-{int(r.excluded_from_previous):,}',
                            ha='center',va='bottom',fontsize=6.7,color=GRAY)

panel('A','Cross-sectional sample',['NHANES','CHARLS','HRS','KLOSA'],[.877,.805,.733,.661],.907,.974,
      ['Age eligible','Symptoms\navailable','CMM status\navailable','Complete-case\nsample'],
      [.16,.38,.60,.82],.145)
panel('B','Prospective sample',['CHARLS','HRS','KLOSA'],[.503,.431,.359],.535,.602,
      ['Age eligible','Score, date and\ndiseases available','CMM-free at\nbaseline','Complete\ncovariates','Valid follow-up\nsample'],
      [.16,.33,.50,.67,.84],.125)
panel('C','Biomarker sample',['NHANES','CHARLS','HRS'],[.201,.135,.069],.233,.300,
      ['Age-eligible\nsource','CRP + WBC\navailable','Symptoms + CMM\navailable','Complete-case\nsample'],
      [.16,.38,.60,.82],.145)
fig.text(.035,.019,'Counts are unweighted. Numbers above arrows indicate exclusions at each step.',fontsize=7,color=TEXT)


fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for rect,text in box_text:
    a=rect.get_window_extent(renderer);b=text.get_window_extent(renderer)
    assert b.x0>=a.x0 and b.x1<=a.x1 and b.y0>=a.y0 and b.y1<=a.y1,text.get_text()
save(fig,'Figure1_sample_flow')
