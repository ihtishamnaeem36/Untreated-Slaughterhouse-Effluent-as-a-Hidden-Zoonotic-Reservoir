"""
FIGURE 4 — Clinical Annotation Barplots (Genus Count + Abundance)
Input : merged_genus_relative_abundance.csv
Output: Figure4_Annotation_Barplot_600dpi.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd, numpy as np

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CONFIRMED={'Arcobacter','Aliarcobacter','Salmonella','Brucella','Campylobacter','Listeria','Yersinia','Leptospira','Helicobacter','Treponema'}
CRITICAL={'Acinetobacter','Pseudomonas'}
OPPORTUNISTIC={'Aeromonas','Escherichia','Klebsiella','Clostridium','Enterococcus','Streptococcus','Vibrio','Fusobacterium','Shewanella','Corynebacterium','Staphylococcus','Stenotrophomonas','Bacteroides'}
COMMENSAL=['Bifidobacterium','Desulfovibrio','Comamonas']

def clin4(name):
    if name in CONFIRMED: return 'Zoonotic Pathogen'
    if name in CRITICAL:  return 'WHO Critical MDR'
    if name in OPPORTUNISTIC: return 'Opportunistic Zoonotic'
    return 'Commensal / Environmental'

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
df['mean']=df[SAMPLES].mean(axis=1)
df_s=df.sort_values('mean',ascending=False).reset_index(drop=True)

all_gen=list(CONFIRMED|CRITICAL|OPPORTUNISTIC)+COMMENSAL
cats=['Zoonotic Pathogen','WHO Critical MDR','Opportunistic Zoonotic','Commensal / Environmental']
counts={'Zoonotic Pathogen':10,'WHO Critical MDR':2,'Opportunistic Zoonotic':13,'Commensal / Environmental':3}
colors=['#C00000','#7030A0','#ED7D31','#70AD47']
abund={}
for cat in cats:
    gen=[g for g in all_gen if clin4(g)==cat]
    abund[cat]=sum(df_s[df_s['name']==g]['mean'].values[0]*100 for g in gen if len(df_s[df_s['name']==g])>0)
count_vals=[counts[c] for c in cats]; abund_vals=[abund[c] for c in cats]
short_lab=['Zoonotic\nPathogen','WHO Critical\nMDR','Opportunistic\nZoonotic','Commensal /\nEnvironmental']

fig,axes=plt.subplots(1,2,figsize=(18,9)); fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.07,right=0.78,wspace=0.35,top=0.88,bottom=0.12)
for ax,vals,ylabel,title_txt in zip(axes,[count_vals,abund_vals],
    ['Number of Genera','Mean Relative Abundance (%)'],
    ['A — Genus Count by Clinical Category\n(n = 28 total annotated genera)',
     'B — Mean Relative Abundance by Clinical Category\n(% of total classified reads)']):
    ax.set_facecolor('#FAFAFA')
    bars=ax.bar(range(4),vals,color=colors,width=0.50,edgecolor='white',linewidth=1.2)
    for bar,v in zip(bars,vals):
        ax.text(bar.get_x()+bar.get_width()/2,v+max(vals)*0.025,
                str(int(v)) if ylabel=='Number of Genera' else f'{v:.2f}%',
                ha='center',va='bottom',fontsize=13,fontweight='bold',color='#212529')
    for bar,v,tot in zip(bars,vals,[28,sum(abund_vals)]):
        if v>max(vals)*0.12:
            ax.text(bar.get_x()+bar.get_width()/2,v/2,f'{v/tot*100:.1f}%\nof total',
                    ha='center',va='center',fontsize=9.5,fontweight='bold',color='white')
    ax.set_xticks(range(4)); ax.set_xticklabels(short_lab,fontsize=11,fontweight='bold',linespacing=1.4)
    for tick,col in zip(ax.get_xticklabels(),colors): tick.set_color(col)
    ax.set_ylabel(ylabel,fontsize=12,labelpad=8); ax.set_ylim(0,max(vals)*1.32)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(axis='y',alpha=0.25,zorder=1); ax.tick_params(axis='x',length=0,pad=10)
    ax.set_title(title_txt,fontsize=11,fontweight='bold',pad=14)
leg4=[mpatches.Patch(color=c,label=l) for c,l in zip(colors,
    ['Zoonotic Pathogen (n=10)','WHO Critical MDR (n=2)',
     'Opportunistic Zoonotic (n=13)','Commensal / Environmental (n=3)'])]
fig.legend(handles=leg4,loc='center right',bbox_to_anchor=(0.99,0.52),fontsize=11,
           frameon=True,framealpha=0.96,edgecolor='#ccc',borderpad=1.0,labelspacing=0.8,
           title='Clinical Classification\n(WHO / WOAH / FAO / CDC)',title_fontsize=11)
plt.savefig('Figure4_Annotation_Barplot_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 4 saved.')
