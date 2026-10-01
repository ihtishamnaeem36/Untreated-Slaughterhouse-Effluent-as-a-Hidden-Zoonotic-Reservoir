"""
FIGURE 9 — Differential Enrichment of Zoonotic Genera by City
Input : merged_genus_relative_abundance.csv
Output: Figure9_ZoonoticEnrichment_600dpi.png
NOTE  : LEfSe not applicable at n=2 per city — descriptive enrichment ratios used
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
ZOONOTIC_ALL=list(CONFIRMED|CRITICAL|OPPORTUNISTIC)
CLIN_COL={'zoonotic':'#C00000','critical':'#7030A0','opportunistic':'#ED7D31'}
CLIN_SYM={'zoonotic':'▶ ','critical':'★ ','opportunistic':'◆ '}

def clin(name):
    if name in CONFIRMED: return 'zoonotic'
    if name in CRITICAL:  return 'critical'
    return 'opportunistic'

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
df['mean']=df[SAMPLES].mean(axis=1)
df['Mardan']  =df[['MSLW1','MSLW2']].mean(axis=1)
df['Peshawar']=df[['PSLW1','PSLW2']].mean(axis=1)
df['Swat']    =df[['SSLW1','SSLW2']].mean(axis=1)
df_s=df.sort_values('mean',ascending=False).reset_index(drop=True)

fig=plt.figure(figsize=(26,15)); fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.04,right=0.98,top=0.93,bottom=0.18,wspace=0.55)
axes9=[fig.add_subplot(1,3,i+1) for i in range(3)]
CITY_CFG=[('Mardan',['Peshawar','Swat'],'#1F4E79'),
          ('Peshawar',['Mardan','Swat'],'#C55A11'),
          ('Swat',['Mardan','Peshawar'],'#375623')]

for ax,(city,others,city_color) in zip(axes9,CITY_CFG):
    ax.set_facecolor('#FAFAFA')
    enrichments=[]
    for g in ZOONOTIC_ALL:
        row=df_s[df_s['name']==g]
        if len(row):
            r=row.iloc[0]; cv=float(r[city]); om=np.mean([float(r[o]) for o in others])
            ratio=cv/om if om>0 else 0
            enrichments.append({'genus':g,'ratio':ratio,'abund':cv*100,'status':clin(g)})
    enrichments.sort(key=lambda x:x['ratio'],reverse=True)
    top10=enrichments[:10]; top10.reverse()
    gen9=[e['genus'] for e in top10]; rat9=[e['ratio'] for e in top10]
    ab9=[e['abund'] for e in top10]; st9=[e['status'] for e in top10]
    col9=[CLIN_COL[s] for s in st9]
    use_log=city=='Peshawar'
    pv=[np.log10(max(r,0.01)) if use_log else r for r in rat9]
    mv=max(pv) if pv else 1
    ax.barh(range(10),pv,color=col9,height=0.55,edgecolor='white',linewidth=0.8,alpha=0.9)
    for i,(p,r,a) in enumerate(zip(pv,rat9,ab9)):
        rt=f'{r:.1f}×' if r<1000 else f'{r:.0f}×'
        at=f'{a:.3f}%' if a>=0.001 else f'{a:.4f}%'
        ax.text(p+mv*0.025,i+0.17,rt,va='center',ha='left',fontsize=8.5,fontweight='bold',color='#212529')
        ax.text(p+mv*0.025,i-0.20,f'({at})',va='center',ha='left',fontsize=7.5,color='#555')
    yl9=[f'{CLIN_SYM.get(s,"  ")}{g}' for g,s in zip(gen9,st9)]
    ax.set_yticks(range(10)); ax.set_yticklabels(yl9,fontsize=10.5,style='italic',fontweight='bold')
    for lbl,s in zip(ax.get_yticklabels(),st9): lbl.set_color(CLIN_COL[s])
    xlabel=('log₁₀ (Enrichment Ratio)\n* log scale — extreme values' if use_log
            else 'Enrichment Ratio\n(city mean ÷ mean of other two cities)')
    ax.set_xlabel(xlabel,fontsize=10,labelpad=8)
    zt=sum(float(df_s[df_s['name']==g][city].values[0])*100 for g in ZOONOTIC_ALL if len(df_s[df_s['name']==g])>0)
    ax.set_title(f'{city}\nTotal zoonotic burden: {zt:.2f}%',fontsize=12,fontweight='bold',color=city_color,pad=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['left'].set_visible(False)
    ax.grid(axis='x',alpha=0.25,zorder=1); ax.tick_params(axis='y',length=0,pad=10)
    ax.set_xlim(0,mv*1.90)

leg9=[mpatches.Patch(color='#C00000',label='▶  Zoonotic Pathogen (WHO/WOAH/FAO/CDC)'),
      mpatches.Patch(color='#7030A0',label='★  WHO Critical-Priority MDR Organism'),
      mpatches.Patch(color='#ED7D31',label='◆  Opportunistic Zoonotic Pathogen'),
      mpatches.Patch(color='none',   label='Bar labels show enrichment fold-change × and (city mean %)'),
      mpatches.Patch(color='none',   label='Environmental genera excluded  |  Descriptive (n=2 per city)')]
fig.legend(handles=leg9,loc='lower center',ncol=3,fontsize=10,frameon=True,framealpha=0.96,
           edgecolor='#ccc',bbox_to_anchor=(0.5,0.01),borderpad=0.9,labelspacing=0.5,
           title='Clinical Classification  |  Methods Notes',title_fontsize=10.5)
plt.savefig('Figure9_ZoonoticEnrichment_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 9 saved.')
