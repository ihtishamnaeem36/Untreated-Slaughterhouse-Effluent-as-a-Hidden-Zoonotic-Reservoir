"""
FIGURE 5 — Zoonotic Genera Abundance by City
Input : merged_genus_relative_abundance.csv
Output: Figure5_ZoonoticCityBar_600dpi.png
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

zoo25=[(g,df_s[df_s['name']==g].iloc[0]['Mardan']*100,
          df_s[df_s['name']==g].iloc[0]['Peshawar']*100,
          df_s[df_s['name']==g].iloc[0]['Swat']*100,clin(g))
       for g in ZOONOTIC_ALL if len(df_s[df_s['name']==g])>0]
zoo25.sort(key=lambda x:(x[1]+x[2]+x[3])/3,reverse=True)
genera5=[d[0] for d in zoo25]; mar_v=[d[1] for d in zoo25]
pes_v=[d[2] for d in zoo25]; swa_v=[d[3] for d in zoo25]; st5=[d[4] for d in zoo25]

fig=plt.figure(figsize=(18,14)); fig.patch.set_facecolor('white')
ax    =fig.add_axes([0.18,0.04,0.55,0.94])
ax_leg=fig.add_axes([0.75,0.60,0.23,0.36])
ax.set_facecolor('#FAFAFA'); ax_leg.axis('off')

y=np.arange(len(genera5)); h=0.24
ax.barh(y+h,mar_v,height=h,color='#1F4E79',alpha=0.92)
ax.barh(y,  pes_v,height=h,color='#C55A11',alpha=0.92)
ax.barh(y-h,swa_v,height=h,color='#375623',alpha=0.92)
ax.set_yticks(y)
ylabels5=[f'{CLIN_SYM.get(s,"  ")}{g}' for g,s in zip(genera5,st5)]
ax.set_yticklabels(ylabels5,fontsize=10.5,style='italic',fontweight='bold')
for lbl,s in zip(ax.get_yticklabels(),st5): lbl.set_color(CLIN_COL[s])
ax.set_xlabel('Mean Relative Abundance (%)',fontsize=12,labelpad=8)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(axis='x',alpha=0.25); ax.tick_params(axis='y',length=0,pad=10)
ax.set_xlim(0,max(max(mar_v),max(pes_v),max(swa_v))*1.15)

leg5=[mpatches.Patch(color='#1F4E79',label='Mardan'),
      mpatches.Patch(color='#C55A11',label='Peshawar'),
      mpatches.Patch(color='#375623',label='Swat'),
      mpatches.Patch(color='none',label=''),
      mpatches.Patch(color='#C00000',label='▶  Zoonotic Pathogen'),
      mpatches.Patch(color='#7030A0',label='★  WHO Critical MDR'),
      mpatches.Patch(color='#ED7D31',label='◆  Opportunistic Zoonotic')]
ax_leg.legend(handles=leg5,loc='upper left',fontsize=10.5,frameon=True,framealpha=0.95,
              edgecolor='#ccc',borderpad=0.9,labelspacing=0.6,
              title='City / Clinical Category',title_fontsize=11)
plt.savefig('Figure5_ZoonoticCityBar_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 5 saved.')
