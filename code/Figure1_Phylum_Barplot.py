"""
============================================================
FIGURE 1 — Phylum-Level Relative Abundance
Slaughterhouse Wastewater — Mardan, Peshawar, Swat, KPK
============================================================
Input : merged_phylum_relative_abundance.csv
Output: Figure1_Phylum_Barplot_600dpi.png
Run   : python Figure1_Phylum_Barplot.py
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pandas as pd, numpy as np

DPI = 600
plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES  = ['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITY_COL = {'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}

ph = pd.read_csv('merged_phylum_relative_abundance.csv')
ph = ph[ph['name'].str.lower()!='unclassified'].copy()
ph['mean'] = ph[SAMPLES].mean(axis=1)
ph_s = ph.sort_values('mean',ascending=False).reset_index(drop=True)
ph_plot = ph_s.head(10)[['name']+SAMPLES].copy()
others  = ph_s.iloc[10:][SAMPLES].sum()
ph_plot = pd.concat([ph_plot,
    pd.DataFrame([['Others']+list(others)],columns=['name']+SAMPLES)],ignore_index=True)

# Campylobacterota = red (zoonotic phylum highlight)
PHYLUM_COLORS = ['#1F4E79','#2E75B6','#4BACC6','#548235','#ED7D31',
                 '#C00000','#A9D18E','#F4B183','#9DC3E6','#B4C7E7','#D9D9D9']

fig = plt.figure(figsize=(20,9)); fig.patch.set_facecolor('white')
ax     = fig.add_axes([0.05,0.11,0.50,0.80])
ax_leg = fig.add_axes([0.57,0.08,0.40,0.88])
ax.set_facecolor('#FAFAFA'); ax_leg.axis('off')

x = np.arange(6); bottoms = np.zeros(6)
for i,row in ph_plot.iterrows():
    vals = row[SAMPLES].values.astype(float)
    ax.bar(x,vals,bottom=bottoms,color=PHYLUM_COLORS[i],width=0.62,edgecolor='white',linewidth=0.6,zorder=3)
    for j,(v,b) in enumerate(zip(vals,bottoms)):
        if v>=0.03:
            dark={'#1F4E79','#548235','#C00000','#2E75B6'}
            tc='white' if PHYLUM_COLORS[i] in dark else '#212529'
            ax.text(x[j],b+v/2,f'{v*100:.1f}%',ha='center',va='center',fontsize=7.5,fontweight='bold',color=tc,zorder=4)
    bottoms+=vals
for j,tot in enumerate(ph_plot[SAMPLES].sum()):
    ax.text(x[j],tot+0.008,f'{tot*100:.1f}%',ha='center',va='bottom',fontsize=10,fontweight='bold',color='#212529')
bracket_y=1.055
for city,idx in [('Mardan',[0,1]),('Peshawar',[2,3]),('Swat',[4,5])]:
    xs=idx[0]-0.38; xe=idx[1]+0.38
    ax.annotate('',xy=(xe,bracket_y),xytext=(xs,bracket_y),
                arrowprops=dict(arrowstyle='-',color=CITY_COL[city],lw=2.2))
    ax.text((xs+xe)/2,bracket_y+0.01,city,ha='center',va='bottom',fontsize=11,fontweight='bold',color=CITY_COL[city])
ax.set_xticks(x); ax.set_xticklabels(SAMPLES,fontsize=11,fontweight='bold')
ax.set_ylabel('Relative Abundance',fontsize=12,labelpad=8)
ax.set_xlabel('Sample',fontsize=12,labelpad=8)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1,decimals=0))
ax.set_xlim(-0.55,5.55); ax.set_ylim(0,1.12)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(axis='y',alpha=0.2,zorder=1)
legend_patches=[]
for i,row in ph_plot.iterrows():
    name=row['name']
    mv=ph_s[ph_s['name']==name]['mean'].values[0]*100 if name in list(ph_s['name']) else 0
    note='  ★ Zoonotic phylum' if name=='Campylobacterota' else ''
    label=f'{name}  ({mv:.1f}%)' if name!='Others' else 'Others (remaining phyla)'
    legend_patches.append(mpatches.Patch(color=PHYLUM_COLORS[i],label=f'{label}{note}'))
ax_leg.legend(handles=legend_patches,loc='upper left',fontsize=10,frameon=True,framealpha=0.95,
              edgecolor='#cccccc',borderpad=1.0,labelspacing=0.75,
              title='Phylum  (mean % shown)',title_fontsize=10.5)
plt.savefig('Figure1_Phylum_Barplot_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 1 saved.')
