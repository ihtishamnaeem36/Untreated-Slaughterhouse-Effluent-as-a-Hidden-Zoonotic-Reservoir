"""
FIGURE 6 — Alpha Diversity (Observed OTUs, Shannon, Simpson, Chao1)
Input : merged_genus_relative_abundance.csv
Output: Figure6_AlphaDiversity_600dpi.png
NOTE  : Chao1 uses pseudo-counts (abundances x1e6) — documented methodological deviation
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd, numpy as np

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
BAR_COLORS=['#1F4E79','#1F4E79','#C55A11','#C55A11','#375623','#375623']

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()

# Calculate alpha diversity from raw data
alpha={}
for s in SAMPLES:
    p=df[s].values; p_pos=p[p>0]
    shannon=-np.sum(p_pos*np.log(p_pos))
    simpson=1-np.sum(p_pos**2)
    observed=int((df[s]>0).sum())
    # Chao1 via pseudo-counts (x1e6 — documented deviation)
    counts=(p*1e6).astype(int)
    singletons=np.sum(counts==1); doubletons=max(np.sum(counts==2),1)
    chao1=observed+(singletons**2)/(2*doubletons)
    alpha[s]={'Shannon':shannon,'Simpson':simpson,'Observed':observed,'Chao1':chao1}

metrics=['Observed','Shannon','Simpson','Chao1']
ylabels6=['Observed OTUs',"Shannon H'",'Simpson (1-D)','Chao1 Richness']

fig,axes=plt.subplots(1,4,figsize=(20,7)); fig.patch.set_facecolor('white')
for ax,met,ylab in zip(axes,metrics,ylabels6):
    ax.set_facecolor('#FAFAFA')
    vals=[alpha[s][met] for s in SAMPLES]
    ax.bar(range(6),vals,color=BAR_COLORS,width=0.62,edgecolor='white',linewidth=0.6)
    for j,v in enumerate(vals):
        label=f'{v:.2f}' if met!='Observed' else str(int(v))
        ax.text(j,v+max(vals)*0.015,label,ha='center',va='bottom',fontsize=8.5,fontweight='bold')
    for city,idx,col in [('Mardan',[0,1],'#1F4E79'),('Peshawar',[2,3],'#C55A11'),('Swat',[4,5],'#375623')]:
        mv=np.mean([vals[i] for i in idx])
        ax.plot([idx[0]-0.35,idx[1]+0.35],[mv,mv],color=col,lw=2,ls='--',zorder=5)
    ax.set_xticks(range(6)); ax.set_xticklabels(SAMPLES,fontsize=8.5,fontweight='bold',rotation=30,ha='right')
    ax.set_ylabel(ylab,fontsize=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.grid(axis='y',alpha=0.25); ax.tick_params(axis='x',length=0)

patches6=[mpatches.Patch(color='#1F4E79',label='Mardan'),
           mpatches.Patch(color='#C55A11',label='Peshawar'),
           mpatches.Patch(color='#375623',label='Swat'),
           mlines.Line2D([],[],color='gray',lw=2,ls='--',label='City mean')]
fig.legend(handles=patches6,loc='lower center',ncol=4,fontsize=10,
           frameon=True,edgecolor='#ccc',bbox_to_anchor=(0.5,-0.04))
plt.tight_layout()
plt.savefig('Figure6_AlphaDiversity_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 6 saved.')
