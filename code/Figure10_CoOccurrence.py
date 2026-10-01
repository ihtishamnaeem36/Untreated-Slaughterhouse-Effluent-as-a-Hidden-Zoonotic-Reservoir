"""
FIGURE 10 — Co-occurrence Network of Zoonotic Genera
Input : merged_genus_relative_abundance.csv
Output: Figure10_CoOccurrence_600dpi.png
Dependencies: pip install adjustText scipy
Verified network stats: 104 edges (65 positive, 39 negative), density 34.7%
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd, numpy as np
from scipy.stats import spearmanr
from adjustText import adjust_text

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
df_s=df.sort_values('mean',ascending=False).reset_index(drop=True)

n=len(ZOONOTIC_ALL)
abund_m={g:df_s[df_s['name']==g].iloc[0][SAMPLES].values.astype(float)
         if len(df_s[df_s['name']==g])>0 else np.zeros(6) for g in ZOONOTIC_ALL}

# Spearman correlation matrix
corr_mat=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i!=j:
            r,_=spearmanr(abund_m[ZOONOTIC_ALL[i]],abund_m[ZOONOTIC_ALL[j]])
            corr_mat[i,j]=r

# Edges at |rho| > 0.6
edges_pos,edges_neg=[],[]
for i in range(n):
    for j in range(i+1,n):
        r=corr_mat[i,j]
        if abs(r)>0.6: (edges_pos if r>0 else edges_neg).append((i,j,r))

n_pos=len(edges_pos); n_neg=len(edges_neg); n_total=n_pos+n_neg
density=n_total/(n*(n-1)/2)*100
print(f'Network: {n_pos} pos + {n_neg} neg = {n_total} total | density={density:.1f}%')

# Fruchterman-Reingold layout
np.random.seed(42)
pos={i:np.array([np.cos(2*np.pi*i/n)*5,np.sin(2*np.pi*i/n)*5]) for i in range(n)}
for iteration in range(400):
    disp={i:np.zeros(2) for i in range(n)}; k=1.5/np.sqrt(n)
    for i in range(n):
        for j in range(n):
            if i!=j:
                d=pos[i]-pos[j]; dist=max(np.linalg.norm(d),0.01)
                disp[i]+=(d/dist)*(k**2/dist)
    for i,j,r in edges_pos+edges_neg:
        d=pos[i]-pos[j]; dist=max(np.linalg.norm(d),0.01); f=dist**2/k
        disp[i]-=(d/dist)*f; disp[j]+=(d/dist)*f
    temp=max(0.05,0.15*(1-iteration/400))
    for i in range(n):
        dm=np.linalg.norm(disp[i])
        if dm>0: pos[i]+=(disp[i]/dm)*min(dm,temp)
all_p=np.array(list(pos.values())); mn,mx=all_p.min(0),all_p.max(0)
for i in range(n): pos[i]=(pos[i]-mn)/(mx-mn+1e-9)*18-9

mean_abunds={g:float(df_s[df_s['name']==g]['mean'].values[0]) if len(df_s[df_s['name']==g])>0 else 1e-6 for g in ZOONOTIC_ALL}
degrees={i:sum(1 for e in edges_pos+edges_neg if i in e[:2]) for i in range(n)}

fig,ax=plt.subplots(figsize=(18,15)); fig.patch.set_facecolor('white'); ax.set_facecolor('#FAFAFA')
for i,j,r in edges_neg:
    x1,y1=pos[i]; x2,y2=pos[j]
    ax.plot([x1,x2],[y1,y2],color='#ED7D31',lw=min(abs(r)*3,3),alpha=0.6,ls='--',zorder=1)
for i,j,r in edges_pos:
    x1,y1=pos[i]; x2,y2=pos[j]; sig=abs(r)>=0.829
    ax.plot([x1,x2],[y1,y2],color='#375623' if sig else '#70AD47',lw=min(r*4,4),alpha=0.75 if sig else 0.4,zorder=1)
for i,g in enumerate(ZOONOTIC_ALL):
    x,y=pos[i]; s=clin(g)
    sz=max(300,np.log10(mean_abunds[g]*100+0.01)*700+900); hub=degrees[i]>=14
    ax.scatter(x,y,s=sz,c=CLIN_COL[s],edgecolors='#212529' if hub else 'white',linewidths=2.8 if hub else 1.0,zorder=3,alpha=0.92)

texts=[]
for i,g in enumerate(ZOONOTIC_ALL):
    x,y=pos[i]; s=clin(g); hub=degrees[i]>=14; prefix=CLIN_SYM.get(s,'  ')
    t=ax.text(x,y,f'{prefix}{g}',ha='center',va='center',fontsize=9 if hub else 8,
              fontweight='bold' if hub else 'normal',color=CLIN_COL[s],style='italic',zorder=5,
              bbox=dict(boxstyle='round,pad=0.18',fc='white',ec='none',alpha=0.8))
    texts.append(t)
try:
    adjust_text(texts,ax=ax,expand_points=(2.0,2.2),expand_text=(1.6,2.0),
                force_points=0.4,arrowprops=dict(arrowstyle='-',color='#aaa',lw=0.5))
except: pass

ax.set_xlim(-13,13); ax.set_ylim(-13,13); ax.axis('off')
leg=[mpatches.Patch(color='none',label='── NODE COLOUR ──'),
     mpatches.Patch(color='#C00000',label='▶  Zoonotic Pathogen (WHO/WOAH/FAO/CDC)'),
     mpatches.Patch(color='#7030A0',label='★  WHO Critical-Priority MDR Organism'),
     mpatches.Patch(color='#ED7D31',label='◆  Opportunistic Zoonotic Pathogen'),
     mpatches.Patch(color='none',label='── NODE SIZE ──')]
for lab,sz in [('Low abundance (<0.01%)',250),('Moderate (~0.2%)',700),('High (~4.5%)',1400)]:
    leg.append(mlines.Line2D([],[],marker='o',linestyle='none',markersize=np.sqrt(sz)*0.52,color='#888',label=f'  {lab}'))
leg+=[mpatches.Patch(color='none',label='── NODE BORDER ──'),
      mlines.Line2D([],[],marker='o',linestyle='none',markersize=13,color='#888',
                    markeredgecolor='#212529',markeredgewidth=2.5,label='  Hub node (degree ≥ 14)'),
      mlines.Line2D([],[],marker='o',linestyle='none',markersize=10,color='#888',
                    markeredgecolor='white',markeredgewidth=1,label='  Standard node'),
      mpatches.Patch(color='none',label='── EDGE TYPE ──'),
      mlines.Line2D([],[],color='#375623',lw=3,label='  Positive co-occurrence, p<0.05  (|ρ|≥0.829)'),
      mlines.Line2D([],[],color='#70AD47',lw=1.5,label='  Positive co-occurrence, p>0.05  (0.6<|ρ|<0.829)'),
      mlines.Line2D([],[],color='#ED7D31',lw=2,ls='--',label='  Mutual exclusion  (ρ < −0.6)'),
      mpatches.Patch(color='none',
          label=f'{n_total} total edges  |  {n_pos} positive  |  {n_neg} negative  |  Density={density:.1f}%  |  Threshold |ρ|>0.6  |  n=6')]
ax.legend(handles=leg,loc='upper right',bbox_to_anchor=(1.0,1.0),fontsize=8.5,frameon=True,
          framealpha=0.97,edgecolor='#ccc',borderpad=1.0,labelspacing=0.42,
          title='Legend — Clinical Classification & Network Parameters',title_fontsize=9.5)
plt.tight_layout()
plt.savefig('Figure10_CoOccurrence_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 10 saved.')
