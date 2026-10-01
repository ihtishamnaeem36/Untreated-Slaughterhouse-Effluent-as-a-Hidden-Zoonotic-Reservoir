"""
FIGURE 7 — PCoA of Bray-Curtis Dissimilarity
Input : merged_genus_relative_abundance.csv
Output: Figure7_PCoA_600dpi.png
NOTE  : Bray-Curtis calculated from raw genus abundances; PCoA via classical MDS
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pandas as pd, numpy as np
from scipy.spatial.distance import braycurtis
from matplotlib.patches import Ellipse

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITIES=['Mardan','Mardan','Peshawar','Peshawar','Swat','Swat']
CITY_COL={'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}
MARKERS=['o','s','o','s','o','s']

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
mat=df[SAMPLES].values.T  # shape (6, n_genera)

# Bray-Curtis matrix
bc_mat=np.zeros((6,6))
for i in range(6):
    for j in range(6):
        bc_mat[i,j]=braycurtis(mat[i],mat[j])

# PCoA via classical MDS
A=-0.5*(bc_mat**2)
G=A-A.mean(axis=1,keepdims=True)-A.mean(axis=0,keepdims=True)+A.mean()
eigvals,eigvecs=np.linalg.eigh(G)
idx=np.argsort(eigvals)[::-1]; eigvals=eigvals[idx]; eigvecs=eigvecs[:,idx]
pos_var=eigvals[eigvals>0].sum()
pct1=eigvals[0]/pos_var*100; pct2=eigvals[1]/pos_var*100
pc1=eigvecs[:,0]*np.sqrt(max(eigvals[0],0))
pc2=eigvecs[:,1]*np.sqrt(max(eigvals[1],0))
sample_colors=[CITY_COL[c] for c in CITIES]

fig=plt.figure(figsize=(14,9)); fig.patch.set_facecolor('white')
ax    =fig.add_axes([0.07,0.09,0.60,0.87])
ax_leg=fig.add_axes([0.70,0.09,0.29,0.87])
ax.set_facecolor('#FAFAFA'); ax_leg.axis('off')

for city,sidx in [('Mardan',[0,1]),('Peshawar',[2,3]),('Swat',[4,5])]:
    cx=np.mean([pc1[i] for i in sidx]); cy=np.mean([pc2[i] for i in sidx])
    rx=max(abs(pc1[i]-cx) for i in sidx)+0.038
    ry=max(abs(pc2[i]-cy) for i in sidx)+0.038
    ax.add_patch(Ellipse((cx,cy),width=rx*2,height=ry*2,
                 facecolor=CITY_COL[city],alpha=0.10,edgecolor=CITY_COL[city],
                 linewidth=1.5,linestyle='--',zorder=1))
ax.axhline(0,color='#BBBBBB',lw=0.8,ls='--',zorder=0)
ax.axvline(0,color='#BBBBBB',lw=0.8,ls='--',zorder=0)
for i in range(6):
    ax.scatter(pc1[i],pc2[i],c=sample_colors[i],marker=MARKERS[i],
               s=250,edgecolors='white',linewidths=1.8,zorder=4)

offsets={'MSLW1':(-0.055,-0.025),'MSLW2':(-0.055,0.025),'PSLW1':(0.010,-0.032),
         'PSLW2':(-0.082,0.012),'SSLW1':(-0.062,0.000),'SSLW2':(0.030,-0.025)}
for i,s in enumerate(SAMPLES):
    ox,oy=offsets[s]
    ax.annotate(s,xy=(pc1[i],pc2[i]),xytext=(pc1[i]+ox,pc2[i]+oy),
                fontsize=10,fontweight='bold',color=sample_colors[i],
                ha='center',va='center',zorder=5,
                arrowprops=dict(arrowstyle='-',color=sample_colors[i],lw=0.9))

ax.set_xlabel(f'PC1  ({pct1:.1f}% of variance explained)',fontsize=12,labelpad=10)
ax.set_ylabel(f'PC2  ({pct2:.1f}% of variance explained)',fontsize=12,labelpad=10)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(alpha=0.15,zorder=0)
ax.set_xlim(pc1.min()-0.09,pc1.max()+0.09)
ax.set_ylim(pc2.min()-0.09,pc2.max()+0.09)
ax.annotate('SSLW2: high\nAcetoanaerobium\n(27.92%)',xy=(pc1[5],pc2[5]),
            xytext=(pc1[5]-0.13,pc2[5]-0.10),fontsize=8.5,color='#375623',style='italic',ha='center',
            arrowprops=dict(arrowstyle='->',color='#375623',lw=1.0),
            bbox=dict(boxstyle='round,pad=0.3',fc='#E2EFDA',ec='#375623',lw=0.8))
ax.annotate('MSLW1 & SSLW1\noverlap (BC=0.100)',
            xy=((pc1[0]+pc1[4])/2,(pc2[0]+pc2[4])/2),
            xytext=(pc1[0]-0.15,pc2[0]-0.10),fontsize=8.5,color='#1F4E79',style='italic',ha='center',
            arrowprops=dict(arrowstyle='->',color='#666',lw=1.0),
            bbox=dict(boxstyle='round,pad=0.3',fc='#D6E4F0',ec='#1F4E79',lw=0.8))

# Legend panel
ax_leg.text(0.05,0.98,'City',fontsize=11,fontweight='bold',color='#212529',va='top',transform=ax_leg.transAxes)
city_p=[mpatches.Patch(color='#1F4E79',label='Mardan'),mpatches.Patch(color='#C55A11',label='Peshawar'),mpatches.Patch(color='#375623',label='Swat')]
leg1=ax_leg.legend(handles=city_p,loc='upper left',bbox_to_anchor=(0.0,0.93),fontsize=10,frameon=True,framealpha=0.95,edgecolor='#cccccc',borderpad=0.8,labelspacing=0.5)
ax_leg.add_artist(leg1)
ax_leg.text(0.05,0.68,'Sample / Symbol',fontsize=11,fontweight='bold',color='#212529',va='top',transform=ax_leg.transAxes)
mh=[mlines.Line2D([],[],marker='o',linestyle='none',color='#555',markersize=10,label='Replicate 1\n(MSLW1/PSLW1/SSLW1)'),
    mlines.Line2D([],[],marker='s',linestyle='none',color='#555',markersize=10,label='Replicate 2\n(MSLW2/PSLW2/SSLW2)'),
    mpatches.Patch(facecolor='#888',alpha=0.18,edgecolor='#888',linestyle='--',label='City grouping ellipse')]
leg2=ax_leg.legend(handles=mh,loc='upper left',bbox_to_anchor=(0.0,0.63),fontsize=9.5,frameon=True,framealpha=0.95,edgecolor='#cccccc',borderpad=0.8,labelspacing=0.7)
ax_leg.add_artist(leg2)
ax_leg.text(0.05,0.33,'Bray-Curtis Dissimilarity',fontsize=11,fontweight='bold',color='#212529',va='top',transform=ax_leg.transAxes)
table_rows=[('Mardan within','0.243','#375623'),('Peshawar within','0.767','#C00000'),
            ('Swat within','0.780','#C00000'),('Mardan–Peshawar','0.720','#C00000'),
            ('Mardan–Swat','0.476','#FFC000'),('Peshawar–Swat','0.661','#C00000')]
for k,(label,val,col) in enumerate(table_rows):
    y=0.27-k*0.053
    ax_leg.text(0.05,y,label+':',fontsize=8.5,color='#212529',va='top',transform=ax_leg.transAxes)
    ax_leg.text(0.82,y,val,fontsize=8.5,fontweight='bold',color=col,va='top',ha='right',transform=ax_leg.transAxes)
ax_leg.text(0.05,-0.01,f'PC1 + PC2 = {pct1+pct2:.1f}% variance explained',fontsize=9,color='#2E75B6',va='top',transform=ax_leg.transAxes,
            bbox=dict(boxstyle='round,pad=0.35',fc='#EBF3FB',ec='#2E75B6',lw=0.8))

plt.savefig('Figure7_PCoA_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 7 saved.')
