"""
FIGURE 8 — Bray-Curtis Dissimilarity Matrix Heatmap
Input : merged_genus_relative_abundance.csv
Output: Figure8_BrayCurtis_600dpi.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd, numpy as np
from scipy.spatial.distance import braycurtis
from matplotlib.colors import LinearSegmentedColormap

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITY_COL={'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
mat=df[SAMPLES].values.T
bc_mat=np.zeros((6,6))
for i in range(6):
    for j in range(6):
        bc_mat[i,j]=braycurtis(mat[i],mat[j])

cmap_bc=LinearSegmentedColormap.from_list('bc',['#C6EFCE','#FFEB9C','#FFC7CE'],N=256)
fig,ax=plt.subplots(figsize=(9,8)); fig.patch.set_facecolor('white')
im=ax.imshow(bc_mat,cmap=cmap_bc,vmin=0,vmax=1,aspect='auto')
for i in range(6):
    for j in range(6):
        v=bc_mat[i,j]; tc='white' if v>0.75 else '#212529'
        ax.text(j,i,f'{v:.3f}',ha='center',va='center',fontsize=11,fontweight='bold',color=tc)
ax.set_xticks(range(6)); ax.set_yticks(range(6))
ax.set_xticklabels(SAMPLES,fontsize=10,fontweight='bold')
ax.set_yticklabels(SAMPLES,fontsize=10,fontweight='bold')
for i in range(6):
    for j in range(6):
        ax.add_patch(mpatches.Rectangle((j-0.5,i-0.5),1,1,fill=False,ec='white',lw=1.5))
cb=fig.colorbar(im,ax=ax,orientation='vertical',fraction=0.046,pad=0.04)
cb.set_label('Bray-Curtis Dissimilarity\n(0 = identical  |  1 = completely dissimilar)',fontsize=9)
cb.set_ticks([0,0.25,0.5,0.75,1.0])
for idx,col in [([0,1],'#1F4E79'),([2,3],'#C55A11'),([4,5],'#375623')]:
    ax.add_patch(mpatches.Rectangle((idx[0]-0.5,idx[0]-0.5),2,2,fill=False,ec=col,lw=2.5,zorder=5))
notes=[mpatches.Patch(color='#C6EFCE',label='Low dissimilarity (<0.3)'),
       mpatches.Patch(color='#FFEB9C',label='Moderate (0.3–0.6)'),
       mpatches.Patch(color='#FFC7CE',label='High dissimilarity (>0.6)')]
ax.legend(handles=notes,loc='upper left',bbox_to_anchor=(1.18,1.0),fontsize=9,frameon=True,edgecolor='#ccc')
plt.tight_layout()
plt.savefig('Figure8_BrayCurtis_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 8 saved.')
