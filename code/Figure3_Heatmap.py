"""
FIGURE 3 — Top 30 Genera Heatmap (Row-Normalised, Reddish Colourmap)
Input : merged_genus_relative_abundance.csv
Output: Figure3_Heatmap_600dpi.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd, numpy as np
from matplotlib.colors import LinearSegmentedColormap

DPI=400  # 400 DPI used to avoid memory limits on large figure; still publication quality
plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.6})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CONFIRMED={'Arcobacter','Aliarcobacter','Salmonella','Brucella','Campylobacter','Listeria','Yersinia','Leptospira','Helicobacter','Treponema'}
CRITICAL={'Acinetobacter','Pseudomonas'}
OPPORTUNISTIC={'Aeromonas','Escherichia','Klebsiella','Clostridium','Enterococcus','Streptococcus','Vibrio','Fusobacterium','Shewanella','Corynebacterium','Staphylococcus','Stenotrophomonas','Bacteroides'}

def clin(name):
    if name in CONFIRMED: return 'zoonotic'
    if name in CRITICAL:  return 'critical'
    if name in OPPORTUNISTIC: return 'opportunistic'
    return 'environmental'

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
df['mean']=df[SAMPLES].mean(axis=1)
df_s=df.sort_values('mean',ascending=False).reset_index(drop=True)
top30=df_s.head(30)[['name']+SAMPLES].copy()
genera=list(top30['name'])
mat=top30[SAMPLES].values.astype(float)
row_max=mat.max(axis=1,keepdims=True); row_max[row_max==0]=1
mat_norm=mat/row_max
cmap_rd=LinearSegmentedColormap.from_list('rd',['#FFFFFF','#FCBBA1','#FB6A4A','#CB181D','#67000D'],N=256)

fig=plt.figure(figsize=(15,24)); fig.patch.set_facecolor('white')
ax_city=fig.add_axes([0.13,0.944,0.58,0.026])
ax_heat=fig.add_axes([0.13,0.100,0.58,0.838])
ax_cb  =fig.add_axes([0.13,0.060,0.58,0.016])
ax_leg =fig.add_axes([0.74,0.510,0.24,0.430])

# City strip
ax_city.set_xlim(0,6); ax_city.set_ylim(0,1); ax_city.axis('off')
for s,e,col,label in [(0,2,'#1F4E79','Mardan'),(2,4,'#C55A11','Peshawar'),(4,6,'#375623','Swat')]:
    ax_city.add_patch(mpatches.FancyBboxPatch((s+0.04,0.06),e-s-0.08,0.88,boxstyle='round,pad=0.02',fc=col,ec='white',lw=1.2))
    ax_city.text((s+e)/2,0.52,label,ha='center',va='center',fontsize=11,fontweight='bold',color='white')

# Heatmap
im=ax_heat.imshow(mat_norm,aspect='auto',cmap=cmap_rd,vmin=0,vmax=1,interpolation='nearest')
for i in range(31): ax_heat.axhline(i-0.5,color='white',lw=0.7)
for j in range(7):  ax_heat.axvline(j-0.5,color='white',lw=0.7)
ax_heat.axvline(1.5,color='#555',lw=1.5,ls='--',alpha=0.7)
ax_heat.axvline(3.5,color='#555',lw=1.5,ls='--',alpha=0.7)

# Cell values only where raw >= 0.5%
for i in range(30):
    for j in range(6):
        v_raw=mat[i,j]; v_norm=mat_norm[i,j]
        if v_raw>=0.005:
            tc='white' if v_norm>0.60 else '#212529'
            ax_heat.text(j,i,f'{v_raw*100:.2f}%',ha='center',va='center',fontsize=6.5,color=tc,
                         fontweight='bold' if v_norm>0.7 else 'normal')

ax_heat.set_yticks(range(30))
y_labels,y_colors,y_weights=[],[],[]
for g in genera:
    s=clin(g)
    if s=='zoonotic':      y_labels.append(f'▶  {g}'); y_colors.append('#C00000'); y_weights.append('bold')
    elif s=='critical':    y_labels.append(f'★  {g}'); y_colors.append('#7030A0'); y_weights.append('bold')
    elif s=='opportunistic': y_labels.append(f'◆  {g}'); y_colors.append('#C55A11'); y_weights.append('bold')
    else:                  y_labels.append(f'    {g}'); y_colors.append('#212529'); y_weights.append('normal')
ax_heat.set_yticklabels(y_labels,fontsize=9.5,style='italic')
for lbl,col,wt in zip(ax_heat.get_yticklabels(),y_colors,y_weights):
    lbl.set_color(col); lbl.set_fontweight(wt)
ax_heat.set_xticks(range(6))
ax_heat.set_xticklabels(['MSLW1\n(Mardan)','MSLW2\n(Mardan)','PSLW1\n(Peshawar)',
                          'PSLW2\n(Peshawar)','SSLW1\n(Swat)','SSLW2\n(Swat)'],
                         fontsize=9.5,fontweight='bold',linespacing=1.5)
ax_heat.tick_params(axis='x',length=0,pad=10); ax_heat.tick_params(axis='y',length=0,pad=5)
ax_heat.xaxis.set_ticks_position('bottom')

# Colorbar
cb=fig.colorbar(im,cax=ax_cb,orientation='horizontal')
cb.set_label('Row-Normalised Relative Abundance  (0 = lowest within genus  |  1 = highest within genus)',fontsize=9,labelpad=7)
cb.ax.tick_params(labelsize=8.5); cb.set_ticks([0,0.2,0.4,0.6,0.8,1.0])

# Legend — top right separate panel
ax_leg.axis('off'); ax_leg.set_xlim(0,1); ax_leg.set_ylim(0,1)
ax_leg.add_patch(mpatches.FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',fc='#F9F9F9',ec='#CCCCCC',lw=1.2))
ax_leg.text(0.5,0.95,'Clinical Classification\nLegend',ha='center',va='top',fontsize=10,fontweight='bold',color='#212529',linespacing=1.5)
ax_leg.plot([0.05,0.95],[0.83,0.83],color='#CCCCCC',lw=0.8)
items=[('▶','#C00000','Zoonotic Pathogen','WHO / WOAH / FAO / CDC listed'),
       ('★','#7030A0','WHO Critical-Priority MDR','Pseudomonas  |  Acinetobacter'),
       ('◆','#C55A11','Opportunistic Zoonotic','AMR or opportunistic risk'),
       ('  ','#444444','Environmental /\nNon-clinical','No zoonotic transmission\nconfirmed')]
for (sym,col,title,desc),y in zip(items,[0.72,0.52,0.32,0.10]):
    ax_leg.text(0.06,y+0.04,sym,ha='left',va='center',fontsize=14,color=col,fontweight='bold')
    ax_leg.text(0.20,y+0.04,title,ha='left',va='center',fontsize=9,color=col,fontweight='bold',style='italic',linespacing=1.3)
    ax_leg.text(0.20,y-0.07,desc,ha='left',va='center',fontsize=7.8,color='#555555',linespacing=1.3)

plt.savefig('Figure3_Heatmap_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 3 saved.')
