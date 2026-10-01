"""
FIGURE ZHI — Zoonotic Hazard Index (Sample-Level + City Risk Ranking)
Input : No CSV needed — ZHI values pre-calculated from genus data
Output: Figure_ZHI_600dpi.png
ZHI formula: [Σ(A × PS × RS × TM)] × (N_detected/25) × 100
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITIES=['Mardan','Mardan','Peshawar','Peshawar','Swat','Swat']
ZHI=[50.51,100.91,192.93,68.76,54.72,106.76]
CITY_COL={'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}
BAR_COLORS=[CITY_COL[c] for c in CITIES]
CITY_MEANS={'Peshawar':130.84,'Swat':80.74,'Mardan':75.71}
CITY_ORDER=['Peshawar','Swat','Mardan']

fig,axes=plt.subplots(1,2,figsize=(15,7)); fig.patch.set_facecolor('white')

# Panel A — sample bars
ax1=axes[0]; ax1.set_facecolor('white')
bars=ax1.bar(range(6),ZHI,color=BAR_COLORS,width=0.6,edgecolor='white',linewidth=0.8,zorder=3)
for b,v in zip(bars,ZHI):
    ax1.text(b.get_x()+b.get_width()/2,v+2.5,f'{v:.1f}',ha='center',va='bottom',fontsize=9.5,fontweight='bold',color='#212529')
bracket_y=208
for city,idx in [('Mardan',[0,1]),('Peshawar',[2,3]),('Swat',[4,5])]:
    xs=idx[0]-0.38; xe=idx[1]+0.38
    ax1.annotate('',xy=(xe,bracket_y),xytext=(xs,bracket_y),arrowprops=dict(arrowstyle='-',color=CITY_COL[city],lw=2.2))
    ax1.text((xs+xe)/2,bracket_y+3.5,city,ha='center',va='bottom',fontsize=11,fontweight='bold',color=CITY_COL[city])
mean_zhi=np.mean(ZHI)
ax1.axhline(y=mean_zhi,color='#888',linestyle='--',linewidth=1.3,zorder=2)
ax1.text(5.42,mean_zhi+2,f'Mean={mean_zhi:.1f}',fontsize=8.5,color='#888',ha='right',va='bottom')
ax1.set_xticks(range(6)); ax1.set_xticklabels(SAMPLES,fontsize=10.5,fontweight='bold')
ax1.set_ylabel('Zoonotic Hazard Index (ZHI)',fontsize=11,labelpad=8); ax1.set_xlabel('Sample',fontsize=11,labelpad=8)
ax1.set_ylim(0,232); ax1.set_xlim(-0.55,5.55)
ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
ax1.grid(axis='y',alpha=0.25,zorder=1)
ax1.set_title('A.  Sample-Level ZHI Scores',fontsize=12,fontweight='bold',pad=12,loc='left')

# Panel B — city ranking
ax2=axes[1]; ax2.set_facecolor('white')
c_names=CITY_ORDER; c_vals=[CITY_MEANS[c] for c in CITY_ORDER]; c_cols=[CITY_COL[c] for c in CITY_ORDER]
ax2.axvspan(0,75,alpha=0.07,color='#70AD47',zorder=0)
ax2.axvspan(75,120,alpha=0.07,color='#FFC000',zorder=0)
ax2.axvspan(120,200,alpha=0.07,color='#FF0000',zorder=0)
hbars=ax2.barh(range(3),c_vals,color=c_cols,height=0.50,edgecolor='white',linewidth=0.8,zorder=3)
for i,(v,name) in enumerate(zip(c_vals,c_names)):
    ax2.text(v+2.5,i,f'ZHI = {v:.1f}  |  Rank #{i+1}',va='center',fontsize=10,fontweight='bold',color='#212529')
ax2.text(37.5,-0.70,'MODERATE\nRISK',ha='center',fontsize=8.5,color='#375623',fontweight='bold',va='top')
ax2.text(97.5,-0.70,'HIGH\nRISK',ha='center',fontsize=8.5,color='#7F6000',fontweight='bold',va='top')
ax2.text(160,-0.70,'VERY HIGH\nRISK',ha='center',fontsize=8.5,color='#C00000',fontweight='bold',va='top')
ax2.set_yticks(range(3)); ax2.set_yticklabels(c_names,fontsize=11,fontweight='bold')
for tick,col in zip(ax2.get_yticklabels(),c_cols): tick.set_color(col)
ax2.set_xlabel('Mean City Zoonotic Hazard Index (ZHI)',fontsize=11,labelpad=8)
ax2.set_xlim(0,200); ax2.set_ylim(-1.0,2.6)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.grid(axis='x',alpha=0.25,zorder=1)
ax2.set_title('B.  City-Level Risk Ranking',fontsize=12,fontweight='bold',pad=12,loc='left')

plt.tight_layout(rect=[0,0.04,1,1])
plt.savefig('Figure_ZHI_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('ZHI figure saved.')
