"""
FIGURE 11 — Virome Composition and Diversity
Input : Slaughterhouse_virome_curated.xlsx  (sheet: 'Virome (all species)')
Output: Figure11_Virome_600dpi.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pandas as pd, numpy as np

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITY_COL={'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}
BAR_COLORS=['#1F4E79','#1F4E79','#C55A11','#C55A11','#375623','#375623']

vdf=pd.read_excel('Slaughterhouse_virome_curated.xlsx',sheet_name='Virome (all species)')
zoo_host=['Aeromonas phage DH1','Yersinia phage Rostov 43']+\
         [x for x in vdf['Viral species'] if 'Streptococcus satellite phage' in str(x)]

def get_cat(row):
    if row['Viral species']=='uncultured human fecal virus': return 'Human Faecal Marker'
    if row['Viral species'] in zoo_host:                    return 'Zoonotic-Host Phage'
    if row['Type']=='Bacteriophage':                        return 'Bacteriophage (other)'
    if row['Type']=='Insect virus':                         return 'Insect Virus'
    if 'Giant virus' in str(row['Type']):                   return 'Giant Virus (protist)'
    return 'Unclassified Virus'

vdf['Category']=vdf.apply(get_cat,axis=1)
CAT_ORDER=['Human Faecal Marker','Zoonotic-Host Phage','Bacteriophage (other)','Insect Virus','Giant Virus (protist)','Unclassified Virus']
CAT_COL={'Human Faecal Marker':'#C00000','Zoonotic-Host Phage':'#7030A0','Bacteriophage (other)':'#1F4E79',
         'Insect Virus':'#70AD47','Giant Virus (protist)':'#FFC000','Unclassified Virus':'#BFBFBF'}
cat_sums=vdf.groupby('Category')[SAMPLES].sum()
RICHNESS=[19,10,23,9,21,8]
SHANNON=[0.956,0.901,0.934,1.186,0.907,0.232]

fig=plt.figure(figsize=(20,16)); fig.patch.set_facecolor('white')
ax_bar =fig.add_axes([0.05,0.52,0.50,0.44])
ax_lega=fig.add_axes([0.57,0.52,0.41,0.44])
ax_rich=fig.add_axes([0.05,0.05,0.22,0.40])
ax_shan=fig.add_axes([0.31,0.05,0.22,0.40])
ax_zoo =fig.add_axes([0.57,0.05,0.41,0.40])

x=np.arange(6)
# Panel A — stacked bar
ax_bar.set_facecolor('#FAFAFA')
bottoms=np.zeros(6)
for cat in CAT_ORDER:
    vals=cat_sums.loc[cat,SAMPLES].values.astype(float) if cat in cat_sums.index else np.zeros(6)
    ax_bar.bar(x,vals,bottom=bottoms,color=CAT_COL[cat],width=0.62,edgecolor='white',linewidth=0.6)
    for j,(v,b) in enumerate(zip(vals,bottoms)):
        if v>=0.018:
            tc='white' if cat in ['Human Faecal Marker','Zoonotic-Host Phage','Bacteriophage (other)'] else '#212529'
            ax_bar.text(x[j],b+v/2,f'{v*100:.1f}%',ha='center',va='center',fontsize=7,fontweight='bold',color=tc)
    bottoms+=vals
for j in range(6):
    ax_bar.text(x[j],bottoms[j]+0.005,f'{bottoms[j]*100:.1f}%',ha='center',va='bottom',fontsize=9,fontweight='bold')
bracket_y=1.07
for city,idx in [('Mardan',[0,1]),('Peshawar',[2,3]),('Swat',[4,5])]:
    xs=idx[0]-0.38; xe=idx[1]+0.38
    ax_bar.annotate('',xy=(xe,bracket_y),xytext=(xs,bracket_y),arrowprops=dict(arrowstyle='-',color=CITY_COL[city],lw=2))
    ax_bar.text((xs+xe)/2,bracket_y+0.01,city,ha='center',va='bottom',fontsize=11,fontweight='bold',color=CITY_COL[city])
ax_bar.set_xticks(x); ax_bar.set_xticklabels(SAMPLES,fontsize=10,fontweight='bold')
ax_bar.set_ylabel('Relative Abundance\n(within viral fraction)',fontsize=10)
ax_bar.set_xlabel('Sample',fontsize=10)
ax_bar.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1,decimals=0))
ax_bar.set_ylim(0,1.17); ax_bar.spines['top'].set_visible(False); ax_bar.spines['right'].set_visible(False)
ax_bar.grid(axis='y',alpha=0.2,zorder=0)
ax_bar.set_title('A.  Viral Community Composition by Functional Category',fontsize=11,fontweight='bold',pad=10,loc='left')

# Legend
ax_lega.axis('off'); ax_lega.set_xlim(0,1); ax_lega.set_ylim(0,1)
ax_lega.add_patch(mpatches.FancyBboxPatch((0,0),1,1,boxstyle='round,pad=0.02',fc='#F9F9F9',ec='#BBBBBB',lw=1))
ax_lega.text(0.5,0.97,'Legend — Viral Categories',ha='center',va='top',fontsize=11,fontweight='bold')
ax_lega.plot([0.03,0.97],[0.90,0.90],color='#BBBBBB',lw=0.8)
leg_items=[('Human Faecal Marker','#C00000','Uncultured human faecal virus\n(6/6 samples; mean 21.84%)\nMixed human-animal effluent marker'),
           ('Zoonotic-Host Phage','#7030A0','Phages of zoonotic bacterial hosts:\nAeromonas phage DH1 (PSLW1)\nYersinia phage Rostov 43 (MSLW2)\n12× Streptococcus satellite phages'),
           ('Bacteriophage (other)','#1F4E79','General bacteriophages\n(dominant; 62.9–95.4% per sample)'),
           ('Insect Virus','#70AD47','Baculoviridae (Naldaviricetes)\nEnvironmental insect-derived input'),
           ('Giant Virus (protist)','#FFC000','Faustovirus, MELD virus,\nPithovirus sibericum'),
           ('Unclassified Virus','#BFBFBF','Unclassified viral sequences')]
ys=np.linspace(0.84,0.05,6)
for (label,col,desc),y in zip(leg_items,ys):
    ax_lega.add_patch(mpatches.Rectangle((0.03,y-0.015),0.07,0.055,fc=col,ec='white',lw=0.5))
    ax_lega.text(0.13,y+0.012,label,fontsize=9,fontweight='bold',color=col,va='center')
    ax_lega.text(0.13,y-0.022,desc,fontsize=7.5,color='#444',va='top',linespacing=1.3)

# Panel B richness
ax_rich.set_facecolor('#FAFAFA')
ax_rich.bar(x,RICHNESS,color=BAR_COLORS,width=0.62,edgecolor='white',linewidth=0.6)
for j,v in enumerate(RICHNESS): ax_rich.text(j,v+0.4,str(int(v)),ha='center',va='bottom',fontsize=9.5,fontweight='bold')
ax_rich.set_xticks(x); ax_rich.set_xticklabels(SAMPLES,fontsize=9,fontweight='bold')
ax_rich.set_ylabel('Viral Species Richness',fontsize=10)
ax_rich.spines['top'].set_visible(False); ax_rich.spines['right'].set_visible(False)
ax_rich.set_ylim(0,28); ax_rich.grid(axis='y',alpha=0.25)
for city,idx,col in [('Mardan',[0,1],'#1F4E79'),('Peshawar',[2,3],'#C55A11'),('Swat',[4,5],'#375623')]:
    mv=np.mean([RICHNESS[i] for i in idx])
    ax_rich.plot([idx[0]-0.35,idx[1]+0.35],[mv,mv],color=col,lw=2,ls='--',zorder=5)
ax_rich.set_title('B.  Viral Species Richness',fontsize=10,fontweight='bold',pad=8,loc='left')

# Panel C shannon
ax_shan.set_facecolor('#FAFAFA')
ax_shan.bar(x,SHANNON,color=BAR_COLORS,width=0.62,edgecolor='white',linewidth=0.6)
for j,v in enumerate(SHANNON): ax_shan.text(j,v+0.02,f'{v:.3f}',ha='center',va='bottom',fontsize=9.5,fontweight='bold')
ax_shan.set_xticks(x); ax_shan.set_xticklabels(SAMPLES,fontsize=9,fontweight='bold')
ax_shan.set_ylabel("Shannon Entropy (H')",fontsize=10)
ax_shan.spines['top'].set_visible(False); ax_shan.spines['right'].set_visible(False)
ax_shan.set_ylim(0,1.50); ax_shan.grid(axis='y',alpha=0.25)
for city,idx,col in [('Mardan',[0,1],'#1F4E79'),('Peshawar',[2,3],'#C55A11'),('Swat',[4,5],'#375623')]:
    mv=np.mean([SHANNON[i] for i in idx])
    ax_shan.plot([idx[0]-0.35,idx[1]+0.35],[mv,mv],color=col,lw=2,ls='--',zorder=5)
ax_shan.set_title("C.  Viral Alpha Diversity (Shannon H')",fontsize=10,fontweight='bold',pad=8,loc='left')

# Panel D zoonotic phages — evenly spaced boxes
ax_zoo.axis('off'); ax_zoo.set_xlim(0,1); ax_zoo.set_ylim(0,1)
ax_zoo.set_title('D.  Zoonotic-Host Phages Detected',fontsize=10,fontweight='bold',pad=8,loc='left')
zoo_items=[('Aeromonas phage DH1','#ED7D31','Host: Aeromonas (opportunistic zoonotic)\nDetected: PSLW1 (Peshawar)\nAbundance: 0.0681% of viral fraction'),
           ('Yersinia phage Rostov 43','#C00000','Host: Yersinia enterocolitica (confirmed zoonotic)\nDetected: MSLW2 (Mardan)\nAbundance: 0.0439% of viral fraction'),
           ('Streptococcus satellite\nphages (12 variants)','#7030A0','Host: Streptococcus (incl. S. suis)\nDetected: MSLW1, PSLW1, SSLW2\nCombined mean: 0.463% of viral fraction'),
           ('Human faecal virus\n(contamination marker)','#C00000','Detected: 6/6 samples\nMean abundance: 21.84% of viral fraction\nHuman-animal mixed effluent confirmed')]
box_h=0.215; gap=(1.0-4*box_h)/5
tops=[1.0-gap-i*(box_h+gap) for i in range(4)]
for (label,col,desc),top in zip(zoo_items,tops):
    ax_zoo.add_patch(mpatches.FancyBboxPatch((0.01,top-box_h),0.98,box_h,boxstyle='round,pad=0.01',fc=col+'18',ec=col,lw=1.5))
    ax_zoo.text(0.04,top-0.028,label,fontsize=9,fontweight='bold',color=col,va='top',linespacing=1.3)
    ax_zoo.text(0.04,top-0.095,desc,fontsize=8,color='#333333',va='top',linespacing=1.4)

plt.savefig('Figure11_Virome_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 11 saved.')
