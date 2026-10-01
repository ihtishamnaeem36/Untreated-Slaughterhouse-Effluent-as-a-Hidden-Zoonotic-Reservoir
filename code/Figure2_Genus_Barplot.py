"""
FIGURE 2 — Top 20 Genus Barplot with Clinical Colour Coding
Input : merged_genus_relative_abundance.csv
Output: Figure2_Genus_Barplot_600dpi.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pandas as pd, numpy as np

DPI=600; plt.rcParams.update({'font.family':'DejaVu Sans','axes.linewidth':0.8})
SAMPLES=['MSLW1','MSLW2','PSLW1','PSLW2','SSLW1','SSLW2']
CITY_COL={'Mardan':'#1F4E79','Peshawar':'#C55A11','Swat':'#375623'}
CONFIRMED={'Arcobacter','Aliarcobacter','Salmonella','Brucella','Campylobacter','Listeria','Yersinia','Leptospira','Helicobacter','Treponema'}
CRITICAL={'Acinetobacter','Pseudomonas'}
OPPORTUNISTIC={'Aeromonas','Escherichia','Klebsiella','Clostridium','Enterococcus','Streptococcus','Vibrio','Fusobacterium','Shewanella','Corynebacterium','Staphylococcus','Stenotrophomonas','Bacteroides'}

def clin(name):
    if name in CONFIRMED: return 'zoonotic'
    if name in CRITICAL:  return 'critical'
    if name in OPPORTUNISTIC: return 'opportunistic'
    return 'environmental'

ENV_PAL=['#2C5F8A','#4488BB','#5BAABB','#548235','#C8860A','#7B9E3E','#843C0C',
         '#375623','#A0522D','#888888','#6D9EC1','#5D9E5D','#9090C0','#8FAF6F',
         '#909090','#2E7B8C','#666666']

df=pd.read_csv('merged_genus_relative_abundance.csv')
df=df[df['name'].str.lower()!='unclassified'].copy()
df['mean']=df[SAMPLES].mean(axis=1)
df_s=df.sort_values('mean',ascending=False).reset_index(drop=True)
top20=df_s.head(20)[['name']+SAMPLES].copy()
others_v=df_s.iloc[20:][SAMPLES].sum()
top20_plot=pd.concat([top20,pd.DataFrame([['Others']+list(others_v)],columns=['name']+SAMPLES)],ignore_index=True)

env_idx=0; bar_colors=[]
for _,row in top20_plot.iterrows():
    s=clin(row['name'])
    if s=='zoonotic': bar_colors.append('#C00000')
    elif s=='critical': bar_colors.append('#7030A0')
    elif s=='opportunistic': bar_colors.append('#ED7D31')
    else: bar_colors.append(ENV_PAL[env_idx%len(ENV_PAL)]); env_idx+=1

fig=plt.figure(figsize=(22,11)); fig.patch.set_facecolor('white')
ax=fig.add_axes([0.04,0.10,0.48,0.82])
ax_leg=fig.add_axes([0.54,0.04,0.44,0.93])
ax.set_facecolor('#FAFAFA'); ax_leg.axis('off')

x=np.arange(6); bottoms=np.zeros(6)
for i,row in top20_plot.iterrows():
    vals=row[SAMPLES].values.astype(float)
    ax.bar(x,vals,bottom=bottoms,color=bar_colors[i],width=0.62,edgecolor='white',linewidth=0.4,zorder=3)
    for j,(v,b) in enumerate(zip(vals,bottoms)):
        if v>=0.026:
            dark={'#C00000','#7030A0','#2C5F8A','#548235','#843C0C','#A0522D','#375623','#2E7B8C','#4488BB','#666666'}
            tc='white' if bar_colors[i] in dark else '#212529'
            ax.text(x[j],b+v/2,f'{v*100:.1f}%',ha='center',va='center',fontsize=6.5,fontweight='bold',color=tc,zorder=4)
    bottoms+=vals
for j,tot in enumerate(top20_plot[SAMPLES].sum()):
    ax.text(x[j],tot+0.007,f'{tot*100:.1f}%',ha='center',va='bottom',fontsize=10,fontweight='bold')
bracket_y=1.085
for city,idx in [('Mardan',[0,1]),('Peshawar',[2,3]),('Swat',[4,5])]:
    xs=idx[0]-0.38; xe=idx[1]+0.38
    ax.annotate('',xy=(xe,bracket_y),xytext=(xs,bracket_y),arrowprops=dict(arrowstyle='-',color=CITY_COL[city],lw=2.2))
    ax.text((xs+xe)/2,bracket_y+0.01,city,ha='center',va='bottom',fontsize=11,fontweight='bold',color=CITY_COL[city])
ax.set_xticks(x); ax.set_xticklabels(SAMPLES,fontsize=11,fontweight='bold')
ax.set_ylabel('Relative Abundance',fontsize=12,labelpad=8); ax.set_xlabel('Sample',fontsize=12,labelpad=8)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1,decimals=0))
ax.set_xlim(-0.55,5.55); ax.set_ylim(0,1.16)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.grid(axis='y',alpha=0.2,zorder=1)

lp=[mpatches.Patch(color='none',label='── ZOONOTIC PATHOGENS ──')]
for _,row in top20_plot.iterrows():
    name=row['name']; s=clin(name)
    if s=='zoonotic':
        mv=df_s[df_s['name']==name]['mean'].values[0]*100 if name in list(df_s['name']) else 0
        lp.append(mpatches.Patch(color='#C00000',label=f'  {name}  ({mv:.2f}%)'))
lp.append(mpatches.Patch(color='none',label='── WHO CRITICAL MDR ──'))
for _,row in top20_plot.iterrows():
    name=row['name']; s=clin(name)
    if s=='critical':
        mv=df_s[df_s['name']==name]['mean'].values[0]*100
        lp.append(mpatches.Patch(color='#7030A0',label=f'  {name}  ({mv:.2f}%)'))
lp.append(mpatches.Patch(color='none',label='── OPPORTUNISTIC ZOONOTIC ──'))
for _,row in top20_plot.iterrows():
    name=row['name']; s=clin(name)
    if s=='opportunistic':
        mv=df_s[df_s['name']==name]['mean'].values[0]*100
        lp.append(mpatches.Patch(color='#ED7D31',label=f'  {name}  ({mv:.2f}%)'))
lp.append(mpatches.Patch(color='none',label='── ENVIRONMENTAL ──'))
ec3=0
for _,row in top20_plot.iterrows():
    name=row['name']; s=clin(name)
    if s=='environmental':
        mv=df_s[df_s['name']==name]['mean'].values[0]*100 if name in list(df_s['name']) else 0
        col=ENV_PAL[ec3%len(ENV_PAL)]
        lbl=f'  {name}  ({mv:.2f}%)' if name!='Others' else '  Others (remaining genera)'
        lp.append(mpatches.Patch(color=col,label=lbl)); ec3+=1
ax_leg.legend(handles=lp,loc='upper left',fontsize=9.5,frameon=True,framealpha=0.95,
              edgecolor='#cccccc',borderpad=1.0,labelspacing=0.55,
              title='Genus | Clinical Classification',title_fontsize=10.5)
plt.savefig('Figure2_Genus_Barplot_600dpi.png',dpi=DPI,bbox_inches='tight',facecolor='white')
plt.close(); print('Figure 2 saved.')
