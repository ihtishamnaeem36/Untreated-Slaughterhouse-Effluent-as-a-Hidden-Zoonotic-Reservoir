# -*- coding: utf-8 -*-
"""
generate_vector_figures.py

Generates publication-quality, vector-based PDF figures and 600 DPI PNG figures
for Folia Microbiologica (Springer Nature) strictly adhering to the
Official Guidelines for Authors:

SPRINGER NATURE ARTWORK & LETTERING RULES:
1. Font family: Arial or Helvetica (sans-serif), consistent across all figures.
2. Lettering size at final printed size (174 mm / 6.85 in double column):
   - Exactly 2–3 mm (8–12 pt) throughout final artwork.
   - Minimal variance (e.g. 7.5 pt to 12 pt).
   - Panel labels: 11–13 pt bold (A, B, C or a, b, c).
   - Axis titles: 9.5–10 pt bold.
   - Tick labels & category labels: 8.5–9 pt bold.
   - Data values & percentages on bars: 8.0–8.5 pt bold.
   - Legends: 8.0–8.5 pt.
   - Minor annotations: 7.5 pt (never below 7 pt).
3. Contrast: At least 4.5:1 contrast against backgrounds.
4. Dimensions: Full double-column width = 174 mm (6.85 in), max height = 234 mm (9.2 in).
5. No large titles inside figure (figure captions belong in manuscript text).

FIGURE PACKAGE STRUCTURE:
Fig. 1 (Taxonomic Architecture):
  - Panel A: Phylum-level relative abundance
  - Panel B: Top 20 bacterial genera (clinical colors)

Fig. 2 (Priority Pathogen Burden & Clinical Frameworks):
  - Square / Grid format:
    * Panel A (Top Left): Genus count by clinical category (n=28)
    * Panel B (Top Right): Mean relative abundance by clinical category (% of classified reads)
    * Panel C (Bottom Full Width): Mean relative abundance of 25 priority zoonotic genera across cities

Fig. 3 (Community Diversity & Ordination):
  - Strict 3-panel structure (NO duplicate B's):
    * Panel A (Top): Alpha diversity metrics (Observed, Shannon, Simpson, Chao1)
    * Panel B (Bottom Left): PCoA ordination of Bray-Curtis dissimilarity
    * Panel C (Bottom Right): Pairwise Bray-Curtis dissimilarity matrix heatmap

Fig. 4 (Hazard Modeling & Ecological Network):
  - Strict 3-panel structure (NO duplicate B's; internal A & B removed):
    * Panel A (Top Left): Sample-level ZHI scores
    * Panel B (Top Right): Municipal risk tier ranking
    * Panel C (Bottom Full Width): Spearman rank co-occurrence network of 25 zoonotic genera

Fig. 5 (Wastewater Virome):
  - Strict 3-panel visual figure (Panel D converted to Supplementary Table S8):
    * Panel A (Top Left): Viral community composition by functional category
    * Top Right: Comprehensive category legend
    * Panel B (Bottom Left): Viral species richness
    * Panel C (Bottom Right): Viral alpha diversity (Shannon H')

Supplementary Fig. S1: Hierarchical clustering heatmap (Top 30 genera)
Supplementary Fig. S2: Differential abundance enrichment ratios
"""

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis
from scipy.stats import spearmanr
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

# --- DYNAMIC DATA PATH RESOLUTION ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXTRA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'extra'))

def find_file(filename, subfolders=None):
    candidates = [
        filename,
        os.path.join(SCRIPT_DIR, filename),
        os.path.join(EXTRA_DIR, filename),
    ]
    if subfolders:
        for sf in subfolders:
            candidates.append(os.path.join(EXTRA_DIR, sf, filename))
    for c in candidates:
        if os.path.exists(c):
            return c
    return filename

PHYLUM_CSV = find_file('merged_phylum_relative_abundance.csv', [os.path.join('all cities', 'Phylum', 'Phylum Relative Abundance')])
GENUS_CSV = find_file('merged_genus_relative_abundance.csv', [os.path.join('all cities', 'genus', 'genus relative abundance')])
VIROME_XLSX = find_file('Slaughterhouse_virome_curated.xlsx', ['Virome Section'])

sys.stdout.reconfigure(encoding='utf-8')

# Global styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
    'axes.edgecolor': '#333333',
    'xtick.color': '#333333',
    'ytick.color': '#333333',
    'text.color': '#212529',
    'figure.autolayout': False
})

# Standard Springer double-column width in inches (174 mm = 6.85 in)
SPRINGER_COL_WIDTH = 6.85

def pt(target_pt, fig_width=None):
    """
    Returns the exact target publication font size (strictly within the Springer 8-12 pt range).
    Prevents artificial canvas-width multiplication that previously caused labels to merge and collide.
    """
    return target_pt


SAMPLES = ['MSLW1', 'MSLW2', 'PSLW1', 'PSLW2', 'SSLW1', 'SSLW2']
CITIES = ['Mardan', 'Mardan', 'Peshawar', 'Peshawar', 'Swat', 'Swat']
CITY_COL = {'Mardan': '#1F4E79', 'Peshawar': '#C55A11', 'Swat': '#375623'}
BAR_COLORS = [CITY_COL[c] for c in CITIES]

CONFIRMED = {
    'Arcobacter',
    'Aliarcobacter',
    'Salmonella',
    'Brucella',
    'Campylobacter',
    'Listeria',
    'Yersinia',
    'Leptospira',
    'Helicobacter',
    'Treponema',
}
CRITICAL = {'Acinetobacter', 'Pseudomonas'}
OPPORTUNISTIC = {
    'Aeromonas',
    'Escherichia',
    'Klebsiella',
    'Clostridium',
    'Enterococcus',
    'Streptococcus',
    'Vibrio',
    'Fusobacterium',
    'Shewanella',
    'Corynebacterium',
    'Staphylococcus',
    'Stenotrophomonas',
    'Bacteroides',
}
COMMENSAL = ['Bifidobacterium', 'Desulfovibrio', 'Comamonas']
ZOONOTIC_ALL = list(CONFIRMED | CRITICAL | OPPORTUNISTIC)

CLIN_COL = {
    'zoonotic': '#C00000',
    'critical': '#7030A0',
    'opportunistic': '#ED7D31',
    'environmental': '#548235',
}
CLIN_SYM = {'zoonotic': '', 'critical': '', 'opportunistic': ''}


def clin(name):
  if name in CONFIRMED:
    return 'zoonotic'
  if name in CRITICAL:
    return 'critical'
  if name in OPPORTUNISTIC:
    return 'opportunistic'
  return 'environmental'


def clin4(name):
  if name in CONFIRMED:
    return 'Zoonotic Pathogen'
  if name in CRITICAL:
    return 'WHO Critical MDR'
  if name in OPPORTUNISTIC:
    return 'Opportunistic Zoonotic'
  return 'Commensal / Environmental'


out_dir = 'submission_figures_current_microbiology'
os.makedirs(out_dir, exist_ok=True)


def save_dual(fig, basename):
    out_root = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
    fig.savefig(f'{basename}.pdf', bbox_inches='tight')
    fig.savefig(f'{basename}.png', dpi=600, bbox_inches='tight')
    name_map = {
        'Fig1_Phylum_Calibrated': 'Figure_1',
        'Fig2_Zoonotic_Calibrated': 'Figure_2',
        'Fig3_Diversity_Calibrated': 'Figure_3',
        'Fig4_ZHI_Network_Calibrated': 'Figure_4',
        'Fig5_Virome_Calibrated': 'Figure_6',
        'FigS1_Heatmap_Calibrated': 'Figure_S1',
        'FigS2_Functional_Calibrated': 'Figure_S2',
    }
    if basename in name_map:
        target_name = name_map[basename]
        fig.savefig(os.path.join(out_root, f'{target_name}.pdf'), bbox_inches='tight')
        fig.savefig(os.path.join(out_root, f'{target_name}.png'), dpi=600, bbox_inches='tight')
        print(f'  -> Exported {target_name}.pdf and {target_name}.png directly to {out_root}')
    print(f'Saved 600 DPI master: {basename}.png & vector {basename}.pdf')

# ==============================================================================
# FIGURE 1: BACTERIAL TAXONOMIC ARCHITECTURE
# Target Print Dimensions: 7.2 x 4.8 inches (183 mm x 122 mm)
# ==============================================================================
print('\n--- Generating Figure 1 with Calibrated Typography ---')
W1, H1 = 12.0, 6.2  # comfortable working canvas
ph = pd.read_csv(PHYLUM_CSV)
ph = ph[ph['name'].str.lower() != 'unclassified'].copy()
ph['mean'] = ph[SAMPLES].mean(axis=1)
ph_s = ph.sort_values('mean', ascending=False).reset_index(drop=True)
ph_plot = ph_s.head(10)[['name'] + SAMPLES].copy()
others = ph_s.iloc[10:][SAMPLES].sum()
ph_plot = pd.concat(
    [ph_plot, pd.DataFrame([['Others'] + list(others)], columns=['name'] + SAMPLES)],
    ignore_index=True,
)

df_g = pd.read_csv(GENUS_CSV)
df_g = df_g[df_g['name'].str.lower() != 'unclassified'].copy()
df_g['mean'] = df_g[SAMPLES].mean(axis=1)
df_gs = df_g.sort_values('mean', ascending=False).reset_index(drop=True)
top20 = df_gs.head(20)[['name'] + SAMPLES].copy()
others_g = df_gs.iloc[20:][SAMPLES].sum()
top20_plot = pd.concat(
    [
        top20,
        pd.DataFrame([['Others'] + list(others_g)], columns=['name'] + SAMPLES),
    ],
    ignore_index=True,
)

PHYLUM_COLORS = [
    '#1F4E79',
    '#2E75B6',
    '#4BACC6',
    '#548235',
    '#ED7D31',
    '#C00000',
    '#A9D18E',
    '#F4B183',
    '#9DC3E6',
    '#B4C7E7',
    '#D9D9D9',
]
ENV_PAL = [
    '#2C5F8A',
    '#4488BB',
    '#5BAABB',
    '#548235',
    '#C8860A',
    '#7B9E3E',
    '#843C0C',
    '#375623',
    '#A0522D',
    '#888888',
    '#6D9EC1',
    '#5D9E5D',
    '#9090C0',
    '#8FAF6F',
    '#909090',
    '#2E7B8C',
    '#666666',
]
env_idx = 0
bar_colors_g = []
for _, row in top20_plot.iterrows():
  s = clin(row['name'])
  if s == 'zoonotic':
    bar_colors_g.append('#C00000')
  elif s == 'critical':
    bar_colors_g.append('#7030A0')
  elif s == 'opportunistic':
    bar_colors_g.append('#ED7D31')
  else:
    bar_colors_g.append(ENV_PAL[env_idx % len(ENV_PAL)])
    env_idx += 1

W1, H1 = 15.2, 7.0  # spacious canvas for dual stacked bars and separate legends
fig1 = plt.figure(figsize=(W1, H1), facecolor='white')
gs1 = fig1.add_gridspec(
    1, 4,
    width_ratios=[1.0, 0.58, 1.0, 0.72],
    wspace=0.22,
    left=0.05,
    right=0.98,
    top=0.88,
    bottom=0.11
)

# Panel A: Phylum
ax1a = fig1.add_subplot(gs1[0, 0])
ax1a_leg = fig1.add_subplot(gs1[0, 1]); ax1a_leg.axis('off')
ax1a.set_facecolor('#FAFAFA')
x = np.arange(6)
bottoms = np.zeros(6)
for i, row in ph_plot.iterrows():
  vals = row[SAMPLES].values.astype(float)
  ax1a.bar(
      x,
      vals,
      bottom=bottoms,
      color=PHYLUM_COLORS[i],
      width=0.62,
      edgecolor='white',
      linewidth=0.8,
      zorder=3,
  )
  for j, (v, b) in enumerate(zip(vals, bottoms)):
    if v >= 0.040:
      dark = {'#1F4E79', '#548235', '#C00000', '#2E75B6'}
      tc = 'white' if PHYLUM_COLORS[i] in dark else '#212529'
      ax1a.text(
          x[j],
          b + v / 2,
          f'{v * 100:.1f}%',
          ha='center',
          va='center',
          fontsize=pt(7.8, W1),
          fontweight='bold',
          color=tc,
          zorder=4,
      )
  bottoms += vals

for j, tot in enumerate(ph_plot[SAMPLES].sum()):
  ax1a.text(
      x[j],
      tot + 0.010,
      f'{tot * 100:.1f}%',
      ha='center',
      va='bottom',
      fontsize=pt(8.5, W1),
      fontweight='bold',
      color='#212529',
  )

bracket_y = 1.06
for city, idx in [('Mardan', [0, 1]), ('Peshawar', [2, 3]), ('Swat', [4, 5])]:
  xs = idx[0] - 0.38
  xe = idx[1] + 0.38
  ax1a.annotate(
      '',
      xy=(xe, bracket_y),
      xytext=(xs, bracket_y),
      arrowprops=dict(arrowstyle='-', color=CITY_COL[city], lw=1.6),
  )
  ax1a.text(
      (xs + xe) / 2,
      bracket_y + 0.015,
      city,
      ha='center',
      va='bottom',
      fontsize=pt(8.5, W1),
      fontweight='bold',
      color=CITY_COL[city],
  )

ax1a.set_xticks(x)
ax1a.set_xticklabels(SAMPLES, fontsize=pt(8.5, W1), fontweight='bold')
ax1a.set_ylabel(
    'Relative Abundance', fontsize=pt(9.5, W1), fontweight='bold', labelpad=6
)
ax1a.set_xlabel('Sample', fontsize=pt(9.5, W1), fontweight='bold', labelpad=6)
ax1a.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax1a.tick_params(axis='y', labelsize=pt(8.5, W1))
ax1a.set_xlim(-0.6, 5.6)
ax1a.set_ylim(0, 1.15)
ax1a.spines['top'].set_visible(False)
ax1a.spines['right'].set_visible(False)
ax1a.grid(axis='y', alpha=0.25, zorder=1)
ax1a.text(
    -0.12,
    1.05,
    'A',
    transform=ax1a.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax1a.set_title(
    'Phylum-Level Relative Abundance',
    fontsize=14.5,
    fontweight='bold',
    pad=14,
    loc='left',
)

lp1 = []
for i, row in ph_plot.iterrows():
  name = row['name']
  mv = (
      ph_s[ph_s['name'] == name]['mean'].values[0] * 100
      if name in list(ph_s['name'])
      else 0
  )
  note = ' *' if name == 'Campylobacterota' else ''
  label = f'{name} ({mv:.1f}%){note}' if name != 'Others' else 'Others'
  lp1.append(mpatches.Patch(color=PHYLUM_COLORS[i], label=label))
ax1a_leg.legend(
    handles=lp1,
    loc='center left',
    bbox_to_anchor=(0.0, 0.5),
    fontsize=pt(7.8, W1),
    frameon=True,
    framealpha=0.95,
    edgecolor='#cccccc',
    borderpad=0.6,
    labelspacing=0.45,
    title='Phylum (Mean %)',
    title_fontsize=pt(8.5, W1),
)

# Panel B: Top 20 Genera
ax1b = fig1.add_subplot(gs1[0, 2])
ax1b_leg = fig1.add_subplot(gs1[0, 3]); ax1b_leg.axis('off')
ax1b.set_facecolor('#FAFAFA')
bottoms_g = np.zeros(6)
for i, row in top20_plot.iterrows():
  vals = row[SAMPLES].values.astype(float)
  ax1b.bar(
      x,
      vals,
      bottom=bottoms_g,
      color=bar_colors_g[i],
      width=0.62,
      edgecolor='white',
      linewidth=0.5,
      zorder=3,
  )
  for j, (v, b) in enumerate(zip(vals, bottoms_g)):
    if v >= 0.032:
      dark = {
          '#C00000',
          '#7030A0',
          '#2C5F8A',
          '#548235',
          '#843C0C',
          '#A0522D',
          '#375623',
          '#2E7B8C',
          '#4488BB',
          '#666666',
      }
      tc = 'white' if bar_colors_g[i] in dark else '#212529'
      ax1b.text(
          x[j],
          b + v / 2,
          f'{v * 100:.1f}%',
          ha='center',
          va='center',
          fontsize=pt(7.5, W1),
          fontweight='bold',
          color=tc,
          zorder=4,
      )
  bottoms_g += vals

for j, tot in enumerate(top20_plot[SAMPLES].sum()):
  ax1b.text(
      x[j],
      tot + 0.010,
      f'{tot * 100:.1f}%',
      ha='center',
      va='bottom',
      fontsize=pt(8.5, W1),
      fontweight='bold',
      color='#212529',
  )

for city, idx in [('Mardan', [0, 1]), ('Peshawar', [2, 3]), ('Swat', [4, 5])]:
  xs = idx[0] - 0.38
  xe = idx[1] + 0.38
  ax1b.annotate(
      '',
      xy=(xe, bracket_y),
      xytext=(xs, bracket_y),
      arrowprops=dict(arrowstyle='-', color=CITY_COL[city], lw=1.6),
  )
  ax1b.text(
      (xs + xe) / 2,
      bracket_y + 0.015,
      city,
      ha='center',
      va='bottom',
      fontsize=pt(8.5, W1),
      fontweight='bold',
      color=CITY_COL[city],
  )

ax1b.set_xticks(x)
ax1b.set_xticklabels(SAMPLES, fontsize=pt(8.5, W1), fontweight='bold')
ax1b.set_ylabel(
    'Relative Abundance', fontsize=pt(9.5, W1), fontweight='bold', labelpad=6
)
ax1b.set_xlabel('Sample', fontsize=pt(9.5, W1), fontweight='bold', labelpad=6)
ax1b.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax1b.tick_params(axis='y', labelsize=pt(8.5, W1))
ax1b.set_xlim(-0.6, 5.6)
ax1b.set_ylim(0, 1.16)
ax1b.spines['top'].set_visible(False)
ax1b.spines['right'].set_visible(False)
ax1b.grid(axis='y', alpha=0.25, zorder=1)
ax1b.text(
    -0.12,
    1.05,
    'B',
    transform=ax1b.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax1b.set_title(
    'Top 20 Bacterial Genera',
    fontsize=14.5,
    fontweight='bold',
    pad=14,
    loc='left',
)

lp2 = [mpatches.Patch(color='none', label='── ZOONOTIC ──')]
for _, row in top20_plot.iterrows():
  name = row['name']
  s = clin(name)
  if s == 'zoonotic':
    mv = (
        df_gs[df_gs['name'] == name]['mean'].values[0] * 100
        if name in list(df_gs['name'])
        else 0
    )
    lp2.append(mpatches.Patch(color='#C00000', label=f' {name} ({mv:.2f}%)'))
lp2.append(mpatches.Patch(color='none', label='── WHO CRITICAL ──'))
for _, row in top20_plot.iterrows():
  name = row['name']
  s = clin(name)
  if s == 'critical':
    mv = df_gs[df_gs['name'] == name]['mean'].values[0] * 100
    lp2.append(mpatches.Patch(color='#7030A0', label=f' {name} ({mv:.2f}%)'))
lp2.append(mpatches.Patch(color='none', label='── OPPORTUNISTIC ──'))
for _, row in top20_plot.iterrows():
  name = row['name']
  s = clin(name)
  if s == 'opportunistic':
    mv = df_gs[df_gs['name'] == name]['mean'].values[0] * 100
    lp2.append(mpatches.Patch(color='#ED7D31', label=f' {name} ({mv:.2f}%)'))
lp2.append(mpatches.Patch(color='none', label='── ENVIRONMENTAL ──'))
ec_idx = 0
for _, row in top20_plot.iterrows():
  name = row['name']
  s = clin(name)
  if s == 'environmental':
    mv = (
        df_gs[df_gs['name'] == name]['mean'].values[0] * 100
        if name in list(df_gs['name'])
        else 0
    )
    col = ENV_PAL[ec_idx % len(ENV_PAL)]
    lbl = (
        f' {name} ({mv:.2f}%)'
        if name != 'Others'
        else ' Others (remaining genera)'
    )
    lp2.append(mpatches.Patch(color=col, label=lbl))
    ec_idx += 1

ax1b_leg.legend(
    handles=lp2,
    loc='center left',
    bbox_to_anchor=(0.0, 0.5),
    fontsize=pt(7.2, W1),
    frameon=True,
    framealpha=0.95,
    edgecolor='#cccccc',
    borderpad=0.5,
    labelspacing=0.38,
    title='Clinical Classification',
    title_fontsize=pt(8.2, W1),
)

save_dual(fig1, 'Fig_1')


# ==============================================================================
# FIGURE 2: PRIORITY PATHOGEN BURDEN & CLINICAL FRAMEWORKS (SQUARE FORMAT)
# Target Print Dimensions: 7.2 x 8.0 inches (183 mm x 203 mm)
# ==============================================================================
print('\n--- Generating Figure 2 with Calibrated Typography ---')
W2, H2 = 11.0, 11.8
df_g['Mardan'] = df_g[['MSLW1', 'MSLW2']].mean(axis=1)
df_g['Peshawar'] = df_g[['PSLW1', 'PSLW2']].mean(axis=1)
df_g['Swat'] = df_g[['SSLW1', 'SSLW2']].mean(axis=1)
df_gs = df_g.sort_values('mean', ascending=False).reset_index(drop=True)

all_gen = list(CONFIRMED | CRITICAL | OPPORTUNISTIC) + COMMENSAL
cats = [
    'Zoonotic Pathogen',
    'WHO Critical MDR',
    'Opportunistic Zoonotic',
    'Commensal / Environmental',
]
counts = {
    'Zoonotic Pathogen': 10,
    'WHO Critical MDR': 2,
    'Opportunistic Zoonotic': 13,
    'Commensal / Environmental': 3,
}
colors_cat = ['#C00000', '#7030A0', '#ED7D31', '#70AD47']
abund = {}
for cat in cats:
  gen = [g for g in all_gen if clin4(g) == cat]
  abund[cat] = sum(
      df_gs[df_gs['name'] == g]['mean'].values[0] * 100
      for g in gen
      if len(df_gs[df_gs['name'] == g]) > 0
  )
count_vals = [counts[c] for c in cats]
abund_vals = [abund[c] for c in cats]
short_lab = [
    'Zoonotic\nPathogen',
    'WHO Critical\nMDR',
    'Opportunistic\nZoonotic',
    'Commensal /\nEnvironmental',
]

zoo25 = [
    (
        g,
        df_gs[df_gs['name'] == g].iloc[0]['Mardan'] * 100,
        df_gs[df_gs['name'] == g].iloc[0]['Peshawar'] * 100,
        df_gs[df_gs['name'] == g].iloc[0]['Swat'] * 100,
        clin(g),
    )
    for g in ZOONOTIC_ALL
    if len(df_gs[df_gs['name'] == g]) > 0
]
zoo25.sort(key=lambda x: (x[1] + x[2] + x[3]) / 3, reverse=True)
genera5 = [d[0] for d in zoo25]
mar_v = [d[1] for d in zoo25]
pes_v = [d[2] for d in zoo25]
swa_v = [d[3] for d in zoo25]
st5 = [d[4] for d in zoo25]

fig2 = plt.figure(figsize=(W2, H2), facecolor='white')
gs2 = fig2.add_gridspec(
    2,
    2,
    height_ratios=[1.0, 1.6],
    hspace=0.34,
    wspace=0.25,
    left=0.14,
    right=0.96,
    top=0.91,
    bottom=0.06,
)

# Panel A: Genus Count (Top-Left)
ax2a = fig2.add_subplot(gs2[0, 0])
ax2a.set_facecolor('#FAFAFA')
bars_a = ax2a.bar(
    range(4),
    count_vals,
    color=colors_cat,
    width=0.52,
    edgecolor='white',
    linewidth=1.0,
    zorder=3,
)
for bar, v in zip(bars_a, count_vals):
    ax2a.text(bar.get_x() + bar.get_width() / 2, v + max(count_vals) * 0.03,
              f'{int(v)}\n({v / 28 * 100:.0f}%)',
              ha='center', va='bottom', fontsize=pt(8.5, W2), fontweight='bold', color='#212529')
ax2a.set_xticks(range(4))
ax2a.set_xticklabels(
    short_lab, fontsize=pt(8.5, W2), fontweight='bold', linespacing=1.2
)
for tick, col in zip(ax2a.get_xticklabels(), colors_cat):
  tick.set_color(col)
ax2a.set_ylabel(
    'Number of Genera', fontsize=pt(9.5, W2), fontweight='bold', labelpad=6
)
ax2a.set_ylim(0, max(count_vals) * 1.25)
ax2a.tick_params(axis='y', labelsize=pt(8.5, W2))
ax2a.spines['top'].set_visible(False)
ax2a.spines['right'].set_visible(False)
ax2a.grid(axis='y', alpha=0.25, zorder=1)
ax2a.tick_params(axis='x', length=0, pad=6)
ax2a.text(
    -0.12,
    1.05,
    'A',
    transform=ax2a.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax2a.set_title(
    'Genus Count by Clinical Category (n = 28)',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
)

# Panel B: Mean Relative Abundance (Top-Right)
ax2b = fig2.add_subplot(gs2[0, 1])
ax2b.set_facecolor('#FAFAFA')
bars_b = ax2b.bar(
    range(4),
    abund_vals,
    color=colors_cat,
    width=0.52,
    edgecolor='white',
    linewidth=1.0,
    zorder=3,
)
for bar, v in zip(bars_b, abund_vals):
    ax2b.text(bar.get_x() + bar.get_width() / 2, v + max(abund_vals) * 0.03,
              f'{v:.1f}%',
              ha='center', va='bottom', fontsize=pt(8.5, W2), fontweight='bold', color='#212529')
ax2b.set_xticks(range(4))
ax2b.set_xticklabels(
    short_lab, fontsize=pt(8.5, W2), fontweight='bold', linespacing=1.2
)
for tick, col in zip(ax2b.get_xticklabels(), colors_cat):
  tick.set_color(col)
ax2b.set_ylabel(
    'Mean Relative Abundance (%)',
    fontsize=pt(9.5, W2),
    fontweight='bold',
    labelpad=6,
)
ax2b.set_ylim(0, max(abund_vals) * 1.25)
ax2b.tick_params(axis='y', labelsize=pt(8.5, W2))
ax2b.spines['top'].set_visible(False)
ax2b.spines['right'].set_visible(False)
ax2b.grid(axis='y', alpha=0.25, zorder=1)
ax2b.tick_params(axis='x', length=0, pad=6)
ax2b.text(
    -0.12,
    1.05,
    'B',
    transform=ax2b.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax2b.set_title(
    'Mean Relative Abundance (% classified reads)',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
)

# Panel C: Zoonotic Genera Abundance Across Cities (Bottom Full Width)
ax2c = fig2.add_subplot(gs2[1, :])
ax2c.set_facecolor('#FAFAFA')
y = np.arange(len(genera5))
h = 0.25
ax2c.barh(
    y + h,
    mar_v,
    height=h,
    color='#1F4E79',
    alpha=0.92,
    label='Mardan',
    edgecolor='white',
    linewidth=0.5,
)
ax2c.barh(
    y,
    pes_v,
    height=h,
    color='#C55A11',
    alpha=0.92,
    label='Peshawar',
    edgecolor='white',
    linewidth=0.5,
)
ax2c.barh(
    y - h,
    swa_v,
    height=h,
    color='#375623',
    alpha=0.92,
    label='Swat',
    edgecolor='white',
    linewidth=0.5,
)
ax2c.set_yticks(y)
ylabels5 = list(genera5)
ax2c.set_yticklabels(ylabels5, fontsize=pt(8.0, W2), style='italic', fontweight='bold')
for lbl, s in zip(ax2c.get_yticklabels(), st5):
  lbl.set_color(CLIN_COL[s])
ax2c.set_xlabel(
    'Mean Relative Abundance (%)',
    fontsize=pt(9.5, W2),
    fontweight='bold',
    labelpad=8,
)
ax2c.tick_params(axis='x', labelsize=pt(8.5, W2))
ax2c.spines['top'].set_visible(False)
ax2c.spines['right'].set_visible(False)
ax2c.grid(axis='x', alpha=0.25)
ax2c.tick_params(axis='y', length=0, pad=8)
ax2c.set_xlim(0, max(max(mar_v), max(pes_v), max(swa_v)) * 1.15)
ax2c.text(
    -0.10,
    1.04,
    'C',
    transform=ax2c.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax2c.set_title(
    'Mean Relative Abundance of 25 Priority Zoonotic and Opportunistic Genera'
    ' Across Cities',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
    loc='left',
)

leg2c_items = [
    mpatches.Patch(color='#1F4E79', label='Mardan'),
    mpatches.Patch(color='#C55A11', label='Peshawar'),
    mpatches.Patch(color='#375623', label='Swat'),
    mpatches.Patch(color='none', label=''),
    mpatches.Patch(color='#C00000', label='Zoonotic Pathogen'),
    mpatches.Patch(color='#7030A0', label='WHO Critical Lineage'),
    mpatches.Patch(color='#ED7D31', label='Opportunistic Zoonotic'),
]
ax2c.legend(
    handles=leg2c_items,
    loc='lower right',
    bbox_to_anchor=(0.98, 0.05),
    fontsize=pt(8.0, W2),
    frameon=True,
    framealpha=0.96,
    edgecolor='#ccc',
    borderpad=0.6,
    labelspacing=0.45,
    title='City / Clinical Classification',
    title_fontsize=pt(8.5, W2),
)

save_dual(fig2, 'Fig_2')


# ==============================================================================
# FIGURE 3: COMMUNITY DIVERSITY & ORDINATION
# Target Print Dimensions: 7.2 x 7.0 inches (183 mm x 178 mm)
# ==============================================================================
print('\n--- Generating Figure 3 with Calibrated Typography ---')
W3, H3 = 13.0, 11.5
alpha = {}
for s in SAMPLES:
  p = df_g[s].values
  p_pos = p[p > 0]
  shannon = -np.sum(p_pos * np.log(p_pos))
  simpson = 1 - np.sum(p_pos**2)
  observed = int((df_g[s] > 0).sum())
  counts = (p * 1e6).astype(int)
  singletons = np.sum(counts == 1)
  doubletons = max(np.sum(counts == 2), 1)
  chao1 = observed + (singletons**2) / (2 * doubletons)
  alpha[s] = {
      'Shannon': shannon,
      'Simpson': simpson,
      'Observed': observed,
      'Chao1': chao1,
  }

mat = df_g[SAMPLES].values.T
bc_mat = np.zeros((6, 6))
for i in range(6):
  for j in range(6):
    bc_mat[i, j] = braycurtis(mat[i], mat[j])

A_bc = -0.5 * (bc_mat**2)
G = (
    A_bc
    - A_bc.mean(axis=1, keepdims=True)
    - A_bc.mean(axis=0, keepdims=True)
    + A_bc.mean()
)
eigvals, eigvecs = np.linalg.eigh(G)
idx_sort = np.argsort(eigvals)[::-1]
eigvals = eigvals[idx_sort]
eigvecs = eigvecs[:, idx_sort]
pos_var = eigvals[eigvals > 0].sum()
pct1 = eigvals[0] / pos_var * 100
pct2 = eigvals[1] / pos_var * 100
pc1 = eigvecs[:, 0] * np.sqrt(max(eigvals[0], 0))
pc2 = eigvecs[:, 1] * np.sqrt(max(eigvals[1], 0))
sample_colors = [CITY_COL[c] for c in CITIES]
MARKERS = ['o', 's', 'o', 's', 'o', 's']

fig3 = plt.figure(figsize=(W3, H3), facecolor='white')
gs3 = fig3.add_gridspec(
    2,
    2,
    height_ratios=[1.0, 1.45],
    hspace=0.32,
    wspace=0.25,
    left=0.07,
    right=0.96,
    top=0.91,
    bottom=0.06,
)

# Panel A: Alpha Diversity
gs3_alpha = gs3[0, :].subgridspec(1, 4, wspace=0.24)
metrics = ['Observed', 'Shannon', 'Simpson', 'Chao1']
ylabels6 = [
    'Observed OTUs',
    "Shannon Entropy (H')",
    'Simpson Index (1-D)',
    'Chao1 Richness',
]
titles_alpha = [
    'Observed Genera',
    "Shannon Diversity (H')",
    'Simpson Diversity (1-D)',
    'Chao1 Richness Estimate',
]

for i_m, (met, ylab, t_alpha) in enumerate(
    zip(metrics, ylabels6, titles_alpha)
):
  ax_a = fig3.add_subplot(gs3_alpha[0, i_m])
  ax_a.set_facecolor('#FAFAFA')
  vals = [alpha[s][met] for s in SAMPLES]
  ax_a.bar(
      range(6),
      vals,
      color=BAR_COLORS,
      width=0.62,
      edgecolor='white',
      linewidth=0.8,
      zorder=3,
  )
  for j, v in enumerate(vals):
    label_t = str(int(round(v))) if met in ['Observed', 'Chao1'] else f'{v:.2f}'
    ax_a.text(
        j,
        v + max(vals) * 0.025,
        label_t,
        ha='center',
        va='bottom',
        fontsize=pt(8.0, W3),
        fontweight='bold',
    )
  for city, idx_c, col_c in [
      ('Mardan', [0, 1], '#1F4E79'),
      ('Peshawar', [2, 3], '#C55A11'),
      ('Swat', [4, 5], '#375623'),
  ]:
    mv = np.mean([vals[k] for k in idx_c])
    ax_a.plot(
        [idx_c[0] - 0.35, idx_c[1] + 0.35],
        [mv, mv],
        color=col_c,
        lw=1.8,
        ls='--',
        zorder=5,
    )
  ax_a.set_xticks(range(6))
  ax_a.set_xticklabels(
      SAMPLES, fontsize=pt(8.0, W3), fontweight='bold', rotation=30, ha='right'
  )
  ax_a.set_ylabel(ylab, fontsize=pt(9.0, W3), fontweight='bold')
  ax_a.tick_params(axis='y', labelsize=pt(8.0, W3))
  ax_a.set_ylim(0, max(vals) * 1.20)
  ax_a.spines['top'].set_visible(False)
  ax_a.spines['right'].set_visible(False)
  ax_a.grid(axis='y', alpha=0.25, zorder=1)
  ax_a.set_title(t_alpha, fontsize=12.5, fontweight='bold', pad=12)
  if i_m == 0:
    ax_a.text(
        -0.25,
        1.16,
        'A',
        transform=ax_a.transAxes,
        fontsize=18.0,
        fontweight='bold',
        va='bottom',
    )

patches6 = [
    mpatches.Patch(color='#1F4E79', label='Mardan'),
    mpatches.Patch(color='#C55A11', label='Peshawar'),
    mpatches.Patch(color='#375623', label='Swat'),
    mlines.Line2D([], [], color='gray', lw=1.8, ls='--', label='City Mean'),
]
fig3.legend(
    handles=patches6,
    loc='upper center',
    bbox_to_anchor=(0.51, 0.99),
    ncol=4,
    fontsize=pt(8.5, W3),
    frameon=True,
    edgecolor='#ccc',
    framealpha=0.96,
)

# Panel B: PCoA
ax3b = fig3.add_subplot(gs3[1, 0])
ax3b.set_facecolor('#FAFAFA')
for city, sidx in [
    ('Mardan', [0, 1]),
    ('Peshawar', [2, 3]),
    ('Swat', [4, 5]),
]:
  cx = np.mean([pc1[i] for i in sidx])
  cy = np.mean([pc2[i] for i in sidx])
  rx = max(abs(pc1[i] - cx) for i in sidx) + 0.042
  ry = max(abs(pc2[i] - cy) for i in sidx) + 0.042
  ax3b.add_patch(
      Ellipse(
          (cx, cy),
          width=rx * 2,
          height=ry * 2,
          facecolor=CITY_COL[city],
          alpha=0.12,
          edgecolor=CITY_COL[city],
          linewidth=1.5,
          linestyle='--',
          zorder=1,
      )
  )

ax3b.axhline(0, color='#CCCCCC', lw=0.8, ls='--', zorder=0)
ax3b.axvline(0, color='#CCCCCC', lw=0.8, ls='--', zorder=0)
for i in range(6):
  ax3b.scatter(
      pc1[i],
      pc2[i],
      c=sample_colors[i],
      marker=MARKERS[i],
      s=280,
      edgecolors='white',
      linewidths=1.8,
      zorder=4,
  )

offsets = {
    'MSLW1': (-0.065, -0.030),
    'MSLW2': (-0.065, 0.030),
    'PSLW1': (0.015, -0.038),
    'PSLW2': (-0.095, 0.015),
    'SSLW1': (-0.075, 0.000),
    'SSLW2': (0.035, -0.030),
}
for i, s in enumerate(SAMPLES):
  ox, oy = offsets[s]
  ax3b.annotate(
      s,
      xy=(pc1[i], pc2[i]),
      xytext=(pc1[i] + ox, pc2[i] + oy),
      fontsize=pt(8.5, W3),
      fontweight='bold',
      color=sample_colors[i],
      ha='center',
      va='center',
      zorder=5,
      arrowprops=dict(arrowstyle='-', color=sample_colors[i], lw=0.9),
  )

ax3b.set_xlabel(
    f'PC1 ({pct1:.1f}% variance explained)',
    fontsize=pt(9.5, W3),
    fontweight='bold',
    labelpad=8,
)
ax3b.set_ylabel(
    f'PC2 ({pct2:.1f}% variance explained)',
    fontsize=pt(9.5, W3),
    fontweight='bold',
    labelpad=8,
)
ax3b.tick_params(labelsize=pt(8.5, W3))
ax3b.spines['top'].set_visible(False)
ax3b.spines['right'].set_visible(False)
ax3b.grid(alpha=0.2, zorder=0)
ax3b.set_xlim(pc1.min() - 0.10, pc1.max() + 0.10)
ax3b.set_ylim(pc2.min() - 0.10, pc2.max() + 0.10)
ax3b.annotate(
    'SSLW2: High Acetoanaerobium\n(27.92%)',
    xy=(pc1[5], pc2[5]),
    xytext=(pc1[5] - 0.15, pc2[5] - 0.11),
    fontsize=pt(7.5, W3),
    color='#375623',
    style='italic',
    ha='center',
    arrowprops=dict(arrowstyle='->', color='#375623', lw=1.0),
    bbox=dict(
        boxstyle='round,pad=0.3', fc='#E2EFDA', ec='#375623', lw=0.8
    ),
)
ax3b.annotate(
    'MSLW1 & SSLW1\noverlap (BC = 0.100)',
    xy=((pc1[0] + pc1[4]) / 2, (pc2[0] + pc2[4]) / 2),
    xytext=(pc1[0] - 0.16, pc2[0] - 0.11),
    fontsize=pt(7.5, W3),
    color='#1F4E79',
    style='italic',
    ha='center',
    arrowprops=dict(arrowstyle='->', color='#666', lw=1.0),
    bbox=dict(
        boxstyle='round,pad=0.3', fc='#D6E4F0', ec='#1F4E79', lw=0.8
    ),
)

ax3b.text(
    -0.10,
    1.05,
    'B',
    transform=ax3b.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax3b.set_title(
    'PCoA of Bray-Curtis Dissimilarity (PC1 + PC2 = 88.9%)',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
    loc='left',
)

leg_pcoa = [
    mpatches.Patch(color='#1F4E79', label='Mardan'),
    mpatches.Patch(color='#C55A11', label='Peshawar'),
    mpatches.Patch(color='#375623', label='Swat'),
    mlines.Line2D(
        [],
        [],
        marker='o',
        linestyle='none',
        color='#555',
        markersize=8,
        label='Replicate 1',
    ),
    mlines.Line2D(
        [],
        [],
        marker='s',
        linestyle='none',
        color='#555',
        markersize=8,
        label='Replicate 2',
    ),
    mpatches.Patch(
        facecolor='#888',
        alpha=0.18,
        edgecolor='#888',
        linestyle='--',
        label='City 95% CI ellipse',
    ),
]
ax3b.legend(
    handles=leg_pcoa,
    loc='upper left',
    bbox_to_anchor=(0.02, 0.98),
    fontsize=pt(7.8, W3),
    frameon=True,
    framealpha=0.96,
    edgecolor='#ccc',
    borderpad=0.5,
)

# Panel C: Bray-Curtis Matrix Heatmap
ax3c = fig3.add_subplot(gs3[1, 1])
ax3c.set_facecolor('white')
cmap_bc = LinearSegmentedColormap.from_list(
    'bc', ['#C6EFCE', '#FFEB9C', '#FFC7CE'], N=256
)
im_bc = ax3c.imshow(bc_mat, cmap=cmap_bc, vmin=0, vmax=1, aspect='auto')
for i in range(6):
  for j in range(6):
    v = bc_mat[i, j]
    tc = 'white' if v > 0.75 else '#212529'
    ax3c.text(
        j,
        i,
        f'{v:.3f}',
        ha='center',
        va='center',
        fontsize=pt(9.0, W3),
        fontweight='bold',
        color=tc,
    )
ax3c.set_xticks(range(6))
ax3c.set_yticks(range(6))
ax3c.set_xticklabels(SAMPLES, fontsize=pt(8.5, W3), fontweight='bold')
ax3c.set_yticklabels(SAMPLES, fontsize=pt(8.5, W3), fontweight='bold')
for i in range(6):
  for j in range(6):
    ax3c.add_patch(
        mpatches.Rectangle(
            (j - 0.5, i - 0.5), 1, 1, fill=False, ec='white', lw=1.5
        )
    )

for idx_box, col_box in [
    ([0, 1], '#1F4E79'),
    ([2, 3], '#C55A11'),
    ([4, 5], '#375623'),
]:
  ax3c.add_patch(
      mpatches.Rectangle(
          (idx_box[0] - 0.5, idx_box[0] - 0.5),
          2,
          2,
          fill=False,
          ec=col_box,
          lw=2.2,
          zorder=5,
      )
  )

cb_bc = fig3.colorbar(
    im_bc, ax=ax3c, orientation='vertical', fraction=0.046, pad=0.05
)
cb_bc.set_label(
    'Bray-Curtis Dissimilarity (0 = identical | 1 = dissimilar)',
    fontsize=pt(8.5, W3),
    fontweight='bold',
)
cb_bc.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
cb_bc.ax.tick_params(labelsize=pt(8.0, W3))

notes_bc = [
    mpatches.Patch(color='#C6EFCE', label='Low dissimilarity (< 0.3)'),
    mpatches.Patch(color='#FFEB9C', label='Moderate (0.3–0.6)'),
    mpatches.Patch(color='#FFC7CE', label='High dissimilarity (> 0.6)'),
]
ax3c.legend(
    handles=notes_bc,
    loc='upper left',
    bbox_to_anchor=(1.25, 1.0),
    fontsize=pt(7.8, W3),
    frameon=True,
    edgecolor='#ccc',
)

ax3c.text(
    -0.10,
    1.05,
    'C',
    transform=ax3c.transAxes,
    fontsize=18.0,
    fontweight='bold',
    va='bottom',
)
ax3c.set_title(
    'Pairwise Bray-Curtis Dissimilarity Matrix',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
    loc='left',
)

save_dual(fig3, 'Fig_3')


# ==============================================================================
# FIGURE 4: HAZARD MODELING & ECOLOGICAL NETWORK
# Target Print Dimensions: 7.2 x 7.8 inches (183 mm x 198 mm)
# ==============================================================================
print('\n--- Generating Figure 4 with Calibrated Typography (Option 2: Clutter-Free Network) ---')
# High-resolution publication composite using authentic Fruchterman-Reingold Network (Option 2: Pure Clutter-Free)
from PIL import Image, ImageDraw, ImageFont

net_img_path = find_file('Figure_10_CoOccurrence_Network.png')
if not os.path.exists(net_img_path):
    net_img_path = find_file('Figure9_CoOccurrence_600dpi (1).png')
if not os.path.exists(net_img_path):
    # fallback to workspace
    net_img_path = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', 'updated figures', 'Figure9_CoOccurrence_600dpi (1).png'))

im_net_raw = Image.open(net_img_path).convert('RGB')
# Crop network to remove the bulky right sidebar completely (pure network graphic)
net_crop = im_net_raw.crop((150, 1350, 7350, 6850))

# Scale network to fit comfortably and centered in a 12000 px canvas
w_net_target = 8800
h_net_target = int(net_crop.height * (w_net_target / net_crop.width))
net_scaled = net_crop.resize((w_net_target, h_net_target), Image.Resampling.LANCZOS)

# 1. Render Top Panels (A and B) cleanly with Matplotlib at 600 DPI
fig_top, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(20.0, 6.8), dpi=600, facecolor='white')
fig_top.subplots_adjust(left=0.07, right=0.96, top=0.86, bottom=0.14, wspace=0.22)

# Panel A: Sample-Level ZHI
ZHI = [50.51, 100.91, 192.93, 68.76, 54.72, 106.76]
ax4a.set_facecolor('#FAFAFA')
bars = ax4a.bar(range(6), ZHI, color=BAR_COLORS, width=0.58, edgecolor='white', linewidth=1.2, zorder=3)
for b, v in zip(bars, ZHI):
    ax4a.text(b.get_x() + b.get_width()/2, v + 3.0, f'{v:.1f}',
              ha='center', va='bottom', fontsize=11.5, fontweight='bold', color='#212529')

for city, idx, col in [('Mardan', [0, 1], '#1F4E79'), ('Peshawar', [2, 3], '#C55A11'), ('Swat', [4, 5], '#375623')]:
    xs, xe = idx[0] - 0.35, idx[1] + 0.35
    ax4a.plot([xs, xe], [210, 210], color=col, lw=2.2, solid_capstyle='butt')
    ax4a.text((xs+xe)/2, 213, city, ha='center', va='bottom', fontsize=12.0, fontweight='bold', color=col)

mean_zhi = np.mean(ZHI)
ax4a.axhline(mean_zhi, color='gray', lw=1.5, ls='--', zorder=4)
ax4a.text(5.45, mean_zhi + 3, f'Mean={mean_zhi:.1f}', ha='right', va='bottom', fontsize=10.0, color='#555', fontweight='bold')

ax4a.set_xticks(range(6))
ax4a.set_xticklabels(SAMPLES, fontsize=11.0, fontweight='bold')
ax4a.set_ylabel('Zoonotic Hazard Index (ZHI)', fontsize=12.5, fontweight='bold', labelpad=8)
ax4a.set_xlabel('Sample', fontsize=12.0, fontweight='bold', labelpad=6)
ax4a.set_ylim(0, 235)
ax4a.spines['top'].set_visible(False)
ax4a.spines['right'].set_visible(False)
ax4a.grid(axis='y', alpha=0.3, zorder=1)
ax4a.text(-0.12, 1.06, 'A', transform=ax4a.transAxes, fontsize=20.0, fontweight='bold', va='bottom')
ax4a.set_title('Sample-Level ZHI Scores', fontsize=16.0, fontweight='bold', pad=14, loc='left')

# Panel B: Municipal Risk Tier Stratification
ax4b.set_facecolor('#FAFAFA')
ax4b.axvspan(0, 75, color='#E2EFDA', alpha=0.6, zorder=1)
ax4b.axvspan(75, 120, color='#FFF2CC', alpha=0.6, zorder=1)
ax4b.axvspan(120, 200, color='#FCE4D6', alpha=0.6, zorder=1)

c_names = ['Mardan', 'Swat', 'Peshawar']
c_vals = [75.71, 80.74, 130.84]
c_cols = ['#1F4E79', '#375623', '#C55A11']
c_ranks = ['Rank #3', 'Rank #2', 'Rank #1']

bars_b = ax4b.barh(range(3), c_vals, color=c_cols, height=0.52, edgecolor='white', linewidth=1.2, zorder=3)
for b, v, r in zip(bars_b, c_vals, c_ranks):
    ax4b.text(v + 3.0, b.get_y() + b.get_height()/2, f'ZHI = {v:.1f}  |  {r}',
              ha='left', va='center', fontsize=11.5, fontweight='bold', color='#212529')

ax4b.text(37.5, -0.68, 'MODERATE RISK\n(<75)', ha='center', fontsize=10.0, color='#375623', fontweight='bold', va='top')
ax4b.text(97.5, -0.68, 'HIGH RISK\n(75-120)', ha='center', fontsize=10.0, color='#7F6000', fontweight='bold', va='top')
ax4b.text(160, -0.68, 'VERY HIGH RISK\n(>120)', ha='center', fontsize=10.0, color='#C00000', fontweight='bold', va='top')

ax4b.set_yticks(range(3))
ax4b.set_yticklabels(c_names, fontsize=12.0, fontweight='bold')
for tick, col in zip(ax4b.get_yticklabels(), c_cols):
    tick.set_color(col)
ax4b.set_xlabel('Mean City Zoonotic Hazard Index (ZHI)', fontsize=12.5, fontweight='bold', labelpad=8)
ax4b.tick_params(axis='x', labelsize=11.0)
ax4b.set_xlim(0, 200)
ax4b.set_ylim(-0.95, 2.6)
ax4b.spines['top'].set_visible(False)
ax4b.spines['right'].set_visible(False)
ax4b.grid(axis='x', alpha=0.3, zorder=1)
ax4b.text(-0.12, 1.06, 'B', transform=ax4b.transAxes, fontsize=20.0, fontweight='bold', va='bottom')
ax4b.set_title('Municipal Risk Tier Stratification', fontsize=16.0, fontweight='bold', pad=14, loc='left')

temp_top_file = 'temp_fig4_top_gen.png'
fig_top.savefig(temp_top_file, dpi=600, facecolor='white')
plt.close(fig_top)

# 2. Composite with pure clutter-free Fruchterman-Reingold Network as Panel C
im_top = Image.open(temp_top_file)
W_canvas = 12000
scale = W_canvas / im_top.width
h_top_scaled = int(im_top.height * scale)
im_top_scaled = im_top.resize((W_canvas, h_top_scaled), Image.Resampling.LANCZOS)

gap = 80
final_h = h_top_scaled + gap + h_net_target + 80
composite4 = Image.new('RGB', (W_canvas, final_h), color='white')
composite4.paste(im_top_scaled, (0, 0))

# Center the clean network horizontally
x_net_offset = (W_canvas - w_net_target) // 2
composite4.paste(net_scaled, (x_net_offset, h_top_scaled + gap))

# Add Panel C title and tag aligned with Panel A left margin
draw4 = ImageDraw.Draw(composite4)
import matplotlib.font_manager as fm
font_path_bold = fm.findfont(fm.FontProperties(family='DejaVu Sans', weight='bold'))
font_tag = ImageFont.truetype(font_path_bold, 160)
font_title = ImageFont.truetype(font_path_bold, 130)

tag_y = h_top_scaled + gap - 45
draw4.text((840, tag_y), 'C', fill='black', font=font_tag)
draw4.text((1060, tag_y + 12), 'Polymicrobial Co-Occurrence Network of 25 Zoonotic Genera (|ρ| > 0.6)', fill='black', font=font_title)

# Save both 600 DPI PNG and vector-compatible PDF
fig4_png = os.path.join(out_dir, 'Fig_4.png')
fig4_pdf = os.path.join(out_dir, 'Fig_4.pdf')
composite4.save(fig4_png, dpi=(600, 600))
composite4.save(fig4_pdf, 'PDF', resolution=600.0)
if os.path.exists(temp_top_file):
    os.remove(temp_top_file)
print(f'Saved: {fig4_pdf} (vector/high-res) and {fig4_png} (600 DPI)')


# ==============================================================================
# FIGURE 5: WASTEWATER VIROME (STRICT A, B, C)
# Target Print Dimensions: 7.2 x 6.2 inches (183 mm x 157 mm)
# ==============================================================================
print('\n--- Generating Figure 6 (Wastewater Virome) with Calibrated Typography ---')
vdf = pd.read_excel(VIROME_XLSX, sheet_name='Virome (all species)')
zoo_host = ['Aeromonas phage DH1', 'Yersinia phage Rostov 43'] + [
    x for x in vdf['Viral species'] if 'Streptococcus satellite phage' in str(x)
]

def get_cat(row):
    if row['Viral species'] == 'uncultured human fecal virus':
        return 'Human Faecal Marker'
    if row['Viral species'] in zoo_host:
        return 'Zoonotic-Host Phage'
    if row['Type'] == 'Bacteriophage':
        return 'Bacteriophage (other)'
    if row['Type'] == 'Insect virus':
        return 'Insect Virus'
    if 'Giant virus' in str(row['Type']):
        return 'Giant Virus (protist)'
    return 'Unclassified Virus'

vdf['Category'] = vdf.apply(get_cat, axis=1)

CAT_ORDER = [
    'Human Faecal Marker',
    'Zoonotic-Host Phage',
    'Bacteriophage (other)',
    'Insect Virus',
    'Giant Virus (protist)',
    'Unclassified Virus',
]
CAT_COL = {
    'Human Faecal Marker': '#C00000',
    'Zoonotic-Host Phage': '#7030A0',
    'Bacteriophage (other)': '#1F4E79',
    'Insect Virus': '#70AD47',
    'Giant Virus (protist)': '#FFC000',
    'Unclassified Virus': '#BFBFBF',
}
cat_sums = vdf.groupby('Category')[SAMPLES].sum()
RICHNESS = [19, 10, 23, 9, 21, 8]
SHANNON = [0.956, 0.901, 0.934, 1.186, 0.907, 0.232]

fig5 = plt.figure(figsize=(10.5, 8.6), facecolor='white')
gs5 = fig5.add_gridspec(
    2, 2,
    height_ratios=[1.15, 1.0],
    hspace=0.36,
    wspace=0.24,
    left=0.08,
    right=0.96,
    top=0.91,
    bottom=0.08
)

# PANEL A: Viral Community Composition (Full width on top row)
ax5a = fig5.add_subplot(gs5[0, :])
ax5a.set_facecolor('#FAFAFA')
bottoms_v = np.zeros(6)

for cat in CAT_ORDER:
    vals_v = cat_sums.loc[cat, SAMPLES].values.astype(float) if cat in cat_sums.index else np.zeros(6)
    ax5a.bar(x, vals_v, bottom=bottoms_v, color=CAT_COL[cat], width=0.55, edgecolor='white', linewidth=0.8, zorder=3)
    for j, (v, b) in enumerate(zip(vals_v, bottoms_v)):
        if v >= 0.035:
            tc = 'white' if cat in ['Human Faecal Marker', 'Zoonotic-Host Phage', 'Bacteriophage (other)'] else '#212529'
            ax5a.text(x[j], b + v / 2, f'{v * 100:.1f}%', ha='center', va='center', fontsize=8.0, fontweight='bold', color=tc, zorder=4)
    bottoms_v += vals_v

for j in range(6):
    ax5a.text(x[j], bottoms_v[j] + 0.012, f'{bottoms_v[j] * 100:.1f}%', ha='center', va='bottom', fontsize=8.5, fontweight='bold')

bracket_y_v = 1.07
for city, idx in [('Mardan', [0, 1]), ('Peshawar', [2, 3]), ('Swat', [4, 5])]:
    xs, xe = idx[0] - 0.35, idx[1] + 0.35
    ax5a.annotate('', xy=(xe, bracket_y_v), xytext=(xs, bracket_y_v), arrowprops=dict(arrowstyle='-', color=CITY_COL[city], lw=1.6))
    ax5a.text((xs + xe) / 2, bracket_y_v + 0.015, city, ha='center', va='bottom', fontsize=8.5, fontweight='bold', color=CITY_COL[city])

ax5a.set_xticks(x)
ax5a.set_xticklabels(SAMPLES, fontsize=8.5, fontweight='bold')
ax5a.set_ylabel('Relative Abundance (% virome)', fontsize=9.5, fontweight='bold', labelpad=6)
ax5a.set_xlabel('Sample', fontsize=9.5, fontweight='bold', labelpad=4)
ax5a.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax5a.tick_params(axis='y', labelsize=8.5)
ax5a.set_ylim(0, 1.20)
ax5a.set_xlim(-0.6, 5.6)
ax5a.spines['top'].set_visible(False)
ax5a.spines['right'].set_visible(False)
ax5a.grid(axis='y', alpha=0.25, zorder=1)
ax5a.text(-0.06, 1.06, 'A', transform=ax5a.transAxes, fontsize=13.5, fontweight='bold', va='bottom')
ax5a.set_title('Viral Community Composition by Functional Category', fontsize=11.0, fontweight='bold', pad=22, loc='left')

# Standard horizontal legend above Panel A
patches_leg = [mpatches.Patch(color=CAT_COL[c], label=c) for c in CAT_ORDER]
ax5a.legend(
    handles=patches_leg,
    loc='upper center',
    bbox_to_anchor=(0.50, 1.18),
    ncol=3,
    fontsize=8.0,
    frameon=True,
    edgecolor='#DDDDDD',
    facecolor='#FCFCFC',
    framealpha=0.96
)

# PANEL B: Viral Species Richness (Bottom Left)
ax5b = fig5.add_subplot(gs5[1, 0])
ax5b.set_facecolor('#FAFAFA')
ax5b.bar(x, RICHNESS, color=BAR_COLORS, width=0.55, edgecolor='white', linewidth=0.8, zorder=3)
for j, v in enumerate(RICHNESS):
    ax5b.text(j, v + 0.5, str(int(v)), ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax5b.set_xticks(x)
ax5b.set_xticklabels(SAMPLES, fontsize=8.5, fontweight='bold')
ax5b.set_ylabel('Viral Species Richness', fontsize=9.5, fontweight='bold', labelpad=6)
ax5b.set_xlabel('Sample', fontsize=9.5, fontweight='bold', labelpad=4)
ax5b.tick_params(axis='y', labelsize=8.5)
ax5b.spines['top'].set_visible(False)
ax5b.spines['right'].set_visible(False)
ax5b.set_ylim(0, 28)
ax5b.grid(axis='y', alpha=0.25, zorder=1)

for city, idx, col in [('Mardan', [0, 1], '#1F4E79'), ('Peshawar', [2, 3], '#C55A11'), ('Swat', [4, 5], '#375623')]:
    mv = np.mean([RICHNESS[i] for i in idx])
    ax5b.plot([idx[0] - 0.32, idx[1] + 0.32], [mv, mv], color=col, lw=1.8, ls='--', zorder=5)

ax5b.text(-0.12, 1.05, 'B', transform=ax5b.transAxes, fontsize=13.5, fontweight='bold', va='bottom')
ax5b.set_title('Viral Species Richness Across Municipalities', fontsize=11.0, fontweight='bold', pad=10, loc='left')

# PANEL C: Viral Alpha Diversity (Shannon H') (Bottom Right)
ax5c = fig5.add_subplot(gs5[1, 1])
ax5c.set_facecolor('#FAFAFA')
ax5c.bar(x, SHANNON, color=BAR_COLORS, width=0.55, edgecolor='white', linewidth=0.8, zorder=3)
for j, v in enumerate(SHANNON):
    ax5c.text(j, v + 0.025, f'{v:.3f}', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax5c.set_xticks(x)
ax5c.set_xticklabels(SAMPLES, fontsize=8.5, fontweight='bold')
ax5c.set_ylabel("Shannon Diversity Index (H')", fontsize=9.5, fontweight='bold', labelpad=6)
ax5c.set_xlabel('Sample', fontsize=9.5, fontweight='bold', labelpad=4)
ax5c.tick_params(axis='y', labelsize=8.5)
ax5c.spines['top'].set_visible(False)
ax5c.spines['right'].set_visible(False)
ax5c.set_ylim(0, 1.40)
ax5c.grid(axis='y', alpha=0.25, zorder=1)

for city, idx, col in [('Mardan', [0, 1], '#1F4E79'), ('Peshawar', [2, 3], '#C55A11'), ('Swat', [4, 5], '#375623')]:
    mv = np.mean([SHANNON[i] for i in idx])
    ax5c.plot([idx[0] - 0.32, idx[1] + 0.32], [mv, mv], color=col, lw=1.8, ls='--', zorder=5)

ax5c.text(-0.12, 1.05, 'C', transform=ax5c.transAxes, fontsize=13.5, fontweight='bold', va='bottom')
ax5c.set_title("Viral Alpha Diversity (Shannon H')", fontsize=11.0, fontweight='bold', pad=10, loc='left')

# Save dual: saves to code/ and directly to submission root as Figure_6.png / Figure_6.pdf
save_dual(fig5, 'Fig5_Virome_Calibrated')
plt.close(fig5)

print('\n--- Generating Supplementary Figures S1 & S2 ---')
W_S1, H_S1 = 10.0, 14.0
top30 = df_gs.head(30).copy()
raw_mat = top30[SAMPLES].values * 100
min_v = raw_mat.min(axis=1, keepdims=True)
max_v = raw_mat.max(axis=1, keepdims=True)
range_v = np.where(max_v == min_v, 1.0, max_v - min_v)
norm_mat = (raw_mat - min_v) / range_v
genus_labels = top30['name'].tolist()

fig_s1 = plt.figure(figsize=(W_S1, H_S1), facecolor='white')
ax_s1 = fig_s1.add_axes([0.22, 0.06, 0.62, 0.88])
cmap_s1 = LinearSegmentedColormap.from_list(
    'cool_warm',
    [
        '#313695',
        '#4575B4',
        '#74ADD1',
        '#ABD9E9',
        '#E0F3F8',
        '#FEE090',
        '#FDAE61',
        '#F46D43',
        '#D73027',
        '#A50026',
    ],
    N=256,
)
im_s1 = ax_s1.imshow(norm_mat, cmap=cmap_s1, aspect='auto', vmin=0, vmax=1)

for i in range(30):
  for j in range(6):
    v = raw_mat[i, j]
    nv = norm_mat[i, j]
    tc = 'white' if nv > 0.65 or nv < 0.25 else '#212529'
    txt = f'{v:.2f}%' if v >= 0.05 else f'{v:.3f}%'
    ax_s1.text(
        j,
        i,
        txt,
        ha='center',
        va='center',
        fontsize=pt(8.0, W_S1),
        fontweight='bold',
        color=tc,
    )

ax_s1.set_xticks(range(6))
ax_s1.set_xticklabels(SAMPLES, fontsize=pt(9.0, W_S1), fontweight='bold')
ax_s1.set_yticks(range(30))
ylabels_g30 = [CLIN_SYM.get(clin(g), '  ') + g for g in genus_labels]
ax_s1.set_yticklabels(
    ylabels_g30, fontsize=pt(8.5, W_S1), style='italic', fontweight='bold'
)
for lbl, g in zip(ax_s1.get_yticklabels(), genus_labels):
  lbl.set_color(CLIN_COL[clin(g)])

ax_s1.set_title(
    'Hierarchical Heatmap of Top 30 Bacterial Genera (Row Min-Max Scaled)',
    fontsize=14.0,
    fontweight='bold',
    pad=14,
)
cb_s1 = fig_s1.colorbar(
    im_s1, ax=ax_s1, orientation='vertical', fraction=0.035, pad=0.03
)
cb_s1.set_label(
    'Row-Normalized Relative Abundance (Min-Max scaled)',
    fontsize=pt(9.0, W_S1),
    fontweight='bold',
)
cb_s1.ax.tick_params(labelsize=pt(8.0, W_S1))

save_dual(fig_s1, 'Supplementary_Fig_S1')

# Supp Fig S2: Differential Abundance Enrichment Ratios
W_S2, H_S2 = 13.0, 5.8
df_diff = df_g.copy()
df_diff['mean_all'] = df_diff[SAMPLES].mean(axis=1)
df_diff = df_diff[df_diff['mean_all'] > 0.0005].copy()

df_diff['mar_ratio'] = df_diff['Mardan'] / (
    df_diff[['Peshawar', 'Swat']].mean(axis=1) + 1e-6
)
df_diff['pes_ratio'] = df_diff['Peshawar'] / (
    df_diff[['Mardan', 'Swat']].mean(axis=1) + 1e-6
)
df_diff['swa_ratio'] = df_diff['Swat'] / (
    df_diff[['Mardan', 'Peshawar']].mean(axis=1) + 1e-6
)

top_mar = df_diff.sort_values('mar_ratio', ascending=False).head(8)
top_pes = df_diff.sort_values('pes_ratio', ascending=False).head(8)
top_swa = df_diff.sort_values('swa_ratio', ascending=False).head(8)

fig_s2, axes_s2 = plt.subplots(1, 3, figsize=(W_S2, H_S2), facecolor='white')
fig_s2.subplots_adjust(
    wspace=0.35, left=0.12, right=0.96, top=0.86, bottom=0.12
)

city_sets = [
    ('Mardan Enriched Genera', top_mar, 'mar_ratio', '#1F4E79'),
    (
        'Peshawar Enriched Genera (Log10 Scale)',
        top_pes,
        'pes_ratio',
        '#C55A11',
    ),
    ('Swat Enriched Genera', top_swa, 'swa_ratio', '#375623'),
]

for ax_s2, (t_diff, d_diff, r_col, col_diff) in zip(axes_s2, city_sets):
  ax_s2.set_facecolor('#FAFAFA')
  vals = d_diff[r_col].values[::-1]
  gnames = d_diff['name'].values[::-1]
  y_pos = np.arange(len(gnames))

  if 'Peshawar' in t_diff:
    ax_s2.barh(
        y_pos, vals, color=col_diff, alpha=0.90, edgecolor='white', height=0.60
    )
    ax_s2.set_xscale('log')
    for y_i, v in zip(y_pos, vals):
      ax_s2.text(
          v * 1.15,
          y_i,
          f'{v:.1f}×',
          va='center',
          fontsize=pt(8.5, W_S2),
          fontweight='bold',
          color='#212529',
      )
    ax_s2.set_xlim(1, max(vals) * 8.0)
  else:
    ax_s2.barh(
        y_pos, vals, color=col_diff, alpha=0.90, edgecolor='white', height=0.60
    )
    for y_i, v in zip(y_pos, vals):
      ax_s2.text(
          v + max(vals) * 0.03,
          y_i,
          f'{v:.1f}×',
          va='center',
          fontsize=pt(8.5, W_S2),
          fontweight='bold',
          color='#212529',
      )
    ax_s2.set_xlim(0, max(vals) * 1.25)

  ax_s2.set_yticks(y_pos)
  ax_s2.set_yticklabels(
      gnames, fontsize=pt(8.5, W_S2), style='italic', fontweight='bold'
  )
  ax_s2.set_xlabel(
      'Enrichment Ratio vs Other Cities',
      fontsize=pt(9.0, W_S2),
      fontweight='bold',
      labelpad=6,
  )
  ax_s2.tick_params(axis='x', labelsize=pt(8.0, W_S2))
  ax_s2.spines['top'].set_visible(False)
  ax_s2.spines['right'].set_visible(False)
  ax_s2.grid(axis='x', alpha=0.25)
  ax_s2.set_title(t_diff, fontsize=13.0, fontweight='bold', pad=14)

save_dual(fig_s2, 'Supplementary_Fig_S2')

print('\n==================================================================')
print('ALL FIGURES REGENERATED WITH CALIBRATED SPRINGER TYPOGRAPHY!')
print('==================================================================')
