# -*- coding: utf-8 -*-
"""
generate_sloan_neutral_and_radar.py

Generates Figure 5 for Journal of Applied Microbiology:
Panel A: Sloan Neutral Community Model (SNCM) Fit (Nonlinear curve fit, 95% binomial CI, key zoonotic/AMR pathogen callouts).
Panel B: Community Assembly Mechanics Partitioning (Horizontal Assembly Breakdown Bar with structured mechanism cards).
Panel C: Linear Multi-Criteria One Health Hazard Profile (Grouped Horizontal Barplot across 6 dimensions for Peshawar, Swat, Mardan).

Outputs 600 DPI PNG and vector PDF.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from scipy.optimize import curve_fit
from scipy.special import betainc

sys.stdout.reconfigure(encoding='utf-8')

# Global publication typography
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9.0
plt.rcParams['axes.linewidth'] = 0.9
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8

# ------------------------------------------------------------------------------
# 1. LOAD DATA & COMPUTE EXACT SLOAN NEUTRAL MODEL
# ------------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXTRA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "extra"))

def find_file(fname):
    candidates = [
        os.path.join(SCRIPT_DIR, fname),
        os.path.join(SCRIPT_DIR, "..", "data", fname),
        os.path.join(SCRIPT_DIR, "..", "extra", fname),
        os.path.join(EXTRA_DIR, fname),
        os.path.join(EXTRA_DIR, "all cities", "genus", "genus relative abundance", fname),
        os.path.join(EXTRA_DIR, "Virome Section", fname),
        fname
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return fname

def find_file(fname):
    candidates = [
        os.path.join(SCRIPT_DIR, fname),
        os.path.join(SCRIPT_DIR, "..", "data", fname),
        os.path.join(SCRIPT_DIR, "..", "extra", fname),
        os.path.join(EXTRA_DIR, fname),
        os.path.join(EXTRA_DIR, "all cities", "genus", "genus relative abundance", fname),
        os.path.join(EXTRA_DIR, "Virome Section", fname),
        fname
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return fname

genus_csv = find_file("merged_genus_relative_abundance.csv")
virome_xlsx = find_file("Slaughterhouse_virome_curated.xlsx")
df_genus = pd.read_csv(genus_csv)
virome_df = pd.read_excel(virome_xlsx)

sample_cols = ['MSLW1', 'MSLW2', 'PSLW1', 'PSLW2', 'SSLW1', 'SSLW2']
df_genus['mean'] = df_genus[sample_cols].mean(axis=1)
df_genus['freq'] = (df_genus[sample_cols] > 0).sum(axis=1) / len(sample_cols)

sub = df_genus[df_genus['mean'] > 0].copy().reset_index(drop=True)
p = sub['mean'].values
freq = sub['freq'].values
N_reads_avg = 134800000 / 6  # average reads per sample = ~22.5 million

def fit_sncm(p_arr, Nm):
    limit = 1.0 / N_reads_avg
    return 1.0 - betainc(Nm * p_arr, Nm * (1.0 - p_arr), limit)

popt, _ = curve_fit(fit_sncm, p, freq, p0=[1000], bounds=(1, 1e7), maxfev=10000)
Nm_est = popt[0]
pred_freq = fit_sncm(p, Nm_est)

# Goodness of fit (R^2)
ss_res = np.sum((freq - pred_freq) ** 2)
ss_tot = np.sum((freq - np.mean(freq)) ** 2)
r2_sncm = 1.0 - (ss_res / ss_tot)

# 95% Confidence intervals
se = np.sqrt(pred_freq * (1.0 - pred_freq) / len(sample_cols))
upper_ci = np.clip(pred_freq + 1.96 * se, 0, 1.0)
lower_ci = np.clip(pred_freq - 1.96 * se, 0, 1.0)

sub['pred_freq'] = pred_freq
sub['upper_ci'] = upper_ci
sub['lower_ci'] = lower_ci
sub['status'] = np.where(freq > upper_ci, 'Above', np.where(freq < lower_ci, 'Below', 'Neutral'))

n_total = len(sub)
n_neut = int((sub['status'] == 'Neutral').sum())
n_above = int((sub['status'] == 'Above').sum())
n_below = int((sub['status'] == 'Below').sum())

pct_neut = n_neut / n_total * 100
pct_above = n_above / n_total * 100
pct_below = n_below / n_total * 100

# ------------------------------------------------------------------------------
# 2. ONE HEALTH HAZARD PROFILE METRICS
# ------------------------------------------------------------------------------
crit_mardan = df_genus[df_genus['name'].isin(['Pseudomonas', 'Acinetobacter'])][['MSLW1', 'MSLW2']].sum().mean() * 100
crit_pesh = df_genus[df_genus['name'].isin(['Pseudomonas', 'Acinetobacter'])][['PSLW1', 'PSLW2']].sum().mean() * 100
crit_swat = df_genus[df_genus['name'].isin(['Pseudomonas', 'Acinetobacter'])][['SSLW1', 'SSLW2']].sum().mean() * 100

ent_taxa = ['Arcobacter', 'Aliarcobacter', 'Campylobacter', 'Salmonella', 'Listeria']
ent_mardan = df_genus[df_genus['name'].isin(ent_taxa)][['MSLW1', 'MSLW2']].sum().mean() * 100
ent_pesh = df_genus[df_genus['name'].isin(ent_taxa)][['PSLW1', 'PSLW2']].sum().mean() * 100
ent_swat = df_genus[df_genus['name'].isin(ent_taxa)][['SSLW1', 'SSLW2']].sum().mean() * 100

zhi_mardan = 75.71
zhi_pesh = 130.84
zhi_swat = 80.74

hfv = virome_df[virome_df['Viral species'].str.contains('human fecal virus', case=False, na=False)]
hfv_mardan = hfv[['MSLW1', 'MSLW2']].mean(axis=1).values[0]
hfv_pesh = hfv[['PSLW1', 'PSLW2']].mean(axis=1).values[0]
hfv_swat = hfv[['SSLW1', 'SSLW2']].mean(axis=1).values[0]

rich_mardan = 22.5
rich_pesh = 24.5
rich_swat = 23.5

phage_rows = virome_df[virome_df['Viral species'].str.contains('phage', case=False, na=False) & ~virome_df['Viral species'].str.contains('human', case=False, na=False)]
phage_mardan = phage_rows[['MSLW1', 'MSLW2']].sum().mean()
phage_pesh = phage_rows[['PSLW1', 'PSLW2']].sum().mean()
phage_swat = phage_rows[['SSLW1', 'SSLW2']].sum().mean()

metrics_display = [
    'Critical AMR\nLineage Burden (%)',
    'Enteric Zoonosis\nShedding (%)',
    'Zoonotic Hazard\nIndex (ZHI)',
    'Domestic Sewage\nInfiltration (%)',
    'Priority Pathogen\nRichness (Count)',
    'Active Viral Predation\nPressure (%)'
]

raw_data = {
    'Peshawar': [crit_pesh, ent_pesh, zhi_pesh, hfv_pesh, rich_pesh, phage_pesh],
    'Swat': [crit_swat, ent_swat, zhi_swat, hfv_swat, rich_swat, phage_swat],
    'Mardan': [crit_mardan, ent_mardan, zhi_mardan, hfv_mardan, rich_mardan, phage_mardan]
}

m_norm, p_norm, s_norm = [], [], []
for i in range(len(metrics_display)):
    vals = [raw_data['Mardan'][i], raw_data['Peshawar'][i], raw_data['Swat'][i]]
    mn, mx = min(vals), max(vals)
    if mx > mn:
        m_norm.append(0.20 + 0.80 * (vals[0] - mn) / (mx - mn))
        p_norm.append(0.20 + 0.80 * (vals[1] - mn) / (mx - mn))
        s_norm.append(0.20 + 0.80 * (vals[2] - mn) / (mx - mn))
    else:
        m_norm.append(1.0)
        p_norm.append(1.0)
        s_norm.append(1.0)

# Colors tailored to Oxford / JAM standards
C_PESH = '#B03A2E'  # Crimson Red
C_SWAT = '#2874A6'  # Slate Blue
C_MARD = '#D35400'  # Ochre Orange
C_NEUT = '#1F77B4'  # Deep Ocean Blue
C_ABOV = '#D62728'  # Carmine Red
C_BELO = '#FF7F0E'  # Amber Orange
C_FIT  = '#2CA02C'  # Forest Green

# ------------------------------------------------------------------------------
# 3. MASTER CANVAS & EXPLICIT SUBPLOT GEOMETRY
# ------------------------------------------------------------------------------
fig = plt.figure(figsize=(16.8, 12.0), facecolor='white')

# Left column: Panel A (Sloan Neutral Model)
ax_a = fig.add_axes([0.07, 0.08, 0.39, 0.84])

# Right column:
# Panel B: 4-Bar Municipal Comparative Breakdown
ax_b = fig.add_axes([0.58, 0.58, 0.39, 0.31])

# Panel C barplot:
ax_c = fig.add_axes([0.58, 0.08, 0.39, 0.42])

# ------------------------------------------------------------------------------
# PANEL A: SLOAN NEUTRAL MODEL SCATTER
# ------------------------------------------------------------------------------
p_neutral = sub[sub['status'] == 'Neutral']
p_above = sub[sub['status'] == 'Above']
p_below = sub[sub['status'] == 'Below']

ax_a.scatter(p_neutral['mean'], p_neutral['freq'], color=C_NEUT, alpha=0.45, s=26, edgecolors='none', label=f'Neutrally Dispersed ({pct_neut:.1f}%, n = {n_neut:,})', zorder=2)
ax_a.scatter(p_above['mean'], p_above['freq'], color=C_ABOV, alpha=0.90, s=42, edgecolors='black', linewidth=0.5, label=f'Positively Selected / Persistent ({pct_above:.1f}%, n = {n_above:,})', zorder=4)
ax_a.scatter(p_below['mean'], p_below['freq'], color=C_BELO, alpha=0.85, s=34, edgecolors='#333333', linewidth=0.4, label=f'Negatively Filtered / Transient ({pct_below:.1f}%, n = {n_below:,})', zorder=3)

sort_idx = np.argsort(p)
p_sorted = p[sort_idx]
pred_sorted = pred_freq[sort_idx]
up_sorted = upper_ci[sort_idx]
low_sorted = lower_ci[sort_idx]

ax_a.plot(p_sorted, pred_sorted, color=C_FIT, lw=2.4, label='Neutral Fit (SNCM)', zorder=5)
ax_a.plot(p_sorted, up_sorted, color=C_FIT, lw=1.2, ls='--', alpha=0.85, zorder=5)
ax_a.plot(p_sorted, low_sorted, color=C_FIT, lw=1.2, ls='--', alpha=0.85, label='95% Binomial CI', zorder=5)

ax_a.set_xscale('log')
ax_a.set_xlabel('Metacommunity Mean Relative Abundance (log₁₀ scale)', fontsize=10.0, fontweight='bold', labelpad=8)
ax_a.set_ylabel('Taxon Occurrence Frequency across Metacommunity (n = 6)', fontsize=10.0, fontweight='bold', labelpad=8)
ax_a.set_ylim(-0.04, 1.12)
ax_a.set_title(f'Sloan Neutral Community Model Fit (R² = {r2_sncm:.3f}, Nm = {Nm_est:,.0f})', fontsize=11.5, fontweight='bold', pad=12, loc='left')
ax_a.legend(loc='lower right', fontsize=8.2, framealpha=0.96, edgecolor='#D0D0D0')
ax_a.grid(True, linestyle=':', alpha=0.45)
ax_a.text(-0.11, 1.03, 'A', transform=ax_a.transAxes, fontsize=15, fontweight='bold', va='bottom')

# Distinct pathogen callouts
callouts = {
    'Campylobacter': dict(offset=(-55, 30), label='Campylobacter\n(Enteric Zoonosis)'),
    'Salmonella':    dict(offset=(-45, 28), label='Salmonella\n(Enteric Zoonosis)'),
    'Acinetobacter': dict(offset=(25, 28),  label='Acinetobacter\n(Critical AMR)'),
    'Pseudomonas':   dict(offset=(-60, -36), label='Pseudomonas\n(Critical AMR)'),
}
for taxon, c_info in callouts.items():
    m = sub[sub['name'] == taxon]
    if not m.empty:
        r = m.iloc[0]
        ax_a.scatter(r['mean'], r['freq'], color='#FFF176', edgecolor='#111111', s=70, zorder=6, linewidth=1.2)
        ax_a.annotate(
            c_info['label'],
            xy=(r['mean'], r['freq']),
            xytext=c_info['offset'],
            textcoords='offset points',
            fontsize=7.8,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#FFFFFF', alpha=0.95, edgecolor='#444444', linewidth=0.6),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.08', color='#222222', lw=0.9),
            zorder=7
        )

# ------------------------------------------------------------------------------
# PANEL B: CITY-STRATIFIED ASSEMBLY BREAKDOWN (4 BARS, PURE DATA FIGURE)
# ------------------------------------------------------------------------------
ax_b.set_facecolor('#FAFAFA')

cities_b = [
    'Mardan\n(n = 1,034)',
    'Swat\n(n = 1,053)',
    'Peshawar\n(n = 1,192)',
    'Metacommunity\n(n = 1,372)'
]
neut_pcts = [73.2, 73.7, 69.0, 69.2]
above_pcts = [8.0, 7.0, 11.7, 12.4]
below_pcts = [18.8, 19.3, 19.3, 18.4]
neut_counts = [757, 776, 822, 949]
above_counts = [83, 74, 140, 170]
below_counts = [194, 203, 230, 253]

y_pos_b = np.arange(len(cities_b))
bar_h_b = 0.52

b_neut = ax_b.barh(y_pos_b, neut_pcts, height=bar_h_b, color=C_NEUT, edgecolor='white', linewidth=1.5, label='Neutrally Dispersed')
b_above = ax_b.barh(y_pos_b, above_pcts, left=neut_pcts, height=bar_h_b, color=C_ABOV, edgecolor='white', linewidth=1.5, label='Positively Selected / Persistent')
b_below = ax_b.barh(y_pos_b, below_pcts, left=np.array(neut_pcts) + np.array(above_pcts), height=bar_h_b, color=C_BELO, edgecolor='white', linewidth=1.5, label='Negatively Filtered / Transient')

for i in range(len(cities_b)):
    ax_b.text(neut_pcts[i] / 2, y_pos_b[i], f"{neut_pcts[i]:.1f}%\n({neut_counts[i]})", ha='center', va='center', color='white', fontweight='bold', fontsize=8.2)
    ax_b.text(neut_pcts[i] + above_pcts[i] / 2, y_pos_b[i], f"{above_pcts[i]:.1f}%\n({above_counts[i]})", ha='center', va='center', color='white', fontweight='bold', fontsize=7.6)
    ax_b.text(neut_pcts[i] + above_pcts[i] + below_pcts[i] / 2, y_pos_b[i], f"{below_pcts[i]:.1f}%\n({below_counts[i]})", ha='center', va='center', color='white', fontweight='bold', fontsize=7.6)

ax_b.set_yticks(y_pos_b)
ax_b.set_yticklabels(cities_b, fontsize=8.6, fontweight='bold')
ax_b.set_xlim(0, 100)
ax_b.set_xlabel('Percentage of Genera Partitioned by Assembly Category (%)', fontsize=9.2, fontweight='bold', labelpad=5)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_color('#888888')
ax_b.spines['bottom'].set_color('#888888')
ax_b.grid(axis='x', linestyle=':', alpha=0.35)

ax_b.set_title('Community Assembly Mechanics Partitioning across Municipalities', fontsize=11.5, fontweight='bold', pad=32, loc='left')
# Standard compact chart legend:
ax_b.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, frameon=False, fontsize=8.0, handletextpad=0.4, columnspacing=0.8)
ax_b.text(-0.06, 1.15, 'B', transform=ax_b.transAxes, fontsize=15, fontweight='bold', va='bottom')

# ------------------------------------------------------------------------------
# PANEL C: LINEAR MULTI-CRITERIA ONE HEALTH HAZARD PROFILE (GROUPED BARPLOT)
# ------------------------------------------------------------------------------
ax_c.set_facecolor('#FAFAFA')
y_pos = np.arange(len(metrics_display))
bar_h = 0.24

clean_labels_rev = metrics_display[::-1]
p_norm_rev = p_norm[::-1]
s_norm_rev = s_norm[::-1]
m_norm_rev = m_norm[::-1]
raw_p_rev = raw_data['Peshawar'][::-1]
raw_s_rev = raw_data['Swat'][::-1]
raw_m_rev = raw_data['Mardan'][::-1]

b_p = ax_c.barh(y_pos + bar_h, p_norm_rev, height=bar_h, color=C_PESH, label='Peshawar (High Risk)', edgecolor='white', linewidth=0.6)
b_s = ax_c.barh(y_pos, s_norm_rev, height=bar_h, color=C_SWAT, label='Swat (Moderate Risk)', edgecolor='white', linewidth=0.6)
b_m = ax_c.barh(y_pos - bar_h, m_norm_rev, height=bar_h, color=C_MARD, label='Mardan (Baseline / Moderate)', edgecolor='white', linewidth=0.6)

ax_c.set_yticks(y_pos)
ax_c.set_yticklabels(clean_labels_rev, fontsize=8.5, fontweight='bold')
ax_c.set_xlim(0, 1.25)
ax_c.set_xlabel('Normalized One Health Hazard Score (0.00 – 1.00 Scale)', fontsize=9.5, fontweight='bold', labelpad=6)
ax_c.grid(axis='x', linestyle=':', alpha=0.45)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.legend(loc='lower right', fontsize=8.4, framealpha=0.95, edgecolor='#D0D0D0')
ax_c.set_title('Linear Multi-Criteria One Health Hazard Profile', fontsize=11.5, fontweight='bold', pad=12, loc='left')
ax_c.text(-0.06, 1.04, 'C', transform=ax_c.transAxes, fontsize=15, fontweight='bold', va='bottom')

# Exact value annotations on bars
for i in range(len(metrics_display)):
    idx_orig = len(metrics_display) - 1 - i
    if idx_orig in [0, 1, 3, 5]:
        lp = f'{raw_p_rev[i]:.2f}%'
        ls = f'{raw_s_rev[i]:.2f}%'
        lm = f'{raw_m_rev[i]:.2f}%'
    else:
        lp = f'{raw_p_rev[i]:.1f}'
        ls = f'{raw_s_rev[i]:.1f}'
        lm = f'{raw_m_rev[i]:.1f}'
    ax_c.text(p_norm_rev[i] + 0.015, y_pos[i] + bar_h, lp, va='center', fontsize=7.4, fontweight='bold', color=C_PESH)
    ax_c.text(s_norm_rev[i] + 0.015, y_pos[i], ls, va='center', fontsize=7.4, fontweight='bold', color=C_SWAT)
    ax_c.text(m_norm_rev[i] + 0.015, y_pos[i] - bar_h, lm, va='center', fontsize=7.4, fontweight='bold', color=C_MARD)

# ------------------------------------------------------------------------------
# 4. SAVE HIGH-RES PNG & VECTOR PDF TO BOTH SUBMISSION ROOT & CODE DIRECTORY
# ------------------------------------------------------------------------------
out_dir = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
png_path = os.path.join(out_dir, "Figure_5.png")
pdf_path = os.path.join(out_dir, "Figure_5.pdf")

local_png = os.path.join(SCRIPT_DIR, "Fig_5.png")
local_pdf = os.path.join(SCRIPT_DIR, "Fig_5.pdf")

# Root submission folder
fig.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(pdf_path, dpi=600, bbox_inches='tight', facecolor='white')

# Code folder
fig.savefig(local_png, dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(local_pdf, dpi=600, bbox_inches='tight', facecolor='white')

plt.close(fig)
print(f"SUCCESS: Generated 600 DPI publication master:")
print(f"  -> {png_path}")
print(f"  -> {pdf_path}")
print(f"  -> {local_png}")
print(f"  -> {local_pdf}")
