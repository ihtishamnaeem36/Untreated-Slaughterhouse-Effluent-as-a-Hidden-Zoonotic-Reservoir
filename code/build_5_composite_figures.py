# -*- coding: utf-8 -*-
"""
build_5_composite_figures.py

Combines the 11 figures into:
- 5 Main Figures for Current Microbiology (strictly complying with the <= 5 figures limit).
- 2 Supplementary Figures for Supplementary Information (Heatmap and Differential Enrichment).
"""

import os, sys
from PIL import Image, ImageDraw, ImageFont

sys.stdout.reconfigure(encoding='utf-8')

src_dir = 'submission_figures'
out_dir = 'submission_figures_current_microbiology'
os.makedirs(out_dir, exist_ok=True)

# Helper function to add Panel Letter (A, B, C, D) at top-left
def add_panel_label(img, label, font_size=180, pad=40):
    draw = ImageDraw.Draw(img)
    # Draw simple bold text or background circle
    # Using default font scaled or basic text
    try:
        # try Arial or standard TTF if available
        import matplotlib.font_manager as fm
        font_path = fm.findfont(fm.FontProperties(family='DejaVu Sans', weight='bold'))
        font = ImageFont.truetype(font_path, font_size)
    except:
        font = ImageFont.load_default()
    
    # Draw background box for clarity
    bbox = draw.textbbox((pad, pad), label, font=font)
    draw.rectangle([bbox[0]-15, bbox[1]-10, bbox[2]+15, bbox[3]+10], fill='white')
    draw.text((pad, pad), label, fill='black', font=font)
    return img

print("1. Generating Main Figure 1 (Taxonomic Composition: Phyla + Top 20 Genera)...")
im1 = Image.open(os.path.join(src_dir, 'Figure_01_Phylum_Abundance.png')).convert('RGB')
im2 = Image.open(os.path.join(src_dir, 'Figure_02_Top_Genera.png')).convert('RGB')

im1 = add_panel_label(im1, "A")
im2 = add_panel_label(im2, "B")

# Resize im2 to have same height as im1
h1 = im1.size[1]
w2_new = int(im2.size[0] * (h1 / im2.size[1]))
im2_resized = im2.resize((w2_new, h1), Image.Resampling.LANCZOS)

# Create side-by-side composite
fig1_composite = Image.new('RGB', (im1.size[0] + w2_new + 60, h1), color='white')
fig1_composite.paste(im1, (0, 0))
fig1_composite.paste(im2_resized, (im1.size[0] + 60, 0))
fig1_composite.save(os.path.join(out_dir, 'Fig_1.png'), dpi=(600, 600))
print("Saved Fig_1.png successfully!")

print("2. Generating Main Figure 2 (Priority Pathogen Annotation: Categories + City Abundance)...")
im4 = Image.open(os.path.join(src_dir, 'Figure_04_Clinical_Annotation.png')).convert('RGB')
im5 = Image.open(os.path.join(src_dir, 'Figure_05_Zoonotic_City_Abundance.png')).convert('RGB')

im4 = add_panel_label(im4, "A")
im5 = add_panel_label(im5, "B")

# Stack vertically or side-by-side: im4 is wide (10338x4809), im5 is tall (6542x7726)
# Resize im5 to match width or side-by-side:
# Let's resize im5 height to match im4 height
h4 = im4.size[1]
w5_new = int(im5.size[0] * (h4 / im5.size[1]))
im5_resized = im5.resize((w5_new, h4), Image.Resampling.LANCZOS)

fig2_composite = Image.new('RGB', (im4.size[0] + w5_new + 60, h4), color='white')
fig2_composite.paste(im4, (0, 0))
fig2_composite.paste(im5_resized, (im4.size[0] + 60, 0))
fig2_composite.save(os.path.join(out_dir, 'Fig_2.png'), dpi=(600, 600))
print("Saved Fig_2.png successfully!")

print("3. Generating Main Figure 3 (Community Diversity: Alpha Diversity + Beta Diversity)...")
im6 = Image.open(os.path.join(src_dir, 'Figure_06_Alpha_Diversity.png')).convert('RGB')
im7 = Image.open(os.path.join(src_dir, 'Figure_07_Beta_Diversity_PCoA.png')).convert('RGB')

im6 = add_panel_label(im6, "A")
# im7 already has 7A and 7B, let's keep or add B
im7 = add_panel_label(im7, "B")

# Stack vertically: both are wide (~10500 width)
w6 = im6.size[0]
h7_new = int(im7.size[1] * (w6 / im7.size[0]))
im7_resized = im7.resize((w6, h7_new), Image.Resampling.LANCZOS)

fig3_composite = Image.new('RGB', (w6, im6.size[1] + h7_new + 60), color='white')
fig3_composite.paste(im6, (0, 0))
fig3_composite.paste(im7_resized, (0, im6.size[1] + 60))
fig3_composite.save(os.path.join(out_dir, 'Fig_3.png'), dpi=(600, 600))
print("Saved Fig_3.png successfully!")

print("4. Generating Main Figure 4 (Hazard Modeling & Ecological Network: ZHI + Network)...")
im9 = Image.open(os.path.join(src_dir, 'Figure_09_Zoonotic_Hazard_Index.png')).convert('RGB')
im10 = Image.open(os.path.join(src_dir, 'Figure_10_CoOccurrence_Network.png')).convert('RGB')

im9 = add_panel_label(im9, "A")
im10 = add_panel_label(im10, "B")

# Stack vertically
w10 = im10.size[0]
h9_new = int(im9.size[1] * (w10 / im9.size[0]))
im9_resized = im9.resize((w10, h9_new), Image.Resampling.LANCZOS)

fig4_composite = Image.new('RGB', (w10, h9_new + im10.size[1] + 60), color='white')
fig4_composite.paste(im9_resized, (0, 0))
fig4_composite.paste(im10, (0, h9_new + 60))
fig4_composite.save(os.path.join(out_dir, 'Fig_4.png'), dpi=(600, 600))
print("Saved Fig_4.png successfully!")

print("5. Generating Main Figure 5 (Wastewater Virome: already 4 panels A, B, C, D)...")
im11 = Image.open(os.path.join(src_dir, 'Figure_11_Virome_Composition.png')).convert('RGB')
im11.save(os.path.join(out_dir, 'Fig_5.png'), dpi=(600, 600))
print("Saved Fig_5.png successfully!")

print("6. Copying Supplementary Figures (Heatmap & Differential Enrichment)...")
im3 = Image.open(os.path.join(src_dir, 'Figure_03_Heatmap_Top30.png')).convert('RGB')
im3.save(os.path.join(out_dir, 'Supplementary_Fig_S1.png'), dpi=(600, 600))

im8 = Image.open(os.path.join(src_dir, 'Figure_08_Differential_Enrichment.png')).convert('RGB')
im8.save(os.path.join(out_dir, 'Supplementary_Fig_S2.png'), dpi=(600, 600))

print("ALL 5 MAIN FIGURES + 2 SUPPLEMENTARY FIGURES GENERATED SUCCESSFULLY!")
