# Untreated Slaughterhouse Effluents Harbor High-Risk Zoonotic Pathobionts and Associated Bacteriophages: A Multi-City One Health Metagenomic Assessment

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23077892.svg)](https://doi.org/10.5281/zenodo.23077892)
[![NCBI BioProject](https://img.shields.io/badge/NCBI%20BioProject-PRJNA1218556-green.svg)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1218556)

This repository contains the computational pipelines, ecological modeling workflows, curated taxonomic abundance tables, and reproducible figure generation scripts for:

> **"Untreated Slaughterhouse Effluents Harbor High-Risk Zoonotic Pathobionts and Associated Bacteriophages: A Multi-City One Health Metagenomic Assessment"**  
> *Journal of Applied Microbiology* (Oxford University Press / Applied Microbiology International)

---

## 📖 Study Overview

Untreated slaughterhouse effluent is an important but understudied environmental pathway linking livestock processing to municipal water networks. In low- and middle-income countries, municipal abattoirs frequently discharge effluent into open drainage channels without prior disinfection. 

This study couples **shotgun metagenomic sequencing** across municipal slaughterhouses in three major urban centers of Khyber Pakhtunkhwa, Pakistan (Mardan, Peshawar, and Swat) with **neutral community assembly modeling** (Sloan Neutral Community Model) and **multi-criteria hazard profiling** to:
1. Describe bacterial community composition and core taxonomic architecture in abattoir wastewater.
2. Profile high-risk zoonotic pathobionts (*Campylobacter*, *Arcobacter*, *Salmonella*, *Listeria*, and the *Brucella/Ochrobactrum* group) alongside genera encompassing WHO critical-priority species (*Acinetobacter* and *Pseudomonas*).
3. Test whether neutral dispersal drift or deterministic environmental filtering shapes metacommunity assembly.
4. Establish an operational **Zoonotic Hazard Index (ZHI)** and municipal-specific **Linear One Health Hazard Profiles**.
5. Characterize the effluent virome and identify bacteriophages belonging to the class *Caudoviricetes* linked to priority bacterial hosts.

---

## 📁 Repository Structure

```text
.
├── README.md                           # Project overview, workflow guide, and reproduction steps
├── requirements.txt                    # Python dependency requirements
├── LICENSE                             # MIT open-source license
│
├── data/                               # Curated abundance matrices and virome profiles
│   ├── merged_phylum_relative_abundance.csv   # Phylum-level relative abundance table (6 libraries)
│   ├── merged_genus_relative_abundance.csv    # Genus-level relative abundance table (1,622 genera)
│   └── Slaughterhouse_virome_curated.xlsx     # Curated virome dataset with host linkages and classifications
│
└── code/                               # Reproducible analysis and publication figure generation
    ├── generate_vector_figures.py      # Master script for Figures 1, 2, 3, 4, 6 and S1, S2
    ├── generate_sloan_neutral_and_radar.py # Master script for Figure 5 (Sloan Model & Linear Hazard Profile)
    ├── build_5_composite_figures.py    # Composite multi-panel figure assembler
    ├── Figure1_Phylum_Barplot.py       # Phylum-level relative abundance (Fig. 1A)
    ├── Figure2_Genus_Barplot.py        # Top 20 bacterial genera (Fig. 1B)
    ├── Figure3_Heatmap.py              # Genus hierarchical clustering heatmap (Fig. S1)
    ├── Figure4_Annotation_Barplot.py   # Pathogen counts and relative abundance by category (Fig. 2A-B)
    ├── Figure5_ZoonoticCityBar.py      # City-level zoonotic pathogen profiles (Fig. 2C)
    ├── Figure6_AlphaDiversity.py       # Bacterial alpha diversity metrics (Fig. 3A)
    ├── Figure7_PCoA.py                 # Principal Coordinates Analysis (PCoA) ordination (Fig. 3B)
    ├── Figure8_BrayCurtis.py           # Pairwise Bray-Curtis dissimilarity matrix (Fig. 3C)
    ├── Figure9_ZoonoticEnrichment.py   # Zoonotic enrichment ratios across cities (Fig. S2)
    ├── Figure10_CoOccurrence.py        # Polymicrobial co-occurrence network (Fig. 4C)
    ├── Figure11_Virome.py              # Wastewater virome composition and diversity (Fig. 6)
    └── FigureZHI.py                    # Zoonotic Hazard Index (ZHI) formulation (Fig. 4A-B)
```

---

## 🚀 Getting Started

### 1. Prerequisites
- Python 3.9, 3.10, or 3.11
- Git

### 2. Clone the Repository
```bash
git clone https://github.com/ihtishamnaeem36/Untreated-Slaughterhouse-Effluent-as-a-Hidden-Zoonotic-Reservoir.git
cd Untreated-Slaughterhouse-Effluent-as-a-Hidden-Zoonotic-Reservoir
```

### 3. Create a Virtual Environment & Install Dependencies
```bash
# On Linux / macOS:
python3 -m venv venv
source venv/bin/activate

# On Windows:
python -m venv venv
venv\Scripts\activate

# Install required packages
pip install -r requirements.txt
```

---

## 🔬 Reproducing Publication Figures

All figures can be reproduced directly by running the dedicated scripts within the `code/` folder:

```bash
cd code
```

### Master Scripts (Complete Multi-Panel Figures)
1. **Figures 1, 2, 3, 4, 6 and S1, S2**:
   ```bash
   python generate_vector_figures.py
   ```
   Generates publication-quality 600 DPI raster (`.png`) and vector (`.pdf`) multi-panel figures.

2. **Figure 5 (Sloan Neutral Community Model & Linear One Health Hazard Profile)**:
   ```bash
   python generate_sloan_neutral_and_radar.py
   ```
   Fits the Sloan Neutral Model via non-linear least squares, computes 95% binomial confidence bounds, and generates the composite radar and linear multi-criteria risk profile.

3. **Composite Assembler**:
   ```bash
   python build_5_composite_figures.py
   ```

### Individual Modular Scripts
| Script | Figure | Description |
| :--- | :--- | :--- |
| `Figure1_Phylum_Barplot.py` | Fig. 1A | Phylum-level relative abundance stacked barplot |
| `Figure2_Genus_Barplot.py` | Fig. 1B | Top 20 dominant bacterial genera across cities |
| `Figure3_Heatmap.py` | Fig. S1 | Hierarchical clustering heatmap of top bacterial genera |
| `Figure4_Annotation_Barplot.py` | Fig. 2A–B | Clinical category distribution and relative abundance |
| `Figure5_ZoonoticCityBar.py` | Fig. 2C | High-consequence zoonotic pathogen abundance by city |
| `Figure6_AlphaDiversity.py` | Fig. 3A | Bacterial alpha diversity (Observed, Shannon, Simpson, Chao1) |
| `Figure7_PCoA.py` | Fig. 3B | PCoA ordination based on Bray-Curtis dissimilarity |
| `Figure8_BrayCurtis.py` | Fig. 3C | Pairwise Bray-Curtis dissimilarity heatmap |
| `Figure9_ZoonoticEnrichment.py` | Fig. S2 | Differential abundance and zoonotic enrichment ratios |
| `Figure10_CoOccurrence.py` | Fig. 4C | Polymicrobial co-occurrence network (Spearman $\|r\| \ge 0.8$, $p < 0.05$) |
| `Figure11_Virome.py` | Fig. 6 | Metagenomic virome composition, diversity, and phage–host links |
| `FigureZHI.py` | Fig. 4A–B | Zoonotic Hazard Index (ZHI) and municipal risk tier rankings |

---

## 📊 Analytical Methodology

### Sloan Neutral Community Model (SNCM)
The model tests whether stochastic neutral drift or deterministic niche filtering drives bacterial community assembly. The occurrence frequency of each genus across municipal abattoir metacommunities is modeled as a function of mean relative abundance ($p$) using a non-linear least-squares fit of the beta distribution:

$$\Phi(p) = I_p(N m p, \, N m (1 - p))$$

where:
- $N$ is metacommunity size
- $m$ is the estimated immigration parameter
- Taxa exceeding the 95% binomial confidence interval represent deterministically selected opportunists (e.g., biofilm formers persisting against hydraulic washout), whereas taxa falling below represent deterministically suppressed lineages.

### Linear One Health Hazard Profile & ZHI
The composite Zoonotic Hazard Index (ZHI) integrates:
1. Relative abundance of high-consequence zoonotic pathobionts
2. Representation of genera carrying WHO critical-priority species
3. Virulence and persistence factors (biofilm formation determinants)
4. Environmental community dissimilarity

---

## 🧬 Data Availability

- **NCBI BioProject Accession**: [PRJNA1218556](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1218556)
- **NCBI BioSample Accessions**: SAMN46789121–SAMN46789126
- **NCBI SRA Run Accessions**: SRR32263152–SRR32263157
- **Zenodo Permanent Archive**: [10.5281/zenodo.23077892](https://doi.org/10.5281/zenodo.23077892)
- Curated taxonomic and virome matrices are provided in the [`data/`](./data/) directory of this repository.

---

## 👥 Authors & Contact

- **Ishaq Khan, Ph.D.** (Co-First & Corresponding Author)  
  Department of Biotechnology, Genetics, Forensic and Microbiology, University of Swat, Pakistan  
  Email: `ishaq@uswat.edu.pk` | ORCID: [0009-0000-8265-7898](https://orcid.org/0009-0000-8265-7898)

- **Shehzad Khan** (Co-First Author)  
  Department of Biotechnology, Genetics, Forensic and Microbiology, University of Swat, Pakistan

- **Muhammad Shafiq, Ph.D.** (Co-Corresponding Author)  
  Department of Clinical Pharmacy, Shantou University Medical College, China  
  Email: `drshafiq@stu.edu.cn` | ORCID: [0000-0002-4346-5903](https://orcid.org/0000-0002-4346-5903)

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
