# SMESH — figure and table scripts

Figures and supplementary tables in *Structure-learning meta-analysis resolves cross-context heterogeneity in microbiome associations* is produced by a script in this folder.

---

## Main figures

| Script | Purpose | Output |
|---|---|---|
| `Figure2.R` | Simulation: SMESH against the benchmark methods — cluster recovery (ARI) and signature precision/recall as the signal fraction varies | `Simulation/Figure/Fig2_{b,c,d}.png` |
| `Figure3.R` | Simulation: SMESH given four different input summary statistics (PALM, ANCOM-BC2, LinDA, MaAsLin3) — how much the choice of back-end matters | `Simulation/Figure/Fig3_{A1,A2,B}.png` |
| `Figure4.R` | Pan-disease (GMrepo, 12 diseases): cluster-level effects, leave-one-context-out consensus, method comparison | `GMrepo_analysis/Figure/Figure{A1,A2,B,C,D,E}.png` |
| `Figure5.R` | Pan-tumor (HGMT, 13 tumour types): the same five panels | `HGMT_analysis/Figure/Fig1_{A1,A2}.png`, `Figure{B,C,D,E}.png` |
| `Figure6.R` | Colorectal neoplasia (26 study-stage strata): the same five panels | `CRC_analysis/Figure/Figure{A1,A2,B,C,D,E}.png` |

## Supplementary figures

| Script | Purpose | Output |
|---|---|---|
| `FigureS1.R` | Simulation: precision and recall separately for shared and context-dependent signatures, benchmark methods | `Simulation/Figure/SFig1_{A,B}.png` |
| `FigureS2.R` | The same split, across the four input summary types | `Simulation/Figure/SFig2_{A,B}.png` |
| `FigureS3.R` | Pan-disease stability: selected G per LOCO round, retention/stability/ARI, and the filtering and taxonomy sensitivity analyses | `GMrepo_analysis/Figure/SFig6_{A,B,C,D}.png` |
| `FigureS4.R` | Pan-tumor stability: the same, without the sensitivity panel | `HGMT_analysis/Figure/SFig6_{A,B,C,D,E}.png` |
| `FigureS5.R` | Colorectal neoplasia stability: as above | `CRC_analysis/Figure/SFig6_{A,B,C,D,E}.png` |

## Supplementary tables

| Script | Purpose | Output |
|---|---|---|
| `Make_TableS1.R` | Study list, pan-disease: project, disease, group sizes, estimable taxa | `Figure/TableS1_pan_disease.{tex,csv}` |
| `Make_TableS2.R` | Study list, pan-tumor | `Figure/TableS2_pan_tumor.{tex,csv}` |
| `Make_TableS3.R` | Study list, colorectal neoplasia, by study and stage | `Figure/TableS3_colorectal.{tex,csv}` |
| `Make_TableS4_S6.R` | SMESH-PALM cluster-level effect estimates, one tab per application | `Figure/TableS4_S6_SMESH_effects.xlsx` |
