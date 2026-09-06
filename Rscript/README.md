# SMESH — figure and table scripts

Every figure and supplementary table in *Structure-learning microbiome meta-analysis reveals
shared and context-dependent microbial signatures* is produced by a script in this folder.

These scripts **only draw**. They read results that already exist; they do not fit any model.
The pipeline that produces those results is in [`../Analysis/`](../Analysis) — see its
[README](../Analysis/README.md) — and must have been run first.

---

## Before you start

**Run everything from the repository root, not from this folder.** All paths are written
relative to the root.

```bash
cd /path/to/MOSIAC          # repository root
Rscript Rscript/Figure4.R
```

Each script is self-contained and takes no arguments — one script, one figure or table.

### Dependencies

```r
install.packages(c("tidyverse", "dplyr", "tidyr", "tibble", "stringr", "purrr",
                   "ggplot2", "patchwork", "ggpubr", "scales", "ggtext", "latex2exp",
                   "clue", "mclust", "openxlsx"))
```

Eight of the fifteen scripts also source [`../utility/heatmap_util.R`](../utility/heatmap_util.R)
— `Figure4`, `Figure5`, `Figure6`, `FigureS3`, `FigureS4`, `FigureS5`, `Make_TableS4_S6`,
`Check_cluster_covariates` — which holds the shared plotting code: the cluster palette, the heat-map
panels, the display ordering rules and the label-matching helpers. Nothing else needs
installing.

The table scripts write into `Figure/` at the repository root; create it if it does not exist.

---

## Main figures

| Script | Purpose | Output |
|---|---|---|
| `Figure2.R` | Simulation: SMESH against the benchmark methods — cluster recovery (ARI) and signature precision/recall as the signal fraction varies | `Simulation/Figure/Fig2_{b,c,d}.png` |
| `Figure3.R` | Simulation: SMESH given four different input summary statistics (PALM, ANCOM-BC2, LinDA, MaAsLin3) — how much the choice of back-end matters | `Simulation/Figure/Fig3_{A1,A2,B}.png` |
| `Figure4.R` | Pan-disease (GMrepo, 12 diseases): cluster-level effects, leave-one-context-out consensus, method comparison | `GMrepo_analysis/Figure/Figure{A1,A2,B,C,D,E}.png` |
| `Figure5.R` | Pan-tumor (HGMT, 15 tumour types): the same five panels | `HGMT_analysis/Figure/Fig1_{A1,A2}.png`, `Figure{B,C,D,E}.png` |
| `Figure6.R` | Colorectal neoplasia (26 study-stage strata): the same five panels | `CRC_analysis/Figure/Figure{A1,A2,B,C,D,E}.png` |

Figure 1 is a conceptual schematic and is not generated from data.

**The five panels of Figures 4–6 are the same in each application:**

- **A1 / A2** — cluster-level effect heat maps: the signatures SMESH selects, grouped by which
  clusters share them.
- **B** — leave-one-context-out consensus matrix: how often two contexts co-cluster when a
  third is held out. Drawn at the fixed G reported in the paper.
- **C** — the same contexts under each competing method, labels aligned to SMESH.
- **D** — number of signatures each method selects, split into all-cluster-shared and
  context-dependent.
- **E** — shared-versus-specific signature counts, with Jaccard overlap against SMESH-PALM.

## Supplementary figures

| Script | Purpose | Output |
|---|---|---|
| `FigureS1.R` | Simulation: precision and recall separately for shared and context-dependent signatures, benchmark methods | `Simulation/Figure/SFig1_{A,B}.png` |
| `FigureS2.R` | The same split, across the four input summary types | `Simulation/Figure/SFig2_{A,B}.png` |
| `FigureS3.R` | Pan-disease stability: selected G per LOCO round, retention/stability/ARI, and the filtering and taxonomy sensitivity analyses | `GMrepo_analysis/Figure/SFig6_{A,B,C,D}.png` |
| `FigureS4.R` | Pan-tumor stability: the same, without the sensitivity panel | `HGMT_analysis/Figure/SFig6_{A,B,C,D,E}.png` |
| `FigureS5.R` | Colorectal stability: as above | `CRC_analysis/Figure/SFig6_{A,B,C,D,E}.png` |

**Panels b and c differ in how G is handled, and it matters when reading them:** panel **b**
re-selects G in every leave-one-context-out round and reports the distribution; panel **c**
holds G fixed at the full-data value so that retention, stability and ARI compare partitions
of the same resolution.

## Supplementary tables

| Script | Purpose | Output |
|---|---|---|
| `Make_TableS1.R` | Study list, pan-disease: project, disease, group sizes, estimable taxa | `Figure/TableS1_pan_disease.{tex,csv}` |
| `Make_TableS2.R` | Study list, pan-tumor | `Figure/TableS2_pan_tumor.{tex,csv}` |
| `Make_TableS3.R` | Study list, colorectal, by study and stage | `Figure/TableS3_colorectal.{tex,csv}` |
| `Make_TableS4_S6.R` | SMESH-PALM cluster-level effect estimates, one tab per application | `Figure/TableS4_S6_SMESH_effects.xlsx` |

`Make_TableS1`–`S3` write a `.tex` fragment ready to `\input`, plus the same content as CSV.

`Make_TableS4_S6` writes **one workbook with exactly three tabs** — "Supplementary table 4"
(pan-disease), "Supplementary table 5" (pan-tumor), "Supplementary table 6" (colorectal) — and
no README tab. Each tab has one row per selected signature, in the order that signature is
drawn along the x axis of panel A of the corresponding main figure, and the columns
`Feature ID`, `Shared type`, then one effect column per cluster (`C1` … `CG`, display order).
A cluster in which no member context could estimate the feature is written as the literal
string `NA`, so it stays distinguishable from a tested-but-zero effect; those are the cells the
figures leave white.

## Supporting checks

`Check_cluster_covariates.R` is **not** a table script — nothing it writes is shipped with the
paper. It tests whether cluster membership is explained by context-level characteristics
(study count, sample size, case fraction, sex, age, library size, country, platform, estimable
features) and backs the statement, made once per application in the main text, that the
clusters are not a restatement of study metadata. It writes
`Figure/cluster_covariates.xlsx`, a multi-sheet workbook with a README tab describing its
conventions.

It is by far the slowest script here. Where a context set is small enough it enumerates every
labelling exactly (924 for pan-disease); otherwise it runs 20,000 Monte Carlo permutations per
variable, which is the case for pan-tumor and colorectal. Expect several minutes.

---

## What each script needs to exist first

| Script | Requires | Produced by |
|---|---|---|
| `Figure2.R`, `FigureS1.R` | `Simulation/Sim_CRC_stage1/` | `Analysis/1_simulation_stage1.R` |
| `Figure3.R`, `FigureS2.R` | `Simulation/Sim_CRC_stage2/` | `Analysis/2_simulation_stage2.R` |
| `Figure4.R`, `FigureS3.R` | `GMrepo_analysis/{Data,GMrepo_loso,Model}/` | steps 0, 3, 4, 6 |
| `Figure5.R`, `FigureS4.R` | `HGMT_analysis/{Data,HGMT_loso,Model}/` | steps 0, 3, 4, 6 |
| `Figure6.R`, `FigureS5.R` | `CRC_analysis/{Data,CRC_loso,Model}/` | steps 0, 3, 4, 6 |
| `Make_TableS1.R` | `GMrepo_analysis/Data/` — `summary_stat.rds`, `covariate.*.rds`, `GMrepo_project_summary.xlsx` | steps 0 and 3 |
| `Make_TableS2.R` | `HGMT_analysis/Data/HGMT_16S.Rdata` | step 0 |
| `Make_TableS3.R` | `CRC_analysis/Data/CRC_strata.Rdata`, `CRC_analysis/CRC_loso/Model2_SMESH_s100.Rdata` | steps 0 and 4 |
| `Make_TableS4_S6.R` | all three `Data/` and `*_loso/` folders | steps 3 and 4 |
| `Check_cluster_covariates.R` | all three `Data/` and `*_loso/` folders | steps 0, 3, 4 |

Step numbers refer to the scripts in [`../Analysis/`](../Analysis): `0_preprocessing`,
`3_summary_statistics`, `4_real_data_analysis`, `5_gmrepo_sensitivity`, `6_comparison_models`.

The `*_loso/` folders hold one file per `(G, method, s)` combination, written by
`Analysis/4_real_data_analysis.R`. **`s = 100` is the fit on all contexts**; `s = 1 … L` are
the leave-one-context-out refits.

`Model/` holds the fixed-effect fits for the two-step competitors (`SKMean_FE_G<G>.Rdata`,
`SHC_FE_G<G>.Rdata`) and `Melody_model.Rdata`, written by `Analysis/6_comparison_models.R`
once the number of clusters has been selected.

`FigureS3.R` additionally reads the pan-disease sensitivity results
`GMrepo_analysis/Sensitivity/res_species_X4_G2.Rdata` and `res_genus_G2.Rdata` for its panel D,
written by `Analysis/5_gmrepo_sensitivity.R`, and the genus lookup
`Data/GMrepo/taxon_id_genus_lookup.csv`.

---

## Two things that will surprise you

**Each script writes one PNG per panel, not one file per figure.** `Figure4.R` produces six
images; the published figure is assembled from them. The letters in the file names are the
panel letters in the paper.

**The output names do not always match the figure numbers.** They are historical:

- all three supplementary scripts write `SFig6_*.png`, whatever figure they actually produce;
- `Figure5.R` writes its panel A as `Fig1_A1.png` / `Fig1_A2.png` rather than `FigureA1/A2.png`.

The tables above give the real mapping. Renaming the outputs would be a one-line change per
script, but it would also break any document that already `\includegraphics` these paths.

---

## How the panels are kept consistent

The applications share one display convention, implemented in `utility/heatmap_util.R` and
applied the same way by the main figure and its supplementary counterpart:

- **Cluster order** comes from `order_clusters_by_selection()` — the most parsimonious cluster
  first — so cluster C1 means the same thing in every panel of an application.
- **Cluster colours are keyed to display position**, not to the raw mixture-component index,
  so the panels agree even if the fit relabels its components.
- **Context order within a cluster** is decided once, in panel A, and reused by every later
  panel. Comparisons across methods align their labels with the Hungarian algorithm
  (`clue::solve_LSAP`) before anything is drawn.
- **A cluster in which no context could estimate a feature is left blank**, rather than shown
  as a selected zero (`mask_unobserved_clusters()`).

If you change an ordering rule, change it in `heatmap_util.R` — the figure scripts read it
from there, and editing one script alone will make its panels disagree with the others.
