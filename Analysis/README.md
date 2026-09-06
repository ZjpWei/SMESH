# SMESH — analysis code

Code for the analyses in *Structure-learning microbiome meta-analysis reveals shared and
context-dependent microbial signatures*.

This folder contains the pipeline: data preprocessing, the simulation studies, and the
real-data model fits. Scripts that draw the figures and build the supplementary tables live
in [`../Rscript/`](../Rscript).

---

## Before you start

**Every script must be run from the repository root, not from this folder.** All paths
inside them are written relative to the root (`./Data/...`, `./utility/...`), and each script
stops with an error if it cannot see `utility/` and `Data/`.

```bash
cd /path/to/MOSIAC          # repository root
Rscript Analysis/0_preprocessing.R
```

### Dependencies

Most dependencies are on CRAN or Bioconductor:

```r
install.packages(c("tidyverse", "dplyr", "tidyr", "tibble", "purrr", "readr", "readxl",
                   "stringr", "ggplot2", "scales", "openxlsx", "jsonlite",
                   "abess", "glmtlp", "glmnet", "scam", "mclust", "cluster",
                   "sparcl", "sClust", "coin", "UpSetR"))

BiocManager::install(c("phyloseq", "ANCOMBC", "maaslin3", "MicrobiomeStat",
                       "S4Vectors", "TreeSummarizedExperiment"))
```

`MASS` and `Matrix` are also used; both ship with R.

Three are installed from source:

| Package | Purpose | Source |
|---|---|---|
| `PALM` | per-study absolute-abundance association summaries | <https://github.com/ZjpWei/PALM> |
| `miMeta` | Melody, the fixed-effect meta-analysis benchmark | <https://github.com/ZjpWei/miMeta> |
| `MIDASim` | semi-synthetic count simulation | <https://github.com/mengyu-he/MIDASim> |

`abess` and `glmtlp`, the sparse solvers behind SMESH, are on CRAN and are in the list above.

SMESH itself is not a package: it is sourced from [`../utility/SMESH.R`](../utility/SMESH.R),
together with `PALM_tune.R` and `ancombc.R` in the same folder.

---

## The pipeline

| Script | What it does | Output |
|---|---|---|
| `0_preprocessing.R` | raw data → the inputs the association models consume | `GMrepo_analysis/Data/`, `HGMT_analysis/Data/`, `CRC_analysis/Data/` |
| `1_simulation_stage1.R` | one simulation replicate: SMESH vs the benchmark methods | `Simulation/Sim_CRC_stage1/` |
| `2_simulation_stage2.R` | one simulation replicate: SMESH across four input summary types | `Simulation/Sim_CRC_stage2/` |
| `3_summary_statistics.R` | per-context association summaries, one method at a time | `<app>/Data/Summary_stat_<method>.Rdata` |
| `4_real_data_analysis.R` | one model fit for one application | `<app>_loso/Model<G>_<method>_s<s>.Rdata` |
| `5_gmrepo_sensitivity.R` | pan-disease sensitivity: taxonomy level, feature filter | `GMrepo_analysis/Sensitivity/res_genus_G2.Rdata`, `res_species_X4_G2.Rdata` |
| `6_comparison_models.R` | comparison-method models at the selected G, all three applications | `<app>/Model/SKMean_FE_G<G>.Rdata`, `SHC_FE_G<G>.Rdata`, `Melody_model.Rdata` |

### 0 — Preprocessing

```bash
Rscript Analysis/0_preprocessing.R
```

Rebuilds all three applications from raw files, in one pass. Every raw input lives under
`Data/`; every processed object is written into the application folder that consumes it. The
three parts are independent and each is wrapped in `local()`, so they can be read and edited
separately.

The GMrepo part begins with **step 0** (folded in from the old `Read.R`), which parses the raw
GMrepo download into `metadata_list.rds`, `abundance_list.rds` and `GMrepo_project_summary.xlsx`.
It is skipped when those three already exist, so a rerun does not repeat the slow parse of
`species_abundance.txt.gz`; delete them to force it.

| Application | Raw input | Result |
|---|---|---|
| Pan-disease (GMrepo) | `Data/GMrepo/` — `species_abundance.txt.gz`, `GMrepo_comparisons_summary.xlsx`, `GMrepo_runs_by_disease/`, lookups | 12 diseases, 32 comparisons, 5,784 samples, 244 species |
| Pan-tumor (HGMT) | `Data/HGMT/` (genus profiles, metadata, SRA run tables, Table S4) | 15 tumour types, 36 comparisons, 8,042 samples, 159 genera |
| Colorectal | `Data/CRC/` — `metadata/`, `metaphlan4_profiles/`, `sequence/` | 26 stage strata, 2,679 samples |

Nothing is overwritten: each part writes to its own output directory, leaving the objects the
paper was run on untouched.

### 1 and 2 — Simulations

Both take one replicate per invocation, so they parallelise across a cluster by looping over
the arguments. They are split because running every method in one job takes a long time.

```bash
# stage 1: SMESH vs True-cluster FE, SKM + FE, SHC + FE and Melody
Rscript Analysis/1_simulation_stage1.R <s> <Ka.per> <pos.pt> <u> <scenario>
Rscript Analysis/1_simulation_stage1.R 1 0.4 0.5 0 3

# stage 2: SMESH on summaries from four different single-study methods
Rscript Analysis/2_simulation_stage2.R <s> <Ka.per> <pos.pt> <u> <scenario> <method>
Rscript Analysis/2_simulation_stage2.R 1 0.4 0.5 0 3 PALM
```

| Argument | Meaning |
|---|---|
| `s` | replicate number, 1–100. Also the random seed (`set.seed(s + 2026)`) |
| `Ka.per` | proportion of features that carry signal: 0.1, 0.2 or 0.4 |
| `pos.pt` | probability an active signature is assigned to the case arm: 0.5, 0.75 or 1 |
| `u` | sequencing-depth unevenness |
| `scenario` | 1–5; scenario 1 has a single homogeneous group, 3 and 5 add subset-shared signatures |
| `method` | stage 2 only: `PALM`, `ANCOMBC2`, `LinDA` or `MaAsLin3` |

Studies are built with MIDASim from the colorectal strata. Each simulated study inherits the
**sample size, library-size distribution and taxon structure of a real stratum** — sample
sizes are not chosen by the simulation. The two arms are `floor(n/2)` each.

In stage 2 the branch on `method` happens *after* the data are generated, so a given
`(s, Ka.per, pos.pt, u, scenario)` yields the same simulated data for all four methods.

### 3 — Summary statistics

Runs one differential-abundance back-end per context and writes the files the
model fits consume. One invocation = one `(application, method)` cell:

```bash
Rscript Analysis/3_summary_statistics.R <application> <method>
Rscript Analysis/3_summary_statistics.R CRC PALM
Rscript Analysis/3_summary_statistics.R HGMT ANCOMBC2
```

| Argument | Values |
|---|---|
| `application` | `GMrepo`, `HGMT`, `CRC` |
| `method` | `PALM`, `ANCOMBC2`, `MaAsLin3`, `LinDA`, `Melody` |

For pan-disease and pan-tumor the per-dataset summaries are pooled within a
disease / tumour type by `PALM::palm.meta.summary()`; for colorectal the strata
are already the contexts, so no pooling happens. `correct = "tune"` throughout.

> **This overwrites the inputs every published result was fitted on.** The
> `Summary_stat_*.Rdata` files now in each `<app>/Data/` predate this script, and
> the retired `Prepare_*.R` files were not consistent about the PALM correction —
> `Prepare_GMrepo.R` never wrote `Summary_stat_SMESH.Rdata` at all. Back the
> current files up and compare before adopting new ones.

### 3 — Real-data model fits

One invocation fits one model to one application:

```bash
Rscript Analysis/4_real_data_analysis.R <application> <G> <method> <s>
Rscript Analysis/4_real_data_analysis.R HGMT 4 SMESH 100
```

| Argument | Values |
|---|---|
| `application` | `GMrepo`, `HGMT`, `CRC` |
| `G` | number of clusters to fit |
| `method` | `SMESH`, `ANCOMBC2`, `MaAsLin3`, `LinDA`, `SKM`, `SHC` |
| `s` | context to hold out — **or `100` to fit every context** |

> **The `s = 100` convention.** The code subsets with `summary.stats[-s]`, and R ignores a
> negative index past the end of a list, so `s = 100` drops nothing and fits all contexts.
> `s = 1 … L` leaves out context *s*, which is how the leave-one-context-out consensus and
> stability results are produced. The saved file names follow the same convention:
> `Model4_SMESH_s100.Rdata` is the full-data fit.

The full result set is produced by looping over `G`, `method` and `s`, for example:

```bash
for G in 2 3 4 5 6; do
  for M in SMESH ANCOMBC2 MaAsLin3 LinDA SKM SHC; do
    for s in $(seq 1 15) 100; do
      Rscript Analysis/4_real_data_analysis.R HGMT $G $M $s
    done
  done
done
```

Settings are not uniform across applications — permutation counts and the names of the saved
summary-statistic objects differ, and they are held in the `APPS` table at the top of the
script exactly as each application was run.

### 5 — Comparison-method models

Run **after** `4_real_data_analysis.R`, once G has been chosen for each
application. It produces the three comparison-method objects the figure scripts
load from `<app>/Model/`:

```bash
Rscript Analysis/6_comparison_models.R          # all three applications
Rscript Analysis/6_comparison_models.R HGMT     # just one
Rscript Analysis/6_comparison_models.R HGMT 3   # override the G to use
```

| Method | What it does | Output object |
|---|---|---|
| SKM + FE | sparse k-means clustering, then a fixed-effect PALM fit per cluster | `detect.signal` |
| SHC + FE | sparse hierarchical clustering, same second step | `detect.signal` |
| Melody | fixed-effect meta-analysis benchmark, no clustering and no G | `Melody_mod` |

SKM and SHC are two-step methods. `4_real_data_analysis.R` already does the
clustering half and saves it as `tab` in `<app>_loso/Model<G>_{SKM,SHC}_s100.Rdata`;
this script does the within-cluster fixed-effect half, which is what the figures
plot. If a clustering is missing it says which command to run rather than failing.

The selected G per application lives in the `APPS` table at the top of the script
(GMrepo 2, HGMT 4, CRC 2). Change it there if the selected model changes.

> **Melody needs its own summary statistics.** `Summary_stat_Melody.Rdata` is
> currently present for GMrepo only, so Melody is regenerated there and skipped
> with an explanatory message for CRC and HGMT — whose existing
> `Model/Melody_model.Rdata` is left untouched. Producing those summaries needs
> `miMeta::melody.get.summary()` on the per-sample data, part of the same
> summary-statistic step noted under Reproducibility below.

### 4 — Pan-disease sensitivity analyses

Two checks on the pan-disease result, each re-running the whole pipeline with one
thing changed and everything else held fixed:

```bash
Rscript Analysis/5_gmrepo_sensitivity.R taxonomy    # genus level instead of species
Rscript Analysis/5_gmrepo_sensitivity.R filtering   # stringent feature filter, X = 4
```

| Argument | What varies | Output |
|---|---|---|
| `taxonomy` | counts aggregated to genus; the dense feature set is the species-level set mapped through the NCBI lineage | `GMrepo_analysis/Sensitivity/res_genus_G2.Rdata` (plus `res_genus.Rdata`) |
| `filtering` | species level, but a species must be estimable in more than `X = 4` of the 12 diseases | `GMrepo_analysis/Sensitivity/res_species_X4_G2.Rdata` |

Both fit all six methods at a fixed `G = 2` and feed panel D of Supplementary
Figure S3, which reports the ARI between each sensitivity clustering and the main
one. The differential-abundance stage takes roughly 20 minutes, so each run is
checkpointed twice and resumes automatically; delete the `*_summary_stats.Rdata`
or `*_smesh_fits.Rdata` file to force a stage to recompute.

**Re-running overwrites the results the published figure uses.** The outputs are
already in `GMrepo_analysis/Sensitivity/`; only re-run if you intend to regenerate them.

---

## Producing the figures and tables

Once the fits exist, everything in the paper is drawn by the scripts in
[`../Rscript/`](../Rscript), also run from the repository root:

```bash
Rscript Rscript/Figure4.R          # pan-disease
Rscript Rscript/FigureS3.R         # pan-disease sensitivity analyses
Rscript Rscript/Make_TableS1.R     # study list, pan-disease
```

`Figure2`–`Figure6` and `FigureS1`–`FigureS5` correspond to the figure numbers in the paper;
`Make_TableS1`–`Make_TableS3` and `Make_TableS4_S6` build the supplementary tables. Figure 1 is a schematic and is
not generated from data.

---

## Reproducibility notes

**Preprocessing.** The GMrepo and HGMT parts reproduce the analysed inputs exactly. The
colorectal part reproduces the design and every stratum count (26 strata, 2,679 samples), but
not the individual control assignment: controls are allocated across a study's stage strata in
proportion to its case counts, and the published allocation was drawn at random with a seed
that was not recorded. `partition_seed` in the script fixes it from here on.

**Simulation template.** `Simulation/Data/CRC.Rdata` is an input the simulations read and
nothing in the repository regenerates: it is the colorectal strata object as it stood before
the current preprocessing (26 strata, 3,931 samples), and the published simulations were run
from it. `CRC_analysis/Data/CRC_strata.Rdata` is the current equivalent (2,679 samples) and is
**not** interchangeable with it.

**Summary statistics.** `4_real_data_analysis.R` reads the `Summary_stat_*.Rdata` files in
each application's `Data/` folder. Those files are inputs here, not outputs: **no script in
this repository regenerates them.** They were produced by the per-application `Prepare_*.R`
scripts, which have been retired to `Backup/prepare/` now that `0_preprocessing.R` covers the
preprocessing. That step runs PALM, ANCOMBC2, MaAsLin3 and LinDA per context; folding it into
this folder is the one remaining gap in the pipeline.

**SMESH implementation.** All scripts source `utility/SMESH.R`, the single maintained
implementation. Earlier generations of the estimator (including `SMESH_v4.R`) were retired to
`Backup/utility/`. The saved results were fitted with those earlier versions, so re-running
will not reproduce them bit for bit.

---

## Repository layout

```
Analysis/            this folder — preprocessing, simulation, real-data analysis
Rscript/             figure and table scripts
utility/             SMESH, PALM tuning and helper functions shared by both
Data/                raw input
GMrepo_analysis/     processed data and results, pan-disease
HGMT_analysis/       processed data and results, pan-tumor
CRC_analysis/        processed data and results, colorectal
Simulation/          simulation results
```
