# SMESH

Code and results for *Structure-learning microbiome meta-analysis reveals shared and
context-dependent microbial signatures*.

SMESH is a two-stage, precision-weighted finite-mixture meta-analysis for microbiome studies.
It takes per-study association summaries, learns how the studies group into a small number of
**clusters**, and simultaneously selects the taxa whose effect is shared by every cluster and
those that are specific to some of them — instead of forcing one pooled effect per taxon.

The repository holds three real-data applications and two simulation studies:

| Application | Contexts | Data | Fitted G |
|---|---|---|---|
| **Pan-disease** | 12 diseases, 32 datasets | GMrepo, 244 species | 2 |
| **Pan-tumor** | 15 tumour types, 36 datasets | HGMT, 159 genera | 4 |
| **Colorectal** | 26 study × stage strata, 2,679 samples | MetaPhlAn 4 profiles, 260 species | 2 |

---

## Where everything lives

```
MOSIAC/
├── Analysis/           the pipeline        — run these, in order
├── Rscript/            figures and tables  — run these afterwards
├── utility/            the estimator and shared helpers
├── Data/               raw inputs, exactly as downloaded
├── GMrepo_analysis/    pan-disease  — processed data, fits, panels
├── HGMT_analysis/      pan-tumor    — processed data, fits, panels
├── CRC_analysis/       colorectal   — processed data, fits, panels
├── Simulation/         simulation template, replicates and panels
└── Figure/             the assembled figures and supplementary tables
```

**Run everything from this directory.** Every path inside every script is written relative to
the repository root (`./Data/...`, `./utility/...`), and the scripts stop with an error if they
cannot see `utility/` and `Data/`.

```bash
cd /path/to/MOSIAC
Rscript Analysis/0_preprocessing.R
```

---

## Code

### `Analysis/` — the pipeline

Seven numbered scripts, from raw files to fitted models. Nothing here draws anything.

| Script | What it does |
|---|---|
| `0_preprocessing.R` | raw data → the inputs the association models consume, all three applications |
| `1_simulation_stage1.R` | one simulation replicate: SMESH against the benchmark methods |
| `2_simulation_stage2.R` | one simulation replicate: SMESH across four input summary types |
| `3_summary_statistics.R` | per-context association summaries, one method at a time |
| `4_real_data_analysis.R` | one model fit — one application, one method, one held-out context |
| `5_gmrepo_sensitivity.R` | pan-disease sensitivity: taxonomy level and feature filter |
| `6_comparison_models.R` | the comparison-method models at the selected G |

Only step 0 runs on its own; every other script takes command-line arguments. Steps 1–4 are
written to be run many times, one invocation per replicate or per fit — that is how they were
run on a cluster. See
[`Analysis/README.md`](Analysis/README.md) for the argument grids, the dependencies, and what
each step reads and writes.

### `Rscript/` — figures and tables

Fifteen scripts, one per figure or table. **These only draw**: they read results that already
exist and fit nothing, so `Analysis/` must have been run first. `Figure2`–`Figure6`,
`FigureS1`–`FigureS5`, `Make_TableS1`–`Make_TableS3` and `Make_TableS4_S6` map onto the numbers
in the paper; `Check_cluster_covariates.R` is a supporting check whose output is not shipped.
See [`Rscript/README.md`](Rscript/README.md).

### `utility/` — the estimator

Sourced by the scripts above; not an R package.

| File | Contents |
|---|---|
| `SMESH.R` | the estimator — `smesh.meta.summary()`, the two-stage fit and the G search |
| `PALM_tune.R` | tuning wrapper for the PALM per-study summaries |
| `ancombc.R` | ANCOM-BC2 wrapper used to produce comparison summary statistics |
| `heatmap_util.R` | everything the figures share: cluster palette, heat-map panels, display ordering, label matching |

The display conventions live in `heatmap_util.R` and nowhere else. Cluster order, cluster
colours, context order and the signature grouping are decided there, which is why the panels of
an application agree with each other and with the supplementary tables.

---

## Data

### `Data/` — raw inputs

Downloaded files, unmodified. Read only by `Analysis/0_preprocessing.R`.

| Folder | Contents |
|---|---|
| `Data/GMrepo/` | GMrepo species abundances, per-disease run lists, project metadata, taxon lookups |
| `Data/HGMT/` | per-project bacterial abundance tables and the study metadata workbook |
| `Data/CRC/` | MetaPhlAn 4 profiles, HUMAnN 3.6 pathways, per-study metadata and sequence lists |

### The three application folders

`GMrepo_analysis/`, `HGMT_analysis/` and `CRC_analysis/` have the same shape. **They hold data
and results only — no code.**

| Subfolder | Contents | Written by |
|---|---|---|
| `Data/` | processed abundances, metadata and covariates, plus `Summary_stat_<method>.Rdata` — the per-context association summaries every model consumes | steps 0 and 3 |
| `<APP>_loso/` | one file per (G, method, held-out context): `Model<G>_<method>_s<s>.Rdata`. **`s = 100` is the fit on all contexts**; `s = 1 … L` are the leave-one-context-out refits | step 4 |
| `<APP>_loso.tar.gz` | the same folder, archived — the form the fits are distributed in | — |
| `Model/` | the comparison models at the selected G: `SKMean_FE_G<G>.Rdata`, `SHC_FE_G<G>.Rdata`, `Melody_model.Rdata` | step 6 |
| `Figure/` | the PNG panels for this application, **one file per panel** | `Rscript/` |

`GMrepo_analysis/` additionally has `Sensitivity/`, the pan-disease genus-level and
stricter-filter refits behind panel D of Figure S3, written by step 5.

### `Simulation/`

| Item | Contents |
|---|---|
| `Data/CRC.Rdata` | the semi-synthetic template the replicates are generated from |
| `Sim_CRC_stage1/` | 991 replicates: SMESH against the benchmark methods |
| `Sim_CRC_stage2/` | 1,000 replicates: SMESH across four input summary types |
| `Sim_CRC_stage{1,2}.tar.gz` | the archived replicates |

`Figure2`, `Figure3`, `FigureS1` and `FigureS2` write their panels into `Simulation/Figure/`,
which the scripts create on first run.

A replicate is a single `Rdata` file named for its scenario and seed
(`res_Ka0.1_pos0.6_u0_scenario1_s1.Rdata`). Regenerating the full set takes a cluster; the
archives are there so the figures can be redrawn without it.

### `Figure/`

The final artefacts, assembled from the per-panel PNGs and the table scripts.

| File | Contents |
|---|---|
| `fig1.pdf` – `fig6.pdf`, `figS1.pdf` – `figS5.pdf` | the figures as they appear in the paper |
| `TableS1_pan_disease.{tex,csv}`, `TableS2_pan_tumor.{tex,csv}`, `TableS3_colorectal.{tex,csv}` | the study lists, as a `\input`-ready fragment and as CSV |
| `TableS4_S6_SMESH_effects.xlsx` | Supplementary Tables S4–S6: SMESH-PALM cluster-level effect estimates, one tab per application |
| `cluster_covariates.xlsx`, `cluster_association_tests.csv`, `context_characteristics.csv` | the cluster-covariate check — supporting output, not part of the paper |

Figure 1 is a schematic and is not generated from data.

---

## Running it end to end

```bash
cd /path/to/MOSIAC

Rscript Analysis/0_preprocessing.R                    # raw  -> processed
Rscript Analysis/3_summary_statistics.R GMrepo PALM   # per-context summaries, one method
Rscript Analysis/4_real_data_analysis.R HGMT 4 SMESH 100   # one fit; the real run is a grid
Rscript Analysis/6_comparison_models.R                # comparison models, all three applications

Rscript Rscript/Figure4.R                             # panels
Rscript Rscript/Make_TableS4_S6.R                     # supplementary tables
```

Steps 1 and 2 are the simulations and are independent of the real-data path. The two folder
READMEs — [`Analysis/README.md`](Analysis/README.md) and [`Rscript/README.md`](Rscript/README.md)
— give the argument grids, the exact inputs each script needs, and the runtimes.

## Contact

Zhoujingpeng Wei, Department of Biostatistics and Medical Informatics, University of
Wisconsin–Madison — <zwei74@wisc.edu>
