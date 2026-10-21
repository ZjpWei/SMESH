# SMESH — analysis code

Code for the analyses in *Structure-learning meta-analysis resolves cross-context heterogeneity in microbiome associations*.

This folder contains the pipeline: data preprocessing, the simulation studies, and the
real-data model fits. Scripts that draw the figures and build the supplementary tables live
in [`../Rscript/`](../Rscript).

---

## Dependencies

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
| `PALM` (0.1.1) | per-study absolute-abundance association summaries | <https://github.com/ZjpWei/PALM_package> |
| `miMeta` (2.1.0) | Melody, the fixed-effect meta-analysis benchmark | <https://github.com/ZjpWei/miMeta> |
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
