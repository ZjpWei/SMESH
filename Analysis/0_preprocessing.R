# ============================================================================
#   0_preprocessing.R  --  build the analysis inputs for all three applications
# ============================================================================
#  Self-contained: raw files under ./Data go in, processed objects come out
#  inputs the association models consume come out.  The three applications are
#  independent; each is wrapped in local() so their helper functions and loop
#  variables cannot collide.
#
#  Run from the project root:  Rscript Analysis/0_preprocessing.R
#
#  Nothing here overwrites the objects the paper was run on; each part writes to
#  its own output directory.  GMrepo and HGMT reproduce the analysed inputs
#  exactly; CRC reproduces the design and every stratum count, but not the
#  individual control assignment, which was drawn at random with a seed that was
#  not recorded.
# ============================================================================

  rm(list = ls())

  library(dplyr)
  library(openxlsx)
  library(jsonlite)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)

  if (!dir.exists("utility") || !dir.exists("Data"))
    stop("run this from the project root, e.g.  Rscript Analysis/0_preprocessing.R")


# ============================================================================
#   PART 1 - PAN-DISEASE (GMrepo)
# ============================================================================
#   Pan-disease application (GMrepo) - data preprocessing
#  Raw input (Data/GMrepo/):
#    species_abundance.txt.gz          GMrepo species-level relative abundance
#    GMrepo_comparisons_summary.xlsx   phenotype comparisons index
#    GMrepo_runs_by_disease/*.tsv      per-disease run tables, incl. the Health pool
#    taxon_id_species_lookup.csv       cached NCBI taxon id -> species name
#
#  Step 0 (from Read.R) turns that download into the per project-disease objects
#  the rest of this part consumes, written alongside the other processed output:
#    metadata_list.rds             per project-disease sample metadata
#    abundance_list.rds            per project-disease relative abundances
#    GMrepo_project_summary.xlsx   assay type, group sizes, candidate confounders
#
#  Output (GMrepo_analysis/Data/):
#    rel.abd.rds             per study count matrix, the PALM input
#    covariate.interest.rds  disease indicator per sample
#    covariate.adjust.rds    covariates to adjust for, per study
#    cluster.rds             repeated-measure grouping where present
#    plus the intermediate objects and the confounder reports
#
#  Steps
#    1. keep metagenomic project-disease pairs among the 12 target diseases
#    2. require >= 10 cases and >= 10 controls, and >= 100 prevalent species
#    3. name species from their NCBI taxon id, recover sequencing depth
#    4. test each candidate confounder case vs control, impute partial gaps
#    5. rescale relative abundance to depth-weighted pseudo-counts
#    6. assemble the PALM inputs
#
#
#  The species name lookup is read from the cached CSV; it only contacts NCBI if
#  a taxon id is missing from that cache.


local({
## selected Disease_HGMT
Dis_lst <- c("Arthritis, Rheumatoid",
             "Colitis, Ulcerative",
             "Crohn Disease",
             "Diabetes Mellitus, Type 1",
             "Diabetes Mellitus, Type 2",
             "Fibromyalgia",
             "Kidney Failure, Chronic",
             "Obesity",
             "Psoriasis",
             "Fatigue Syndrome, Chronic",
             "Irritable Bowel Syndrome",
             "Non-alcoholic Fatty Liver Disease")

## `in_dir` is the raw GMrepo download; everything produced here - including
## step 0's parsed intermediate - goes to the application folder.
in_dir  <- "./Data/GMrepo"
out_dir <- "./GMrepo_analysis/Data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================================
## Step 0 (from Read.R): raw GMrepo download -> per project-disease objects
## =========================================================================
## Skipped when its three outputs already exist; delete them to force a rerun.

local({
  in_dir_gm  <- "./Data/GMrepo"
  out_dir_gm <- "./GMrepo_analysis/Data"
  dir.create(out_dir_gm, showWarnings = FALSE, recursive = TRUE)

  produced <- file.path(out_dir_gm, c("metadata_list.rds", "abundance_list.rds",
                                      "GMrepo_project_summary.xlsx"))
  if (all(file.exists(produced))) {
    message("step 0: outputs present, skipping the raw parse")
    return(invisible(NULL))
  }
  message("step 0: parsing the raw GMrepo download")

  health_mesh   <- "D006262"                                     # "Health" phenotype
  health_file   <- file.path(runs_dir, paste0("all_runs_associated_with_", health_mesh, ".tsv"))
  prevalence_thresh <- 0.10

  dis_lst <- c("Arthritis, Rheumatoid",
               "Colitis, Ulcerative",
               "Crohn Disease",
               "Diabetes Mellitus, Type 1",
               "Diabetes Mellitus, Type 2",
               "Fibromyalgia",
               "Kidney Failure, Chronic",
               "Obesity",
               "Psoriasis",
               "Fatigue Syndrome, Chronic",
               "Irritable Bowel Syndrome",
               "Non-alcoholic Fatty Liver Disease")

  ## metadata columns worth flagging as confounders if not entirely NA/empty
  confounder_cols <- c(host_age                    = "Age",
                       sex                         = "Gender",
                       BMI                         = "BMI",
                       country                     = "Country",
                       diet                        = "Diet",
                       disease_stage               = "Disease_Stage",
                       Recent_Antibiotics_Use      = "Recent_Antibiotics_Use",
                       antibiotics_used            = "Antibiotics_Used",
                       Antibiotics_Dose            = "Antibiotics_Dose",
                       Days_Without_Antibiotics_Use = "Days_Without_Antibiotics_Use")

  ## ---------------------------------------------------------------------
  ## 1. Find the project IDs related to the diseases in the list
  ## ---------------------------------------------------------------------

  GMrepo_comparisons_summary <- read_excel(file.path(in_dir_gm, "GMrepo_comparisons_summary.xlsx"))

  file_sub <- GMrepo_comparisons_summary %>%
    dplyr::filter(phenotype2_term %in% dis_lst) %>%
    dplyr::distinct(phenotype2, phenotype2_term)

  MESH_lst <- file_sub$phenotype2

  meta_file_dir <- file.path(runs_dir,
                              paste0(gsub("[^A-Za-z0-9]+", "_", file_sub$phenotype2_term),
                                     "_",
                                     file_sub$phenotype2,
                                     ".tsv"))

  ## ---------------------------------------------------------------------
  ## All species-level relative abundance (one row per run/species)
  ## ---------------------------------------------------------------------

  species_abundance <- fread(cmd = paste("gunzip -c", shQuote(file.path(in_dir_gm, "species_abundance.txt.gz"))))
  species_abundance <- species_abundance[taxon_rank_level == "species" & ncbi_taxon_id != -1]
  species_abundance[, taxon_col := paste0("taxon_", ncbi_taxon_id)]
  setkey(species_abundance, accession_id)

  ## Shared pool of healthy/control runs, subset per project below
  health_meta_all <- fread(health_file, colClasses = "character", na.strings = "")

  ## ---------------------------------------------------------------------
  ## 2-5. For each disease -> each project: metadata, abundance matrix, summary stats
  ## ---------------------------------------------------------------------

  metadata_list  <- list()
  abundance_list <- list()
  summary_rows   <- list()

  for (i in seq_len(nrow(file_sub))) {

    mesh <- file_sub$phenotype2[i]
    term <- file_sub$phenotype2_term[i]
    meta_path <- meta_file_dir[i]

    if (!file.exists(meta_path)) {
      warning("Metadata file not found, skipping: ", meta_path)
      next
    }

    case_meta <- fread(meta_path, colClasses = "character", na.strings = "")
    case_meta[, sample_group := "disease"]

    projects <- unique(case_meta$project_id)

    for (proj in projects) {

      proj_case <- case_meta[project_id == proj]
      proj_ctrl <- health_meta_all[project_id == proj]

      if (nrow(proj_ctrl) > 0) {
        proj_ctrl[, sample_group := "health"]
        proj_meta <- rbind(proj_case, proj_ctrl, fill = TRUE)
      } else {
        proj_meta <- proj_case
      }

      key <- paste0(proj, "_", mesh)

      ## ---- extract corresponding species-level abundance for this project's runs ----
      run_ids  <- unique(proj_meta$run_id)
      abund_sub <- species_abundance[.(run_ids), on = "accession_id", nomatch = 0]

      ## data.frame with 0 rows (not NULL) so this key stays present in
      ## abundance_list, aligned 1:1 with metadata_list, even with no data
      mat <- data.frame(accession_id = character(0))
      n_species      <- 0
      n_species_prev <- 0
      n_disease      <- 0
      n_health       <- 0

      if (nrow(abund_sub) > 0) {
        wide <- tidyr::pivot_wider(abund_sub,
                                    id_cols     = accession_id,
                                    names_from  = taxon_col,
                                    values_from = relative_abundance,
                                    values_fill = 0,
                                    values_fn   = sum)
        mat <- as.data.frame(wide)

        n_species      <- ncol(mat) - 1
        prevalence     <- colMeans(mat[, -1, drop = FALSE] > 0)
        n_species_prev <- sum(prevalence >= prevalence_thresh)

        ## Not every run in the metadata has a QC-passed abundance profile
        ## (species_abundance.txt.gz only covers QCstatus == 1 runs, and not
        ## even all of those - e.g. amplicon runs with no species-level calls).
        ## Align metadata rows to exactly the samples present in `mat`, same order,
        ## so metadata_list and abundance_list always have matching sample sets.
        proj_meta <- proj_meta[match(mat$accession_id, proj_meta$run_id)]

        n_disease <- sum(proj_meta$sample_group == "disease")
        n_health  <- sum(proj_meta$sample_group == "health")
      } else {
        proj_meta <- proj_meta[0]
      }

      ## ---- non-empty confounders present in this project's metadata ----
      avail_conf <- character(0)
      for (col in names(confounder_cols)) {
        if (col %in% names(proj_meta)) {
          vals <- proj_meta[[col]]
          if (any(!is.na(vals) & vals != "")) {
            avail_conf <- c(avail_conf, confounder_cols[[col]])
          }
        }
      }

      ## ---- sequencing/data type(s) present (Amplicon, Metagenomics, ...) ----
      exp_types <- proj_meta$experiment_type
      exp_types <- exp_types[!is.na(exp_types) & exp_types != ""]
      data_type <- paste(sort(unique(exp_types)), collapse = "; ")

      metadata_list[[key]]  <- proj_meta
      abundance_list[[key]] <- mat

      summary_rows[[key]] <- data.frame(
        project_id                  = proj,
        disease_name                = term,
        disease_mesh                = mesh,
        data_type                   = data_type,
        n_species                   = n_species,
        n_species_prevalence_filtered = n_species_prev,
        n_disease_samples           = n_disease,
        n_health_samples            = n_health,
        confounders                 = paste(avail_conf, collapse = "; "),
        stringsAsFactors            = FALSE
      )
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)

  ## ---------------------------------------------------------------------
  ## Save outputs
  ## ---------------------------------------------------------------------

  saveRDS(metadata_list,  file.path(out_dir_gm, "metadata_list.rds"))
  saveRDS(abundance_list, file.path(out_dir_gm, "abundance_list.rds"))
  openxlsx::write.xlsx(summary_df, file.path(out_dir_gm, "GMrepo_project_summary.xlsx"))

})

## =========================================================================
## Part 1 (from Analysis.R): Metagenomics-only filter, sample-size and
## species-prevalence filters, species naming, sequencing depth recovery
## =========================================================================

## ---------------------------------------------------------------------
## Load the objects step 0 just wrote
## ---------------------------------------------------------------------

metadata_list  <- readRDS(file.path(out_dir, "metadata_list.rds"))
abundance_list <- readRDS(file.path(out_dir, "abundance_list.rds"))
project_summary <- read.xlsx(file.path(out_dir, "GMrepo_project_summary.xlsx"))

## ---------------------------------------------------------------------
## Keep only project-disease pairs with Metagenomics (species-level) data
## ---------------------------------------------------------------------

keep_summary <- project_summary %>%
  dplyr::filter(disease_name %in% Dis_lst,
                data_type == "Metagenomics",
                n_species > 0)

keep_keys <- paste0(keep_summary$project_id, "_", keep_summary$disease_mesh)

metadata_list_mgx  <- metadata_list[keep_keys]
abundance_list_mgx <- abundance_list[keep_keys]

## ---------------------------------------------------------------------
## Save filtered lists
## ---------------------------------------------------------------------

saveRDS(metadata_list_mgx,  file.path(out_dir, "metadata_list_metagenomics.rds"))
saveRDS(abundance_list_mgx, file.path(out_dir, "abundance_list_metagenomics.rds"))

## ---------------------------------------------------------------------
## Species name lookup (abundance columns are "taxon_<ncbi_taxon_id>";
## GMrepo does not ship a name table, so resolve names via NCBI Entrez)
## ---------------------------------------------------------------------

taxon_ids <- unique(unlist(lapply(abundance_list_mgx, function(x) {
  sub("^taxon_", "", setdiff(colnames(x), "accession_id"))
})))

fetch_taxon_names <- function(taxon_ids, batch_size = 200, delay = 0.4) {
  chunks <- split(taxon_ids, ceiling(seq_along(taxon_ids) / batch_size))
  do.call(rbind, lapply(chunks, function(ids) {
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                  "?db=taxonomy&retmode=json&id=", paste(ids, collapse = ","))
    res  <- jsonlite::fromJSON(url)
    recs <- res$result[as.character(ids)]
    Sys.sleep(delay)
    data.frame(
      ncbi_taxon_id   = as.integer(ids),
      scientific_name = vapply(recs, function(r) r$scientificname %||% NA_character_, character(1)),
      rank            = vapply(recs, function(r) r$rank %||% NA_character_, character(1)),
      stringsAsFactors = FALSE
    )
  }))
}

lookup_path     <- file.path(in_dir,  "taxon_id_species_lookup.csv")
lookup_out_path <- file.path(out_dir, "taxon_id_species_lookup.csv")

if (file.exists(lookup_path)) {
  taxon_lookup <- read.csv(lookup_path, stringsAsFactors = FALSE)
  missing_ids  <- setdiff(taxon_ids, taxon_lookup$ncbi_taxon_id)
  if (length(missing_ids) > 0) {
    taxon_lookup <- rbind(taxon_lookup, fetch_taxon_names(missing_ids))
  }
} else {
  taxon_lookup <- fetch_taxon_names(taxon_ids)
}

write.csv(taxon_lookup, lookup_out_path, row.names = FALSE)

## ---------------------------------------------------------------------
## 1. Drop project-disease pairs without >= 10 case AND >= 10 control samples
## ---------------------------------------------------------------------

min_group_n <- 10

keep_summary2 <- keep_summary %>%
  dplyr::filter(n_disease_samples >= min_group_n, n_health_samples >= min_group_n)

keep_keys2 <- paste0(keep_summary2$project_id, "_", keep_summary2$disease_mesh)

metadata_list_f  <- metadata_list_mgx[keep_keys2]
abundance_list_f <- abundance_list_mgx[keep_keys2]

## ---------------------------------------------------------------------
## 2. Drop species with zero prevalence, then drop project-disease pairs
##    left with fewer than 100 species
## ---------------------------------------------------------------------

min_species_n <- 100

abundance_list_f2 <- list()
for (key in names(abundance_list_f)) {
  mat <- abundance_list_f[[key]]
  species_cols <- setdiff(colnames(mat), "accession_id")
  prevalent <- species_cols[colSums(mat[, species_cols, drop = FALSE] > 0) > 0]

  if (length(prevalent) >= min_species_n) {
    abundance_list_f2[[key]] <- mat[, c("accession_id", prevalent), drop = FALSE]
  }
}

## ---------------------------------------------------------------------
## 3. Keep only metadata entries matching the filtered abundance data
## ---------------------------------------------------------------------

final_keys <- names(abundance_list_f2)

metadata_list_final  <- metadata_list_f[final_keys]
abundance_list_final <- abundance_list_f2[final_keys]

## ---------------------------------------------------------------------
## 3b. Recover sequencing depth (nr_reads_sequenced, read as character in
##     Read.R) as numeric, on the final filtered/matched samples only
## ---------------------------------------------------------------------

for (key in names(metadata_list_final)) {
  metadata_list_final[[key]]$nr_reads_sequenced <-
    as.numeric(metadata_list_final[[key]]$nr_reads_sequenced)
}

## ---------------------------------------------------------------------
## 4. Attach species names to the abundance data column names
##    ("Scientific name [taxon_id]"; falls back to "Unclassified [taxon_id]")
## ---------------------------------------------------------------------

label_species <- function(taxon_ids) {
  idx  <- match(as.integer(taxon_ids), taxon_lookup$ncbi_taxon_id)
  name <- taxon_lookup$scientific_name[idx]
  ifelse(is.na(name),
         paste0("Unclassified [", taxon_ids, "]"),
         paste0(name, " [", taxon_ids, "]"))
}

for (key in names(abundance_list_final)) {
  mat <- abundance_list_final[[key]]
  species_cols <- setdiff(colnames(mat), "accession_id")
  taxon_ids_col <- sub("^taxon_", "", species_cols)
  colnames(mat)[colnames(mat) %in% species_cols] <- label_species(taxon_ids_col)
  abundance_list_final[[key]] <- mat
}

## ---------------------------------------------------------------------
## 5. Move accession_id to rownames so the abundance data is numeric-only
## ---------------------------------------------------------------------

for (key in names(abundance_list_final)) {
  mat <- abundance_list_final[[key]]
  rownames(mat) <- mat$accession_id
  abundance_list_final[[key]] <- as.matrix(mat[, setdiff(colnames(mat), "accession_id"), drop = FALSE])
}

## ---------------------------------------------------------------------
## Save final filtered, species-named lists
## ---------------------------------------------------------------------

saveRDS(metadata_list_final,  file.path(out_dir, "metadata_list_final.rds"))
saveRDS(abundance_list_final, file.path(out_dir, "abundance_list_final.rds"))

## =========================================================================
## Part 2 (from Analysis2.R): confounder detection, testing, imputation,
## and depth-scaled pseudo-counts
## =========================================================================

## Minimum non-missing samples required *in each group* (case, control) before
## a confounder is even tested. Below this, a test result would be driven by
## 1-2 observations and isn't trustworthy - such confounders are reported as
## "not tested" rather than forced through a test.
min_group_n <- 5

## ---------------------------------------------------------------------
## 1. Which confounders are present in the metadata
## ---------------------------------------------------------------------

confounder_candidates <- c("host_age", "sex", "BMI", "country", "diet", "disease_stage",
                            "Recent_Antibiotics_Use", "antibiotics_used", "Antibiotics_Dose",
                            "Days_Without_Antibiotics_Use")

all_meta <- dplyr::bind_rows(metadata_list_final, .id = "project_disease")

is_present <- function(x) any(!is.na(x) & x != "")

confounders_present <- confounder_candidates[
  sapply(confounder_candidates, function(col) col %in% names(all_meta) && is_present(all_meta[[col]]))
]

cat("Confounders present anywhere in the final data:\n")
cat(paste(" -", confounders_present), sep = "\n")

## per-project non-missing counts, for reference
confounder_availability <- t(sapply(metadata_list_final, function(md) {
  sapply(confounders_present, function(col) sum(!is.na(md[[col]]) & md[[col]] != ""))
}))
confounder_availability <- data.frame(project_disease = rownames(confounder_availability),
                                       confounder_availability, row.names = NULL,
                                       check.names = FALSE)

write.csv(confounder_availability, file.path(out_dir, "confounder_availability_by_project.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------
## 2. Case vs. control test per confounder, per project-disease pair
##    - continuous  -> Wilcoxon rank-sum test
##    - binary      -> Fisher's exact test (2x2)
##    - categorical -> Fisher's exact test with simulated p-value (>2 levels)
##    Missing values are handled by complete-case exclusion per confounder
##    (see write-up for why), and a confounder is only tested if both groups
##    retain >= min_group_n non-missing observations.
## ---------------------------------------------------------------------

classify_confounder <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("empty")
  n_levels <- length(unique(x))
  if (n_levels < 2) return("no_variance")
  numeric_x <- suppressWarnings(as.numeric(x))
  if (!any(is.na(numeric_x)) && n_levels > 2) return("continuous")
  if (n_levels == 2) return("binary")
  return("categorical")
}

test_confounder <- function(values, group, min_group_n) {
  keep <- !is.na(values) & values != "" & !is.na(group)
  values <- values[keep]
  group  <- group[keep]

  type <- classify_confounder(values)
  if (type %in% c("empty", "no_variance")) {
    return(list(p_value = NA_real_, type = type, note = "no usable variance"))
  }

  group_n <- table(group)
  if (length(group_n) < 2 || any(group_n < min_group_n)) {
    return(list(p_value = NA_real_, type = type, note = "insufficient non-missing samples in one or both groups"))
  }

  p_value <- switch(type,
    continuous = tryCatch(wilcox.test(as.numeric(values) ~ group, exact = FALSE)$p.value, error = function(e) NA_real_),
    binary     = tryCatch(fisher.test(table(group, values))$p.value, error = function(e) NA_real_),
    categorical = tryCatch(fisher.test(table(group, values), simulate.p.value = TRUE, B = 10000)$p.value,
                            error = function(e) NA_real_)
  )

  list(p_value = p_value, type = type, note = "tested")
}

confounder_test_results <- setNames(
  lapply(names(metadata_list_final), function(key) {
    md    <- metadata_list_final[[key]]
    group <- md$sample_group
    setNames(
      lapply(confounders_present, function(col) test_confounder(md[[col]], group, min_group_n)),
      confounders_present
    )
  }),
  names(metadata_list_final)
)

saveRDS(confounder_test_results, file.path(out_dir, "confounder_test_results.rds"))

## ---------------------------------------------------------------------
## 3. Per project-disease list of significant confounders (p <= 0.05);
##    NULL when none are significant (or none were testable)
## ---------------------------------------------------------------------

sig_threshold <- 0.05

significant_confounders <- setNames(
  lapply(names(metadata_list_final), function(key) {
    res    <- confounder_test_results[[key]]
    pvals  <- sapply(res, function(r) r$p_value)
    sig    <- names(pvals)[!is.na(pvals) & pvals <= sig_threshold]
    if (length(sig) == 0) NULL else sig
  }),
  names(metadata_list_final)
)

saveRDS(significant_confounders, file.path(out_dir, "significant_confounders.rds"))

## ---------------------------------------------------------------------
## 4. Confounder test report (xlsx): one row per project-disease x confounder,
##    with its p-value, whether it's significant, and non-missing proportion
## ---------------------------------------------------------------------

disease_lookup <- unique(project_summary[, c("project_id", "disease_mesh", "disease_name")])

confounder_report <- dplyr::bind_rows(lapply(names(metadata_list_final), function(key) {
  md      <- metadata_list_final[[key]]
  res     <- confounder_test_results[[key]]
  n_total <- nrow(md)

  parts        <- strsplit(key, "_")[[1]]
  proj         <- parts[1]
  mesh         <- parts[2]
  disease_name <- disease_lookup$disease_name[disease_lookup$project_id == proj & disease_lookup$disease_mesh == mesh]
  if (length(disease_name) == 0) disease_name <- NA_character_

  dplyr::bind_rows(lapply(names(res), function(col) {
    r         <- res[[col]]
    non_na_n  <- sum(!is.na(md[[col]]) & md[[col]] != "")
    data.frame(
      project_id        = proj,
      disease_name       = disease_name,
      disease_mesh       = mesh,
      confounder         = col,
      p_value            = r$p_value,
      significant        = !is.na(r$p_value) & r$p_value <= sig_threshold,
      non_na_proportion  = round(non_na_n / n_total, 3),
      stringsAsFactors   = FALSE
    )
  }))
}))

write.xlsx(confounder_report, file.path(out_dir, "confounder_test_report.xlsx"))

## ---------------------------------------------------------------------
## 5. Impute partially-missing confounders (median for continuous, mode for
##    binary/categorical) with a companion "<col>_imputed" missingness flag,
##    so adjustment models don't have to drop samples over one covariate.
##    Confounders that are entirely missing, or have no missing values, are
##    left untouched (nothing to learn an imputed value from / nothing to fix).
## ---------------------------------------------------------------------

impute_confounder <- function(x) {
  missing <- is.na(x) | x == ""

  if (length(x[!missing]) == 0 || !any(missing)) {
    return(list(imputed = x, missing_flag = as.integer(missing), method = NA_character_))
  }

  x_nonmiss <- x[!missing]
  type <- classify_confounder(x)

  if (type == "continuous") {
    fill   <- median(as.numeric(x_nonmiss))
    method <- "median"
  } else {
    fill   <- names(sort(table(x_nonmiss), decreasing = TRUE))[1]
    method <- "mode"
  }

  x[missing] <- as.character(fill)
  list(imputed = x, missing_flag = as.integer(missing), method = method)
}

imputation_info <- setNames(
  lapply(names(metadata_list_final), function(key) {
    md <- metadata_list_final[[key]]
    setNames(lapply(confounders_present, function(col) impute_confounder(md[[col]])), confounders_present)
  }),
  names(metadata_list_final)
)

metadata_list_imputed <- setNames(
  lapply(names(metadata_list_final), function(key) {
    md   <- metadata_list_final[[key]]
    info <- imputation_info[[key]]
    for (col in confounders_present) {
      md[[col]] <- info[[col]]$imputed
      partially_missing <- any(info[[col]]$missing_flag == 1) && !all(info[[col]]$missing_flag == 1)
      if (partially_missing) {
        md[[paste0(col, "_imputed")]] <- info[[col]]$missing_flag
      }
    }
    md
  }),
  names(metadata_list_final)
)

saveRDS(metadata_list_imputed, file.path(out_dir, "metadata_list_imputed.rds"))

## ---------------------------------------------------------------------
## 6. Report confounders after imputation: re-run the case-vs-control test
##    on the imputed data (now at full N), alongside the original
##    (pre-imputation) result and whether missingness itself differed by
##    case/control group (a cheap check for missing-not-at-random bias).
## ---------------------------------------------------------------------

test_missingness_association <- function(missing_flag, group) {
  keep <- !is.na(group)
  missing_flag <- missing_flag[keep]
  group <- group[keep]
  if (length(unique(missing_flag)) < 2 || length(unique(group)) < 2) return(NA_real_)
  tryCatch(fisher.test(table(group, missing_flag))$p.value, error = function(e) NA_real_)
}

confounder_test_results_imputed <- setNames(
  lapply(names(metadata_list_imputed), function(key) {
    md    <- metadata_list_imputed[[key]]
    group <- md$sample_group
    setNames(
      lapply(confounders_present, function(col) test_confounder(md[[col]], group, min_group_n)),
      confounders_present
    )
  }),
  names(metadata_list_imputed)
)

saveRDS(confounder_test_results_imputed, file.path(out_dir, "confounder_test_results_imputed.rds"))

confounder_report_imputed <- dplyr::bind_rows(lapply(names(metadata_list_final), function(key) {
  md_before <- metadata_list_final[[key]]
  group     <- md_before$sample_group
  res_before <- confounder_test_results[[key]]
  res_after  <- confounder_test_results_imputed[[key]]
  info       <- imputation_info[[key]]
  n_total    <- nrow(md_before)

  parts        <- strsplit(key, "_")[[1]]
  proj         <- parts[1]
  mesh         <- parts[2]
  disease_name <- disease_lookup$disease_name[disease_lookup$project_id == proj & disease_lookup$disease_mesh == mesh]
  if (length(disease_name) == 0) disease_name <- NA_character_

  dplyr::bind_rows(lapply(confounders_present, function(col) {
    n_missing <- sum(info[[col]]$missing_flag)
    data.frame(
      project_id                     = proj,
      disease_name                   = disease_name,
      disease_mesh                   = mesh,
      confounder                     = col,
      imputation_method              = info[[col]]$method,
      n_imputed                      = n_missing,
      non_na_proportion_before       = round((n_total - n_missing) / n_total, 3),
      missingness_assoc_with_group_p = test_missingness_association(info[[col]]$missing_flag, group),
      p_value_before                 = res_before[[col]]$p_value,
      significant_before             = !is.na(res_before[[col]]$p_value) & res_before[[col]]$p_value <= sig_threshold,
      p_value_after                  = res_after[[col]]$p_value,
      significant_after              = !is.na(res_after[[col]]$p_value) & res_after[[col]]$p_value <= sig_threshold,
      stringsAsFactors                = FALSE
    )
  }))
}))

write.xlsx(confounder_report_imputed, file.path(out_dir, "confounder_report_after_imputation.xlsx"))

## ---------------------------------------------------------------------
## 7. Depth-scaled pseudo-counts: relative_abundance(%)/100 * nr_reads_sequenced
##
##    Caveat (see write-up): this is an approximation, not a reconstruction
##    of true per-species read counts. MetaPhlAn's relative abundance is
##    derived from marker-gene coverage, not a raw read tally, so there is
##    no exact linear map back to counts; and nr_reads_sequenced is the RAW
##    pre-QC, pre-host-removal read total (GMrepo doesn't expose a
##    post-filtering depth), so it overstates the true denominator by a
##    sample-varying amount. Treat these as approximate counts for tools
##    that require count-shaped input, not as ground truth.
## ---------------------------------------------------------------------

abundance_list_counts <- setNames(
  lapply(names(abundance_list_final), function(key) {
    mat   <- abundance_list_final[[key]]
    md    <- metadata_list_final[[key]]
    depth <- md$nr_reads_sequenced[match(rownames(mat), md$run_id)]
    round(sweep(mat, MARGIN = 1, STATS = depth / 100, FUN = "*"))
  }),
  names(abundance_list_final)
)

saveRDS(abundance_list_counts, file.path(out_dir, "abundance_list_counts.rds"))

## =========================================================================
## Part 3 (from Analysis3.R): build rel.abd / covariate.interest /
## covariate.adjust / cluster in the target meta-analysis format
## =========================================================================

count_list <- abundance_list_counts
study_ids  <- names(count_list)

## ---------------------------------------------------------------------
## 1. rel.abd: one sample x feature count matrix per study (study = one
##    project-disease pair), no missing values by construction.
## ---------------------------------------------------------------------

rel.abd <- setNames(lapply(study_ids, function(key) count_list[[key]]), study_ids)

stopifnot(all(sapply(rel.abd, function(m) !anyNA(m))))

## ---------------------------------------------------------------------
## 2. covariate.interest: one sample x 1 numeric matrix per study (disease
##    status: 1 = disease/case, 0 = health/control), row order matched to
##    rel.abd for that study.
## ---------------------------------------------------------------------

covariate.interest <- setNames(lapply(study_ids, function(key) {
  md           <- metadata_list_imputed[[key]]
  sample_order <- rownames(rel.abd[[key]])
  group        <- md$sample_group[match(sample_order, md$run_id)]

  matrix(as.numeric(group == "disease"), ncol = 1,
         dimnames = list(sample_order, "disease"))
}), study_ids)

stopifnot(all(sapply(covariate.interest, function(m) !anyNA(m))))

## ---------------------------------------------------------------------
## 3. covariate.adjust: one sample x confounder data frame per study, using
##    only the confounders found significant for that study (Part 2 above),
##    pulled from the imputed (no-missing) metadata so the "no NA"
##    requirement holds. Continuous confounders -> numeric, categorical ->
##    factor. Studies with no significant confounder get a 0-column data
##    frame (still row-aligned to rel.abd) rather than NULL, so downstream
##    code can always assume one data frame per study exists.
## ---------------------------------------------------------------------

covariate.adjust <- setNames(lapply(study_ids, function(key) {
  md           <- metadata_list_imputed[[key]]
  sample_order <- rownames(rel.abd[[key]])
  md           <- md[match(sample_order, md$run_id), ]

  cols <- significant_confounders[[key]]
  if (is.null(cols)) {
    return(data.frame(row.names = sample_order))
  }

  df <- as.data.frame(lapply(cols, function(col) {
    x         <- md[[col]]
    numeric_x <- suppressWarnings(as.numeric(x))
    if (!anyNA(numeric_x)) numeric_x else factor(x)
  }))
  names(df)    <- cols
  rownames(df) <- sample_order
  df
}), study_ids)

stopifnot(all(sapply(covariate.adjust, function(df) ncol(df) == 0 || !anyNA(df))))

## ---------------------------------------------------------------------
## 4. cluster: correlated-sample groups per study.
##
##    GMrepo's own metadata has no patient/subject ID field, so as a
##    baseline, samples are clustered by `sample_id` (the BioSample
##    accession) - if the same BioSample appears on more than one run,
##    those runs are the *same physical sample* resequenced, definitively
##    correlated.
##
##    For projects manually confirmed (via their NCBI SRA RunTable, in
##    Project_meta/) to have real subject-level repeated sampling - i.e. the
##    same *person* contributing multiple, distinct BioSamples over time -
##    that subject ID is used instead, since it captures the true clustering
##    the BioSample check structurally cannot see. Every other project's
##    RunTable was checked and showed no such column, or a candidate column
##    that turned out to be 1:1 with rows (e.g. PRJNA751448's
##    Participant_ID) or a mislabeled group name rather than a subject ID
##    (e.g. PRJNA793776's "isolate" was actually a case/control description) -
##    so the BioSample fallback is what's actually used for those.
##
##    Per the target format, an entry is only included for studies that
##    actually have correlated samples; other studies are simply absent from
##    the list, matching "default NULL = all samples independent".
## ---------------------------------------------------------------------

subject_id_sources <- list(
  PRJNA993675 = list(file = "Project_meta/SraRunTable_PRJNA993675.csv", column = "subject"),
  PRJEB42155  = list(file = "Project_meta/SraRunTable_PRJEB42155.csv",  column = "host_subject_id"),
  PRJNA398089 = list(file = "Project_meta/SraRunTable_PRJNA398089.csv", column = "host_subject_id"),
  PRJNA604850 = list(file = "Project_meta/SraRunTable_PRJNA604850.csv", column = "host_subject_id")
)

subject_lookup <- lapply(subject_id_sources, function(src) {
  sra <- read.csv(file.path(in_dir, src$file), stringsAsFactors = FALSE, check.names = FALSE)
  setNames(sra[[src$column]], sra$Run)
})

cluster_candidates <- lapply(study_ids, function(key) {
  md           <- metadata_list_imputed[[key]]
  sample_order <- rownames(rel.abd[[key]])
  proj         <- strsplit(key, "_")[[1]][1]

  if (proj %in% names(subject_lookup)) {
    run_ids     <- md$run_id[match(sample_order, md$run_id)]
    cluster_ids <- subject_lookup[[proj]][run_ids]
  } else {
    cluster_ids <- md$sample_id[match(sample_order, md$run_id)]
  }

  if (anyDuplicated(cluster_ids) == 0) return(NULL)
  setNames(cluster_ids, sample_order)
})
names(cluster_candidates) <- study_ids

cluster <- cluster_candidates[!sapply(cluster_candidates, is.null)]

cat("Studies with correlated samples:\n")
for (key in names(cluster)) {
  n_runs   <- length(cluster[[key]])
  n_groups <- length(unique(cluster[[key]]))
  source   <- if (strsplit(key, "_")[[1]][1] %in% names(subject_lookup)) "subject ID" else "BioSample dup"
  cat(sprintf(" - %-25s %d runs -> %d clusters  [%s]\n", key, n_runs, n_groups, source))
}
if (length(cluster) == 0) cat(" (none found)\n")

## ---------------------------------------------------------------------
## Save
## ---------------------------------------------------------------------

saveRDS(rel.abd,             file.path(out_dir, "rel.abd.rds"))
saveRDS(covariate.interest,  file.path(out_dir, "covariate.interest.rds"))
saveRDS(covariate.adjust,    file.path(out_dir, "covariate.adjust.rds"))
saveRDS(cluster,             file.path(out_dir, "cluster.rds"))

})


# ============================================================================
#   PART 2 - PAN-TUMOR (HGMT)
# ============================================================================
#   Pan-tumor application (HGMT) - data preprocessing
#  Raw input:
#    Data/HGMT/13059_2025_3865_MOESM1_ESM.xlsx   Table S4 of the HGMT paper:
#                                                one row per project-batch, with
#                                                the assay type and case/control
#                                                phenotypes
#    Data/HGMT/Bacteria_<batch>.txt              genus profiles, downloaded from
#                                                https://mai.fudan.edu.cn/hgmt
#    Data/HGMT/selected_project_<batch>.txt      per-sample metadata
#    Data/HGMT/<project>.csv                     SRA run table, for read depth
#
#  Output (written under HGMT_analysis/Data/):
#    HGMT_16S.Rdata            otu_final, meta_final, covariate.adjust
#    HGMT_analysis.Rdata       otu_filter, meta_filter, covariate.adjust_filter,
#                              feature_ID
#
#  Steps
#    1. select the non-WGS batches from Table S4, one batch per project
#    2. per project: read the genus profile and metadata, rescale relative
#       abundance to depth-weighted pseudo-counts, split into one case-control
#       comparison per tumour phenotype
#    3. keep comparisons with at least 10 cases and 10 controls
#    4. drop samples with fewer than 2,000 assigned reads
#    5. choose the covariates to adjust for, per comparison, by a within-study
#       imbalance test
#    6. keep genera reaching 10% prevalence in at least 6 tumour types
#
#


local({
  ## ---- configuration ---------------------------------------------------------
  hgmt_dir     <- "./Data/HGMT"
  table_s4     <- "./Data/HGMT/13059_2025_3865_MOESM1_ESM.xlsx"
  out_raw      <- "./HGMT_analysis/Data/HGMT_16S.Rdata"
  out_filtered <- "./HGMT_analysis/Data/HGMT_analysis.Rdata"

  min_group_n  <- 10     # minimum cases and controls per comparison
  depth_cut    <- 2000   # minimum assigned reads per sample
  prev_cut     <- 0.1    # within-study prevalence a genus must reach
  min_types    <- 6      # number of tumour types it must reach it in
  adjust_alpha <- 0.05   # imbalance p-value below which a covariate is adjusted

  ## MeSH-style phenotype names -> the labels used throughout the analysis.
  phenotype_label <- c(
    "Lung Neoplasms"                                = "Lung Cancer",
    "Adenomatous Polyps"                            = "Colorectal Polyps",
    "Colonic Polyps"                                = "Colorectal Polyps",
    "Breast Neoplasms"                              = "Breast Cancer",
    "Brain Neoplasms"                               = "Brain Metastasis",
    "Carcinoma, Hepatocellular"                     = "Hepatocellular Carcinoma",
    "Carcinoma, Non-Small-Cell Lung"                = "Non-Small-Cell Lung Cancer",
    "Carcinoma, Pancreatic Ductal"                  = "Pancreatic Ductal Adenocarcinoma",
    "Colorectal Neoplasms"                          = "Colorectal Cancer",
    "Endometrial Neoplasms"                         = "Endometrial Cancer",
    "Pancreatic Neoplasms"                          = "Pancreatic Cancer",
    "Precursor Cell Lymphoblastic Leukemia-Lymphoma" = "Acute Lymphoblastic Leukemia",
    "Stomach Neoplasms"                             = "Gastric Cancer",
    "Thyroid Neoplasms"                             = "Thyroid Cancer")

  ## Tumour types excluded from the reported analysis.
  drop_types <- c("Lung Cancer", "Breast Diseases")


# =============================================================================
#  1. Batches to read
# =============================================================================
#  Table S4 lists one row per project-batch.  WGS batches are dropped, and where
#  a project still has several batches only the first is kept, so each project
#  contributes once.

  table_S4 <- read_excel(table_s4, sheet = "Table S4")
  tab_S4   <- table_S4 %>% filter(Assay_type != "WGS")

  proj_ids  <- unique(tab_S4$batch)
  base_proj <- str_extract(proj_ids, "PRJ[A-Z0-9]+")
  proj_ids_single <- proj_ids[!duplicated(base_proj)]
  message("Table S4: ", nrow(table_S4), " rows -> ", nrow(tab_S4),
          " non-WGS -> ", length(proj_ids_single), " projects")


# =============================================================================
#  2. Read each project and split into case-control comparisons
# =============================================================================

  meta_lst <- list(); otu_lst <- list()

  for (pid in proj_ids_single) {

    if (pid == "PRJNA531273-PRJNA397112") {
      ## This project ships its disease and health arms as separate files.
      Bacteria <- rbind(read.delim(file.path(hgmt_dir, "Bacteria_PRJNA531273-PRJNA397112_disease.txt")),
                        read.delim(file.path(hgmt_dir, "Bacteria_PRJNA531273-PRJNA397112_Health.txt")))
      meta <- rbind(
        read_tsv(file.path(hgmt_dir, "selected_project_PRJNA531273-PRJNA397112_disease.txt"),
                 skip = 1, show_col_types = FALSE),
        read_tsv(file.path(hgmt_dir, "selected_project_PRJNA531273-PRJNA397112_Health.txt"),
                 skip = 1, show_col_types = FALSE))
      seqdat <- rbind(
        read_csv(file.path(hgmt_dir, "PRJNA531273-PRJNA397112_disease.csv"), show_col_types = FALSE) %>%
          transmute(Run, depth = Bases / AvgSpotLen),
        read_csv(file.path(hgmt_dir, "PRJNA531273-PRJNA397112_Health.csv"), show_col_types = FALSE) %>%
          transmute(Run, depth = Bases / AvgSpotLen)) %>%
        column_to_rownames("Run")

    } else {
      Bacteria <- read.delim(file.path(hgmt_dir, paste0("Bacteria_", pid, ".txt")))

      ## The SRA run table is per project, so strip any batch suffix.
      seqdat <- read_csv(file.path(hgmt_dir, paste0(strsplit(pid, "_")[[1]][1], ".csv")),
                         show_col_types = FALSE)
      seqdat <- if ("AvgSpotLen" %in% colnames(seqdat)) {
        if (pid %in% c("PRJEB6070_16S", "PRJEB6070_PE", "PRJEB6070_SE")) {
          ## AvgSpotLen is unusable for this project; Bytes tracks depth instead.
          seqdat %>% transmute(Run, depth = Bytes) %>% column_to_rownames("Run")
        } else {
          seqdat %>% transmute(Run, depth = Bases / AvgSpotLen) %>% column_to_rownames("Run")
        }
      } else {
        seqdat %>% transmute(Run, depth = Bases) %>% column_to_rownames("Run")
      }
      seqdat <- seqdat %>% filter(!is.na(depth))

      meta <- read_tsv(file.path(hgmt_dir, paste0("selected_project_", pid, ".txt")),
                       skip = 1, show_col_types = FALSE)

      ## Some batch files mix assays; keep the majority one so a comparison is
      ## not built from two platforms at once.
      if (length(unique(meta$`Assay type`)) > 1) {
        tmp <- table(meta$`Assay type`)
        meta <- meta[meta$`Assay type` == names(tmp)[which.max(tmp)], ]
      }
    }

    ## ---- align metadata and profiles
    common_runs  <- intersect(meta$`Run ID`, unique(Bacteria$Run.ID))
    meta_sub     <- meta %>% filter(`Run ID` %in% common_runs)
    bacteria_sub <- Bacteria[Bacteria$Run.ID %in% common_runs, ]

    ## ---- collapse the lineage string to genus and pivot to a matrix
    genus_mat <- bacteria_sub %>%
      mutate(genus = str_remove(str_extract(Taxa, "g__[^|]+"), "g__"),
             genus = ifelse(is.na(genus), "Unknown_genus", genus)) %>%
      group_by(Run.ID, genus) %>%
      summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = genus, values_from = Abundance, values_fill = 0) %>%
      column_to_rownames("Run.ID") %>%
      as.matrix()

    ## ---- relative abundance -> depth-weighted pseudo-counts
    shared_ID <- intersect(rownames(genus_mat), rownames(seqdat))
    genus_mat <- genus_mat[shared_ID, ] * seqdat[shared_ID, "depth"]
    meta_sub  <- meta_sub[meta_sub$`Run ID` %in% shared_ID, ]

    ## One project carries a set of very deep runs that are not comparable with
    ## the rest of its samples.
    if (pid == "PRJEB26531_16S") {
      keep_ids  <- setdiff(rownames(genus_mat),
                           shared_ID[seqdat[shared_ID, "depth"] >= 1e7])
      genus_mat <- genus_mat[keep_ids, ]
      meta_sub  <- meta_sub[meta_sub$`Run ID` %in% keep_ids, ]
    }

    ## ---- one case-control comparison per tumour phenotype
    if (!"Health" %in% meta_sub$`Phenotype name`) next
    for (pdis in setdiff(unique(meta_sub$`Phenotype name`), "Health")) {
      plabel <- if (pdis %in% names(phenotype_label)) phenotype_label[[pdis]] else pdis
      if (plabel %in% drop_types) next

      rows <- meta_sub$`Phenotype name` %in% c(pdis, "Health")
      key  <- paste0(pid, "_", plabel)
      meta_lst[[key]] <- meta_sub[rows, ]
      otu_lst[[key]]  <- genus_mat[meta_sub$`Run ID`[rows], ]
    }
  }
  message("case-control comparisons built: ", length(otu_lst))


# =============================================================================
#  3-4. Group-size and depth filters
# =============================================================================

  meta_lst <- keep(meta_lst, function(df) {
    tab <- table(df$`Phenotype name`)
    n_health  <- ifelse(is.na(tab["Health"]), 0, tab["Health"])
    n_disease <- sum(tab[names(tab) != "Health"])
    n_health >= min_group_n && n_disease >= min_group_n
  })
  otu_lst <- otu_lst[names(meta_lst)]
  message("after the >=", min_group_n, "/", min_group_n, " filter: ", length(otu_lst))

  otu_final <- map(otu_lst, function(X) {
    X <- as.matrix(X)
    if (nrow(X) && ncol(X)) X <- X[rowSums(X, na.rm = TRUE) >= depth_cut, , drop = FALSE]
    X
  })
  otu_final  <- keep(otu_final, ~ ncol(.x) > 0 && nrow(.x) > 0)
  meta_final <- mapply(function(m, o) m[m$`Run ID` %in% rownames(o), , drop = FALSE],
                       meta_lst[names(otu_final)], otu_final, SIMPLIFY = FALSE)
  message("after the depth filter: ", length(otu_final), " comparisons, ",
          sum(vapply(otu_final, nrow, integer(1))), " samples")


# =============================================================================
#  5. Covariates to adjust for, chosen per comparison
# =============================================================================
#  A covariate is adjusted for only where it is recorded and differs between
#  cases and controls; adjusting on a covariate that is balanced adds noise, and
#  one that is entirely missing cannot be used at all.

  check_confounders <- function(meta,
                                vars_quant = c("age", "BMI"),
                                vars_qual  = c("country", "sex")) {
    meta <- meta %>%
      filter(!is.na(`Phenotype name`)) %>%
      mutate(DiseaseStatus = ifelse(`Phenotype name` == "Health", "Health", "Disease"))

    bind_rows(lapply(c(vars_quant, vars_qual), function(v) {
      blank <- function(x) { x <- as.character(x); x[x %in% c("", " ", "NA", "N/A", "na")] <- NA; x }
      if (!v %in% colnames(meta))
        return(tibble(variable = v, missing_status = "not_present", p_value = NA_real_, adjust = FALSE))
      x <- blank(meta[[v]])
      if (all(is.na(x)))
        return(tibble(variable = v, missing_status = "all_NA", p_value = NA_real_, adjust = FALSE))

      st <- ifelse(any(is.na(x)), "partial_NA", "complete")
      d  <- meta %>% mutate(tmp = x) %>% filter(!is.na(tmp))
      if (n_distinct(d$DiseaseStatus) < 2)
        return(tibble(variable = v, missing_status = st, p_value = NA_real_, adjust = FALSE))

      p <- if (v %in% vars_qual) {
        tryCatch(fisher.test(table(d$DiseaseStatus, d$tmp))$p.value, error = function(e) NA_real_)
      } else {
        d <- d %>% mutate(tmp = suppressWarnings(as.numeric(tmp))) %>% filter(!is.na(tmp))
        if (!nrow(d) || n_distinct(d$DiseaseStatus) < 2) NA_real_
        else tryCatch(wilcox.test(tmp ~ DiseaseStatus, data = d)$p.value, error = function(e) NA_real_)
      }
      tibble(variable = v, missing_status = st, p_value = p,
             adjust = !is.na(p) && p <= adjust_alpha)
    }))
  }

  conf <- imap(meta_final, ~ check_confounders(.x) %>% mutate(Study = .y, .before = 1)) %>%
    bind_rows()

  covariate.adjust <- list()
  for (d in names(otu_final)) {
    ci <- conf %>% filter(Study == d)
    ## A partially missing covariate is dropped: imputing it would invent the
    ## very imbalance the test just detected.
    use <- ci$variable[ci$adjust & ci$missing_status == "complete"]
    if (!length(use)) next

    cv <- meta_final[[d]] %>% column_to_rownames("Run ID")
    cv <- as.data.frame(cv[, use, drop = FALSE])

    ## A numeric covariate can still hold a value that does not coerce - one BMI
    ## entry in PRJNA290926 is such a case - which would otherwise drop that
    ## sample from the association model.  Fill those with the column mean.
    for (v in intersect(use, c("age", "BMI"))) {
      cv[[v]] <- suppressWarnings(as.numeric(cv[[v]]))
      cv[[v]][is.na(cv[[v]])] <- mean(cv[[v]], na.rm = TRUE)
    }
    if ("country" %in% use) cv$country <- as.factor(cv$country)
    covariate.adjust[[d]] <- cv
  }
  message("comparisons with covariate adjustment: ", length(covariate.adjust),
          " of ", length(otu_final))

  save(otu_final, meta_final, covariate.adjust, file = out_raw)


# =============================================================================
#  6. Genus prevalence filter, applied across tumour types
# =============================================================================

  tumour <- str_remove(names(otu_final), "^PRJ[A-Z0-9-]+(?:_(?:16S|WGS|PE|SE))?_")
  per_type <- lapply(unique(tumour), function(l)
    unique(unlist(lapply(otu_final[tumour == l],
                         function(d) colnames(d)[colMeans(d != 0) >= prev_cut]))))
  feature_ID <- names(which(table(unlist(per_type)) >= min_types))
  message("tumour types: ", length(unique(tumour)), " | genera retained: ", length(feature_ID))

  otu_filter <- lapply(otu_final, function(d) {
    x <- d[, intersect(colnames(d), feature_ID), drop = FALSE]
    x[rowSums(x) > 0, , drop = FALSE]
  })
  meta_filter <- lapply(names(otu_filter), function(d) {
    m <- meta_final[[d]]; rownames(m) <- m$`Run ID`; m[rownames(otu_filter[[d]]), ]
  })
  names(meta_filter) <- names(otu_filter)

  covariate.adjust_filter <- lapply(names(otu_filter), function(l) {
    if (is.null(covariate.adjust[[l]])) NULL
    else covariate.adjust[[l]][rownames(otu_filter[[l]]), , drop = FALSE]
  })
  names(covariate.adjust_filter) <- names(otu_filter)
  covariate.adjust_filter <- compact(covariate.adjust_filter)

  save(otu_filter, meta_filter, covariate.adjust_filter, feature_ID, file = out_filtered)

  message("\nwritten:")
  message("  ", out_raw)
  message("  ", out_filtered, "  (", length(otu_filter), " comparisons x ",
          length(feature_ID), " genera)")

})


# ============================================================================
#   PART 3 - COLORECTAL NEOPLASIA
# ============================================================================
#   Colorectal neoplasia application - data preprocessing
#  Raw input (all under ./Data/):
#    Data/CRC/sequence/            per-study SRA / ENA run tables, used for read depth
#    Data/CRC/metadata/            per-study sample metadata (17 studies)
#    Data/CRC/metaphlan4_profiles/ per-study MetaPhlAn 4 relative-abundance profiles
#
#  Output:
#    CRC_analysis/Data/CRC_strata.Rdata   otu_data_sp, meta_data_sp (26 strata)
#    CRC_analysis/Data/CRC_analysis.Rdata otu_filter, meta_filter, depth
#
#  Steps
#    1. read sequencing depth per sample
#    2. read metadata and attach the read count
#    3. read MetaPhlAn profiles, rescale to depth-weighted pseudo-counts, keep
#       species-level features
#    4. split each study into Adenoma / early / late strata against its controls
#    5. keep strata with at least 10 samples in every disease category
#    6. partition each study's controls across its surviving strata
#    7. keep species present in >=10% of samples in at least 23 strata
#
#
#  ---------------------------------------------------------------------------
#  WHAT THIS REPRODUCES.  The stratum set, the case counts and the control
#  counts all match the published analysis exactly:
#  26 strata, 2,679 samples.
#
#  What it cannot reproduce is WHICH controls land in which stratum.  The
#  published partition was drawn at random with a seed that was not recorded, so
#  roughly three quarters of each stratum's samples coincide with the published
#  one and the rest differ.  Because within-stratum prevalence shifts slightly as
#  a result, step 7 keeps 276 species here against the 260 analysed.
#
#  Outputs therefore use NEW file names and never overwrite the analysed
#  objects.  Re-point `out_dir` only if you intend to replace them.


local({
  ## ---- configuration ---------------------------------------------------------
  seq_dir   <- "./Data/CRC/sequence"
  meta_dir  <- "./Data/CRC/metadata"
  prof_dir  <- "./Data/CRC/metaphlan4_profiles"
  out_dir   <- "./CRC_analysis/Data"

  min_group_n <- 10     # minimum samples per disease category within a stratum
  pool_prev_cut <- 0.2  # prevalence across ALL pooled samples, before strata
  prev_cut    <- 0.1    # within-stratum prevalence a species must reach
  min_strata  <- 23     # number of strata a species must reach it in

  early_stage <- c("0", "I", "II")
  late_stage  <- c("III", "IV")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
#  1. Sequencing depth
# =============================================================================
#  The run tables come from two archives with different column names, and a few
#  studies key their samples differently again; each branch maps that study's
#  identifier onto the sample id used in the metadata.

  seqdata <- NULL
  for (d in dir(seq_dir)) {
    p <- file.path(seq_dir, d)

    if (length(strsplit(d, "_tsv")[[1]]) == 2) {          # ENA file report
      tmp <- read.delim(p, row.names = 1)
      tmp <- if (d == "filereport_read_run_PRJEB7774_tsv.txt") {
        tmp %>% transmute(sample_ID = paste0("SID", sample_alias), reads = read_count)
      } else {
        tmp %>% transmute(sample_ID = sample_alias, reads = read_count)
      }
    } else {                                              # SRA run table
      tmp <- read.csv(p)
      tmp <- if (d %in% c("SraRunTable_PRJNA531273.csv", "SraRunTable_PRJNA397112.csv")) {
        tmp %>% transmute(sample_ID = paste0("GupDM_", Sample.Name), reads = round(Bases / AvgSpotLen))
      } else if (d %in% c("SraRunTable_PRJNA447983.csv", "SraRunTable_PRJNA1167935.csv")) {
        tmp %>% transmute(sample_ID = Sample.Name, reads = round(Bases / AvgSpotLen))
      } else if (d == "SraRunTable_PRJNA1237248.csv") {
        tmp %>% transmute(sample_ID = Library.Name, reads = round(Bases / AvgSpotLen))
      } else {
        tmp %>% transmute(sample_ID = BioSample, reads = round(Bases / AvgSpotLen))
      }
    }
    seqdata <- rbind(seqdata, tmp)
  }
  message("depth records: ", nrow(seqdata))


# =============================================================================
#  2. Metadata, with the read count attached
# =============================================================================
#  A sample without a depth record cannot be rescaled to counts, so the inner
#  join here is also the first sample-level filter.

  metadata <- list()
  for (d in dir(meta_dir)) {
    id  <- sub(".tsv", "", sub(".*__", "", d))
    tmp <- read.delim(file.path(meta_dir, d), row.names = 1)

    ## ThomasAM_2018b keys its profiles by subject, not by the metadata rowname.
    if (d == "Public_study__ThomasAM_2018b.tsv")
      rownames(tmp) <- sub("^.*_(SBJ[^_]*)_.*$", "\\1", rownames(tmp))

    metadata[[id]] <- tmp %>%
      rownames_to_column("sample_ID") %>%
      inner_join(seqdata %>% group_by(sample_ID) %>%
                   summarise(read_count = sum(reads), .groups = "drop"),
                 by = "sample_ID") %>%
      column_to_rownames("sample_ID")
  }
  message("studies with metadata: ", length(metadata))


# =============================================================================
#  3. MetaPhlAn profiles -> species-level pseudo-counts
# =============================================================================
#  MetaPhlAn reports relative abundance; multiplying by the sample's read count
#  puts every study on a common count-like scale for the association models.

  countdata_species <- list()
  for (d in dir(prof_dir)) {
    id  <- sub(".tsv", "", sub(".*__", "", d))
    cnt <- read.delim(file.path(prof_dir, d), row.names = 1)

    if (d == "Public_study__ThomasAM_2018b.tsv")
      colnames(cnt) <- sub("^.*_(SBJ[^_]*)_.*$", "\\1", colnames(cnt))
    colnames(cnt) <- gsub("\\.", "-", colnames(cnt))

    match_id <- rownames(metadata[[id]])
    cnt <- (t(cnt) / colSums(cnt))[match_id, ] * metadata[[id]][, "read_count"]

    sp <- cnt[, str_which(colnames(cnt), "\\|s__[^|]+$"), drop = FALSE]
    colnames(sp) <- paste0("s_", sub(".*s__", "", colnames(sp)))
    countdata_species[[id]] <- sp
  }
  message("studies with profiles: ", length(countdata_species))

# =============================================================================
#  3b. Pooled prevalence filter
# =============================================================================
#  Every sample from every study is stacked into one matrix and a species kept
#  only if present in at least `pool_prev_cut` of them.  This runs BEFORE the
#  strata are formed, so all strata share one candidate species set.

  uni.features_d <- unique(unlist(lapply(countdata_species, colnames)))
  pool_data <- NULL
  for (d in names(countdata_species)) {
    tmp_count <- matrix(0, nrow = nrow(countdata_species[[d]]), ncol = length(uni.features_d),
                        dimnames = list(rownames(countdata_species[[d]]), uni.features_d))
    tmp_count[, colnames(countdata_species[[d]])] <- countdata_species[[d]]
    pool_data <- rbind(pool_data, tmp_count)
  }
  uni.features <- names(which(colMeans(pool_data != 0) >= pool_prev_cut))
  message("pooled prevalence >= ", pool_prev_cut, ": ", length(uni.features),
          " of ", length(uni.features_d), " species")

  for (d in names(countdata_species))
    countdata_species[[d]] <- countdata_species[[d]][
      , intersect(colnames(countdata_species[[d]]), uni.features), drop = FALSE]
  rm(pool_data)


# =============================================================================
#  4. Study-stage strata
# =============================================================================
#  A stratum is one neoplasia stage within one study, compared against that
#  study's controls.  Studies without a control group, or without AJCC staging,
#  contribute nothing.

  otu_data_sp <- list(); meta_data_sp <- list()

  for (d in names(metadata)) {
    md <- metadata[[d]]
    if (!"Control" %in% md$Disease) next          # no control group

    take <- function(keep, tag) {
      m <- md %>%
        transmute(SampleID = Sample.Source.ID, Disease, Age, BMI, Sex,
                  Group = as.numeric(Disease != "Control"),
                  Tumor.Staging.AJCC, Primary.Tumor.Location, Study.name) %>%
        filter(keep)
      key <- paste0(d, ":", tag)
      meta_data_sp[[key]] <<- m
      otu_data_sp[[key]]  <<- countdata_species[[d]][m$SampleID, , drop = FALSE]
    }

    if ("Adenoma" %in% md$Disease)
      take(md$Disease %in% c("Adenoma", "Control"), "Adenoma")

    ## NB: the original code tested `early_stage` for BOTH branches; kept as is
    ## so the stratum set matches the published analysis.
    if (any(early_stage %in% md$Tumor.Staging.AJCC)) {
      take(md$Disease %in% c("colorectal carcinoma", "Control") &
             (md$Tumor.Staging.AJCC %in% early_stage | md$Disease == "Control"), "early")
      take(md$Disease %in% c("colorectal carcinoma", "Control") &
             (md$Tumor.Staging.AJCC %in% late_stage  | md$Disease == "Control"), "late")
    }
  }


# =============================================================================
#  5. Drop strata with fewer than `min_group_n` in any disease category
# =============================================================================
#  Run BEFORE the control partition below, so a study whose other stages are too
#  small gives its whole control set to the stage that survives.

  keep <- vapply(meta_data_sp, function(m) all(table(m$Disease) >= min_group_n), logical(1))
  meta_data_sp <- meta_data_sp[keep]
  otu_data_sp  <- otu_data_sp[keep]
  message("strata retained: ", length(otu_data_sp), " of ", length(keep))


# =============================================================================
#  6. Partition each study's controls across its surviving strata
# =============================================================================
#  Step 4 gives every stratum its study's full control set, which would make the
#  stage strata of one study share controls and so cease to be independent.  The
#  controls are therefore split: each is used in exactly one stratum, allocated
#  in proportion to that stratum's case count, with a floor of `min_group_n` so
#  no stratum is left with fewer controls than the inclusion rule requires.
#
#  Seed 20252025 is the one the published data was produced with, so this
#  reproduces the original assignment sample for sample.

  ## ---- per-stratum prevalence filter, BEFORE the partition -------------------
  ## Prevalence is computed on each study's FULL control set.  Filtering after the
  ## split would use the smaller per-stratum control groups and let more species
  ## through, so the order matters.

  for (d in names(otu_data_sp)) {
    tmp_data <- otu_data_sp[[d]]
    otu_data_sp[[d]] <- tmp_data[, colMeans(tmp_data != 0) >= prev_cut, drop = FALSE]
  }

  ## ---- allocate the controls -------------------------------------------------
  ## Reproduces the published assignment exactly: same seed, same rounding
  ## (`round`, floor of `min_group_n`, remainder taken off the largest share) and
  ## the same study order, which matters because one sample() call is made per
  ## study and each consumes the RNG stream in turn.

  partition_seed <- 20252025
  set.seed(partition_seed)

  multi_study <- names(which(table(sub(":.*$", "", names(meta_data_sp))) > 1))
  for (d in multi_study) {
    tmp_id <- which(sub(":.*$", "", names(meta_data_sp)) == d)
    cont   <- sapply(meta_data_sp[tmp_id], function(x) table(x$Disease))[1, ]
    cont_c <- sapply(meta_data_sp[tmp_id], function(x) table(x$Disease))[2, 1]
    prop   <- round(cont / sum(cont) * cont_c)
    prop   <- pmax(prop, min_group_n)
    prop[which.max(prop)] <- prop[which.max(prop)] - (sum(prop) - cont_c)

    Control_ID <- sample(meta_data_sp[[tmp_id[1]]]$SampleID[
      meta_data_sp[[tmp_id[1]]]$Disease == "Control"])
    for (dd in seq_along(prop)) {
      cases <- meta_data_sp[[tmp_id[dd]]]$SampleID[
        meta_data_sp[[tmp_id[dd]]]$Disease != "Control"]
      ctrls <- if (dd == 1) Control_ID[seq_len(prop[1])]
               else Control_ID[(1 + sum(prop[1:(dd-1)])):sum(prop[1:dd])]
      IDs <- c(cases, ctrls)
      m <- meta_data_sp[[tmp_id[dd]]]
      rownames(m) <- as.character(m$SampleID)
      meta_data_sp[[tmp_id[dd]]] <- m[as.character(IDs), ]
      otu_data_sp[[tmp_id[dd]]]  <- otu_data_sp[[tmp_id[dd]]][as.character(IDs), , drop = FALSE]
    }
  }

  message("samples after partition: ", sum(vapply(meta_data_sp, nrow, integer(1))),
          " (each control used once)")

  save(otu_data_sp, meta_data_sp, file = file.path(out_dir, "CRC_strata.Rdata"))


# =============================================================================
#  7. Species prevalence filter, applied across strata
# =============================================================================

  feature_ID <- names(which(table(unlist(lapply(otu_data_sp, function(d)
    colnames(d)[colMeans(d != 0) >= prev_cut]))) >= min_strata))
  message("species retained: ", length(feature_ID))

  otu_filter <- lapply(otu_data_sp, function(d) {
    x <- d[, intersect(colnames(d), feature_ID), drop = FALSE]
    x[rowSums(x) > 0, , drop = FALSE]
  })

  meta_filter <- list(); depth <- list()
  for (d in names(otu_filter)) {
    meta_filter[[d]] <- meta_data_sp[[d]][rownames(otu_filter[[d]]), ]
    depth[[d]]       <- rowSums(otu_data_sp[[d]][rownames(otu_filter[[d]]), , drop = FALSE])
  }

  save(otu_filter, meta_filter, depth, file = file.path(out_dir, "CRC_analysis.Rdata"))

  message("\nwritten to ", out_dir, ":")
  message("  CRC_strata.Rdata    ", length(otu_data_sp), " strata")
  message("  CRC_analysis.Rdata  ", length(otu_filter), " strata x ", length(feature_ID), " species")

})
