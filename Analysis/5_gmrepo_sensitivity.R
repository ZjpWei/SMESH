# =============================================================================
#   5_gmrepo_sensitivity.R  --  pan-disease sensitivity analyses
# =============================================================================
#  Two checks on the pan-disease (GMrepo) result, both re-running the whole
#  pipeline with one thing changed and everything else held fixed:
#
#    taxonomy   counts aggregated to GENUS instead of species.  Tests whether
#               the clustering survives a coarser taxonomic resolution.
#    filtering  species level, but with a STRINGENT feature filter: a species
#               must be estimable in more than X = 4 of the 12 diseases.  Tests
#               whether the result depends on the feature-inclusion cutoff.
#
#  Both feed panel D of Supplementary Figure S3 (Rscript/FigureS3.R), which
#  reports the ARI between each sensitivity clustering and the main one.
#
#  This merges Analysis5.R (taxonomy) and Analysis6.R (filtering), which are now
#  retired to Backup/gmrepo_scripts/.
#  The two were structurally identical apart from their inputs, their feature-set
#  rule and their output names, so those differences live in the ANALYSES table
#  below and nothing else varies.  No method parameter has been changed.
#
#  Pipeline, per differential-abundance back-end:
#    1. per project-disease summary statistics (beta-hat, SE)
#    2. PALM meta-analysis within each disease -> context-level summaries
#    3. mask features outside the dense feature set
#    4. structure-learning meta-analysis at a FIXED G = 2
#
#  Methods: SMESH-PALM, SMESH-ANCOMBC2, SMESH-MaAsLin3, SMESH-LinDA (all four
#  are SMESH run on different summary statistics), plus SKM+FE and SHC+FE,
#  which cluster the same context-level summaries with sparse k-means / sparse
#  hierarchical clustering and then fit a fixed-effect model per cluster.
#
#  Usage, from the project root:
#      Rscript Analysis/5_gmrepo_sensitivity.R <analysis>
#      Rscript Analysis/5_gmrepo_sensitivity.R taxonomy
#      Rscript Analysis/5_gmrepo_sensitivity.R filtering
#
#  Runtime is dominated by the differential-abundance stage (~20 min), so the
#  run is checkpointed twice.  Both checkpoints are resumed automatically when
#  present; delete them to force a stage to recompute.
# =============================================================================

  rm(list = ls())

  library(dplyr)
  library(tibble)
  library(sparcl)
  library(Matrix)
  ## utility/PALM_tune.R and utility/SMESH.R call abess() and glmtlp()
  ## unqualified, so both have to be attached here.
  library(abess)
  library(glmtlp)

  set.seed(2026)

  ## ---- per-analysis settings --------------------------------------------------
  ## These are the only differences between the two original scripts.  The
  ## permutation counts are NOT uniform - they are carried over exactly as each
  ## analysis was run.
  ANALYSES <- list(
    taxonomy = list(
      unit        = "genera",
      inputs      = c(rel.abd            = "rel.abd_genus.rds",
                      covariate.interest = "covariate.interest_genus.rds",
                      covariate.adjust   = "covariate.adjust_genus.rds",
                      cluster            = "cluster_genus.rds"),
      out_all     = "res_genus_G2.Rdata",
      ## second, two-object file that Species_genus_consistency.R still reads
      out_comp    = "res_genus.Rdata",
      ckpt_stats  = "res_genus_summary_stats.Rdata",
      ckpt_fits   = "res_genus_smesh_fits.Rdata",
      maaslin_out = "MaAsLin3_genus",
      doc_prefix  = "smesh_genus_",
      nperm_smesh = 30,
      X           = NULL),
    filtering = list(
      unit        = "species",
      inputs      = c(rel.abd            = "rel.abd.rds",
                      covariate.interest = "covariate.interest.rds",
                      covariate.adjust   = "covariate.adjust.rds",
                      cluster            = "cluster.rds"),
      out_all     = "res_species_X4_G2.Rdata",
      out_comp    = NULL,
      ckpt_stats  = "res_species_X4_summary_stats.Rdata",
      ckpt_fits   = "res_species_X4_smesh_fits.Rdata",
      maaslin_out = "MaAsLin3_X4",
      doc_prefix  = "smesh_X4_",
      nperm_smesh = 20,
      X           = 4)   # keep species estimable in MORE than X diseases
  )

  ## ---- inputs -----------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1)
    stop("need <analysis>, one of: ", paste(names(ANALYSES), collapse = ", "), "\n",
         "  Rscript Analysis/5_gmrepo_sensitivity.R taxonomy")
  analysis <- args[1]
  if (!analysis %in% names(ANALYSES))
    stop("unknown analysis '", analysis, "'; choose ",
         paste(names(ANALYSES), collapse = ", "))

  if (!dir.exists("utility") || !dir.exists("GMrepo_analysis"))
    stop("run this from the project root, e.g. Rscript Analysis/5_gmrepo_sensitivity.R taxonomy")

  cfg <- ANALYSES[[analysis]]

  base_dir <- "./GMrepo_analysis/Sensitivity"  # sensitivity results and checkpoints
  proc_dir <- "./GMrepo_analysis/Data"  # processed objects from 0_preprocessing.R
  raw_dir  <- "./Data/GMrepo"           # raw download and lookup tables
  out_all  <- file.path(base_dir, cfg$out_all)
  out_comp <- if (is.null(cfg$out_comp)) NULL else file.path(base_dir, cfg$out_comp)

  G           <- 2      # fixed for every method
  prev_cut    <- 0.1    # within-dataset prevalence filter
  nperm_smesh <- cfg$nperm_smesh
  X           <- cfg$X  # NULL for the taxonomy check, which uses a lookup instead
  maaslin_out <- file.path(base_dir, cfg$maaslin_out)

  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(maaslin_out, showWarnings = FALSE, recursive = TRUE)

  ckpt_stats <- file.path(base_dir, cfg$ckpt_stats)
  ckpt_fits  <- file.path(base_dir, cfg$ckpt_fits)

  message(analysis, " sensitivity | unit = ", cfg$unit, " | G = ", G,
          " | nperm = ", nperm_smesh,
          if (!is.null(X)) paste0(" | X = ", X) else "")

  rel.abd            <- readRDS(file.path(proc_dir, cfg$inputs[["rel.abd"]]))
  covariate.interest <- readRDS(file.path(proc_dir, cfg$inputs[["covariate.interest"]]))
  covariate.adjust   <- readRDS(file.path(proc_dir, cfg$inputs[["covariate.adjust"]]))
  cluster            <- readRDS(file.path(proc_dir, cfg$inputs[["cluster"]]))

  studies <- names(rel.abd)

  ## study key "<project_id>_<disease_mesh>" -> disease name, so studies of the
  ## same disease can be meta-analysed together.
  project_summary <- openxlsx::read.xlsx(file.path(proc_dir, "GMrepo_project_summary.xlsx"))
  mesh_to_disease <- unique(project_summary[, c("disease_mesh", "disease_name")])
  disease_info <- setNames(
    mesh_to_disease$disease_name[match(sub("^.*_", "", studies), mesh_to_disease$disease_mesh)],
    studies
  )
  stopifnot(!anyNA(disease_info))
  diseases <- unique(disease_info)


# =============================================================================
#  Helpers
# =============================================================================

  ## Counts and metadata for one project-disease, after the prevalence filter.
  dataset_of <- function(d) {
    Y    <- rel.abd[[d]]
    keep <- colMeans(Y != 0) >= prev_cut & apply(Y, 2, var) > 0
    Y    <- Y[, keep, drop = FALSE]

    meta <- data.frame(disease = as.integer(covariate.interest[[d]][, 1]),
                       row.names = rownames(covariate.interest[[d]]))
    adj <- as.data.frame(covariate.adjust[[d]])
    if (ncol(adj) > 0) meta <- cbind(meta, adj[rownames(meta), , drop = FALSE])

    list(Y = Y, meta = meta[rownames(Y), , drop = FALSE])
  }

  ## An empty (est, stderr) pair on the full feature universe, so every dataset
  ## returns the same shape regardless of which features it could fit.
  empty_sums <- function(features, n) {
    m <- matrix(NA_real_, length(features), 1, dimnames = list(features, "disease"))
    list(est = m, stderr = m, n = n)
  }

  ## PALM meta-analysis of the per-dataset summaries, within each disease.
  meta_by_disease <- function(sums, n_by_study) {
    out <- list()
    for (dty in diseases) {
      keep <- names(sums)[disease_info[names(sums)] == dty]
      res  <- PALM::palm.meta.summary(summary.stats = sums[keep])
      out[[dty]] <- list(
        est    = matrix(res$disease$coef, ncol = 1,
                        dimnames = list(res$disease$feature, "disease")),
        stderr = matrix(res$disease$stderr, ncol = 1,
                        dimnames = list(res$disease$feature, "disease")),
        n      = sum(n_by_study[keep])
      )
    }
    out
  }

  ## Blank out features outside the dense set, so every method is fitted on the
  ## same features.
  mask_to <- function(meta_sums, features) {
    lapply(meta_sums, function(s) {
      drop <- setdiff(rownames(s$est), features)
      if (length(drop)) {
        s$est[drop, "disease"]    <- NA
        s$stderr[drop, "disease"] <- NA
      }
      s
    })
  }

  ## --- the dense feature set: the one genuine methodological difference ------

  ## filtering: species estimable in more than X diseases.  Counted on the
  ## UNMASKED context-level summaries, since masking is what this decides.
  dense_features_of <- function(meta_sums, X) {
    feats <- unique(unlist(lapply(meta_sums, function(s) rownames(s$est))))
    n_dis <- vapply(feats, function(f) {
      sum(vapply(meta_sums, function(s) {
        f %in% rownames(s$est) && !is.na(s$est[f, "disease"])
      }, logical(1)))
    }, numeric(1))
    names(n_dis)[n_dis > X]
  }

  ## taxonomy: the genera implied by the species-level dense set, mapped through
  ## the NCBI lineage (not by parsing species names - "[Clostridium] symbiosum"
  ## is genus Otoolea).  Using the species-level set keeps the two resolutions
  ## comparable rather than each picking its own features.
  dense_genera_of <- function() {
    species_env <- new.env()
    load(file.path(proc_dir, "species_lst.Rdata"), envir = species_env)

    genus_lookup   <- read.csv(file.path(raw_dir, "taxon_id_genus_lookup.csv"),
                               stringsAsFactors = FALSE, colClasses = "character")
    taxid_to_genus <- setNames(genus_lookup$genus_name, genus_lookup$species_taxon_id)

    dense_genus <- taxid_to_genus[sub(".*\\[(\\d+)\\]$", "\\1", species_env$dense_feature)]
    dense_genus[is.na(dense_genus) | dense_genus == ""] <- "Unclassified"
    sort(unique(unname(dense_genus)))
  }


# =============================================================================
#  1-2. Summary statistics and the dense feature set
# =============================================================================
#  The differential-abundance stage is the expensive part, so it is checkpointed
#  at the end of section 2 and skipped when that file already exists.

if (!file.exists(ckpt_stats)) {

  n_by_study       <- vapply(studies, function(d) nrow(rel.abd[[d]]), numeric(1))
  feature_universe <- sort(unique(unlist(lapply(rel.abd, colnames))))

  ## ---- PALM ------------------------------------------------------------------
  message("[1/4] PALM summary statistics")
  null_model <- PALM::palm.null.model(rel.abd = rel.abd,
                                      covariate.adjust = covariate.adjust,
                                      prev.filter = prev_cut)

  summary_stat_PALM <- PALM::palm.get.summary(null.obj = null_model,
                                              covariate.interest = covariate.interest,
                                              cluster = cluster,
                                              correct = "tune")

  ## ---- ANCOM-BC2 -------------------------------------------------------------
  message("[2/4] ANCOM-BC2 summary statistics")
  source("./utility/ancombc.R")
  summary_stat_ANCOMBC2 <- lapply(studies, function(d) {
    dat <- dataset_of(d)
    s   <- empty_sums(feature_universe, nrow(dat$Y))
    fit <- try(ancombc.fun(feature.table = t(dat$Y), meta = dat$meta,
                           formula = paste(colnames(dat$meta), collapse = " + "),
                           adjust.method = "fdr", group = NULL, subject = NULL,
                           method = "ancombc2"), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      ok <- fit$res$passed_ss_disease
      s$est[fit$res$taxon[ok], 1]    <- fit$res$lfc_disease[ok]
      s$stderr[fit$res$taxon[ok], 1] <- fit$res$se_disease[ok]
    }
    s
  })
  names(summary_stat_ANCOMBC2) <- studies

  ## ---- MaAsLin3 --------------------------------------------------------------
  message("[3/4] MaAsLin3 summary statistics")
  summary_stat_MaAsLin3 <- lapply(studies, function(d) {
    dat <- dataset_of(d)
    s   <- empty_sums(feature_universe, nrow(dat$Y))
    fit <- try(maaslin3::maaslin3(
      input_data = t(dat$Y), input_metadata = dat$meta,
      output = file.path(maaslin_out, d),
      formula = paste("~", paste(colnames(dat$meta), collapse = " + ")),
      normalization = "TSS", transform = "LOG", augment = TRUE,
      standardize = FALSE, max_significance = 0.1,
      median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
      max_pngs = 0, save_models = FALSE,
      plot_summary_plot = FALSE, plot_associations = FALSE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      r <- fit$fit_data_abundance$results %>% filter(!is.na(coef), !is.na(stderr))
      s$est[r$feature, 1]    <- r$coef
      s$stderr[r$feature, 1] <- r$stderr
    }
    s
  })
  names(summary_stat_MaAsLin3) <- studies

  ## ---- LinDA -----------------------------------------------------------------
  message("[4/4] LinDA summary statistics")
  summary_stat_LinDA <- lapply(studies, function(d) {
    dat <- dataset_of(d)
    s   <- empty_sums(feature_universe, nrow(dat$Y))
    fit <- try(MicrobiomeStat::linda(
      feature.dat = t(dat$Y), meta.dat = dat$meta,
      formula = paste("~", paste(colnames(dat$meta), collapse = " + ")),
      feature.dat.type = "count", prev.filter = 0,
      adaptive = TRUE, alpha = 0.05), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      o <- fit$output$disease
      s$est[rownames(o), 1]    <- o$log2FoldChange
      s$stderr[rownames(o), 1] <- o$lfcSE
    }
    s
  })
  names(summary_stat_LinDA) <- studies

  ## ---- context-level summaries and the mask ----------------------------------
  message("Meta-analysing within disease")
  meta_PALM     <- meta_by_disease(summary_stat_PALM,     n_by_study)
  meta_ANCOMBC2 <- meta_by_disease(summary_stat_ANCOMBC2, n_by_study)
  meta_MaAsLin3 <- meta_by_disease(summary_stat_MaAsLin3, n_by_study)
  meta_LinDA    <- meta_by_disease(summary_stat_LinDA,    n_by_study)

  ## The cutoff is defined on the reference (PALM) summaries and then applied to
  ## every method, so all six are fitted on the same features.
  if (analysis == "filtering") {
    feature_ID <- dense_features_of(meta_PALM, X)
    message("species estimable in > ", X, " diseases: ", length(feature_ID),
            " of ", length(unique(unlist(lapply(meta_PALM, function(s) rownames(s$est))))))
  } else {
    feature_ID <- dense_genera_of()
    message("dense genera: ", length(feature_ID), " of ", length(feature_universe))
  }

  summary_stat_meta_filter          <- mask_to(meta_PALM,     feature_ID)
  summary_stat_meta_filter_ancombc2 <- mask_to(meta_ANCOMBC2, feature_ID)
  summary_stat_meta_filter_MaAsLin3 <- mask_to(meta_MaAsLin3, feature_ID)
  summary_stat_meta_filter_Linda    <- mask_to(meta_LinDA,    feature_ID)

  ## Checkpoint before the model fits.  Delete this file to force a rerun.
  ## `X` is only meaningful for the filtering check, so it is only stored there -
  ## this keeps each checkpoint identical to the one its original script wrote.
  save(list = c("summary_stat_PALM", "summary_stat_ANCOMBC2", "summary_stat_MaAsLin3",
                "summary_stat_LinDA", "summary_stat_meta_filter",
                "summary_stat_meta_filter_ancombc2", "summary_stat_meta_filter_MaAsLin3",
                "summary_stat_meta_filter_Linda", "feature_ID", "disease_info",
                if (analysis == "filtering") "X"),
       file = ckpt_stats)

} else {
  message("resuming from ", basename(ckpt_stats))
  load(ckpt_stats)
}


# =============================================================================
#  3. Structure-learning meta-analysis at G = 2
# =============================================================================

  source("./utility/SMESH.R")

if (!file.exists(ckpt_fits)) {

  fit_smesh <- function(sums, tag) {
    message("SMESH fit: ", tag)
    smesh.meta.summary(summary.stats = sums, G = G, nperm = nperm_smesh,
                      doc = file.path(base_dir, paste0(cfg$doc_prefix, tag)),
                      verbose = TRUE)
  }

  SMESH_model    <- fit_smesh(summary_stat_meta_filter,          "PALM")
  ANCOMBC2_model <- fit_smesh(summary_stat_meta_filter_ancombc2, "ANCOMBC2")
  MaAsLin2_model <- fit_smesh(summary_stat_meta_filter_MaAsLin3, "MaAsLin3")
  LinDA_model    <- fit_smesh(summary_stat_meta_filter_Linda,    "LinDA")

  save(SMESH_model, ANCOMBC2_model, MaAsLin2_model, LinDA_model, file = ckpt_fits)

} else {
  message("resuming from ", basename(ckpt_fits))
  load(ckpt_fits)
}


# =============================================================================
#  4. SKM + FE and SHC + FE
# =============================================================================
#  Both cluster the contexts on their z-score profiles with a sparse method,
#  then fit a fixed-effect PALM model within each cluster.  K is fixed at G.

  source("./utility/PALM_tune.R")

  contexts <- names(summary_stat_meta_filter)
  zmat <- vapply(contexts, function(d) {
    s <- summary_stat_meta_filter[[d]]
    z <- s$est[feature_ID, 1] / s$stderr[feature_ID, 1]
    z[!is.finite(z)] <- 0
    z
  }, numeric(length(feature_ID)))
  Xmat <- t(zmat)                     # contexts x features
  storage.mode(Xmat) <- "double"

  ## Sparse k-means, with the L1 bound chosen by the package's permutation rule.
  w_skm  <- KMeansSparseCluster.permute(Xmat, K = G, nperms = 50, silent = TRUE)$bestw
  SKM_tab <- data.frame(
    cluster = as.integer(KMeansSparseCluster(Xmat, K = G, wbounds = w_skm)[[1]]$Cs),
    study   = contexts, row.names = contexts
  )

  ## Sparse hierarchical clustering, cut at G.
  w_shc  <- HierarchicalSparseCluster.permute(Xmat, nperms = 50)$bestw
  SHC_tab <- data.frame(
    cluster = as.integer(cutree(HierarchicalSparseCluster(Xmat, wbound = w_shc)$hc, k = G)),
    study   = contexts, row.names = contexts
  )

  ## Fixed-effect PALM within each cluster; one column of coefficients per
  ## cluster.  palm_tune() returns coefficients over whatever feature set the
  ## supplied summaries carry, so align by name rather than assuming a length.
  fe_signal <- function(tab) {
    cols <- lapply(sort(unique(tab$cluster)), function(g) {
      cf <- palm_tune(summary.stats = summary_stat_meta_filter[tab$study[tab$cluster == g]])$disease$coef
      if (is.null(names(cf))) names(cf) <- rownames(summary_stat_meta_filter[[1]]$est)
      cf
    })
    feats <- Reduce(union, lapply(cols, names))
    out <- vapply(cols, function(cf) {
      v <- setNames(rep(0, length(feats)), feats)
      v[names(cf)] <- cf
      v
    }, numeric(length(feats)))
    rownames(out) <- feats
    out
  }

  SKM_detect <- fe_signal(SKM_tab)
  SHC_detect <- fe_signal(SHC_tab)


# =============================================================================
#  5. Save
# =============================================================================

  save(list = c(
    ## per project-disease summary statistics
    "summary_stat_PALM", "summary_stat_ANCOMBC2", "summary_stat_MaAsLin3", "summary_stat_LinDA",
    ## context-level summaries, masked to the dense feature set
    "summary_stat_meta_filter", "summary_stat_meta_filter_ancombc2",
    "summary_stat_meta_filter_MaAsLin3", "summary_stat_meta_filter_Linda",
    ## fitted models at G = 2
    "SMESH_model", "ANCOMBC2_model", "MaAsLin2_model", "LinDA_model",
    "SKM_tab", "SHC_tab", "SKM_detect", "SHC_detect",
    ## bookkeeping
    "feature_ID", "disease_info", "G",
    if (analysis == "filtering") "X"),
    file = out_all
  )

  ## The taxonomy check also refreshes the two-object file that the existing
  ## Species_genus_consistency.R expects.
  if (!is.null(out_comp)) {
    save(summary_stat_meta_filter, SMESH_model, file = out_comp)
    message("saved: ", out_all, " and ", out_comp)
  } else {
    message("saved: ", out_all)
  }
