# =============================================================================
#   3_summary_statistics.R  --  per-context association summaries
# =============================================================================
#  The missing link between preprocessing and model fitting.  Takes the
#  processed data from 0_preprocessing.R, runs one differential-abundance
#  back-end per context, and writes the Summary_stat_*.Rdata files that
#  4_real_data_analysis.R consumes.
#
#  This is merged from the summary-statistic halves of the retired
#  Prepare_CRC.R, Prepare_HGMT.R and Prepare_GMrepo.R (now in Backup/prepare/).
#  Every method call is carried over verbatim; the per-application differences
#  live in the APPS table below.
#
#  Usage, from the project root:
#      Rscript Analysis/3_summary_statistics.R <application> <method>
#      Rscript Analysis/3_summary_statistics.R CRC PALM
#      Rscript Analysis/3_summary_statistics.R HGMT ANCOMBC2
#
#    application  GMrepo | HGMT | CRC
#    method       PALM | ANCOMBC2 | MaAsLin3 | LinDA | Melody
#
#  One invocation = one (application, method) cell.  Each is slow - MaAsLin3 and
#  ANCOMBC2 run a model per context - so they are meant to be run separately,
#  which is also how they were run originally.
#
#  Output:  <app>/Data/Summary_stat_<method>.Rdata     context-level (pooled)
#           <app>/Data/summary_stat.rds                 per-dataset, PALM
#           <app>/Data/summary_stat_<method>.rds        per-dataset, other back-ends
#
#  The per-dataset files are written only where a meta step happens (pan-disease
#  and pan-tumor).  Table S1 needs them: it reports the taxa each individual
#  study could estimate, which the pooled object no longer carries.
#
#  -------------------------------------------------------------------------
#  `correct = "tune"` is used for all three applications.  The retired scripts
#  were not consistent on this - Prepare_CRC.R used "median" and the surviving
#  pan-disease block used NULL - so "tune" is applied uniformly here, matching
#  Prepare_HGMT.R and the pan-disease sensitivity code (Analysis5.R/Analysis6.R).
#
#  Because of that, and because Prepare_GMrepo.R never wrote
#  Summary_stat_SMESH.Rdata at all, re-running this script will OVERWRITE the
#  Summary_stat_*.Rdata files that every published result was fitted on, and it
#  may not reproduce them exactly.  Back them up and compare before adopting the
#  new versions.
# =============================================================================

  rm(list = ls())

  library(dplyr)
  library(tibble)
  library(stringr)

  ## ---- per-application settings ----------------------------------------------
  ## `context_of` maps a per-dataset key to the context the paper reports.  For
  ## colorectal the strata ARE the contexts, so there is no meta-analysis step;
  ## for the other two, PALM::palm.meta.summary() pools the datasets of one
  ## disease / tumour type into a single context-level summary.
  APPS <- list(
    GMrepo = list(
      data_dir  = "./GMrepo_analysis/Data",
      source    = "rds",                       # rel.abd.rds + covariate.* + cluster
      min_ctx   = 8,        # estimable in MORE than 8 diseases, i.e. >= 9 (-> 244)
      mask_post = TRUE,     # the feature set is decided AFTER the meta step
      prev_cut  = 0.1,
      null_prev = 0,
      correct   = "tune",
      use_clust = TRUE,     # repeated-measure grouping passed to PALM
      meta      = TRUE,     # pool datasets within disease
      obj       = c(PALM = "summary_stat_meta_filter",
                    ANCOMBC2 = "summary_stat_meta_filter_ancombc2",
                    MaAsLin3 = "summary_stat_meta_filter_MaAsLin3",
                    LinDA    = "summary_stat_meta_filter_Linda",
                    Melody   = "summary_stat_melody")),
    HGMT = list(
      data_dir  = "./HGMT_analysis/Data",
      source    = "HGMT_analysis.Rdata",       # otu_filter, meta_filter, feature_ID
      min_ctx   = 6,        # recorded in HGMT_analysis.Rdata as feature_ID (159)
      mask_post = FALSE,
      prev_cut  = 0.1,
      null_prev = 0.1,
      correct   = "tune",
      use_clust = FALSE,
      meta      = TRUE,
      obj       = c(PALM = "summary_stat_meta_filter",
                    ANCOMBC2 = "summary_stat_meta_filter_ancombc2",
                    MaAsLin3 = "summary_stat_meta_filter_MaAsLin3",
                    LinDA    = "summary_stat_meta_filter_Linda",
                    Melody   = "summary_stat_melody")),
    CRC = list(
      data_dir  = "./CRC_analysis/Data",
      source    = "CRC_analysis.Rdata",        # otu_filter, meta_filter, depth
      min_ctx   = NA,       # feature set already fixed by 0_preprocessing.R (260)
      mask_post = FALSE,
      prev_cut  = 0.1,
      null_prev = 0,
      correct   = "tune",
      use_clust = FALSE,
      meta      = FALSE,    # strata are the contexts
      obj       = c(PALM = "summary_stat_filter",
                    ANCOMBC2 = "summary_stat_filter_ancombc2",
                    MaAsLin3 = "summary_stat_filter_MaAsLin3",
                    LinDA    = "summary_stat_filter_LinDA",
                    Melody   = "summary_stat_melody"))
  )

  OUT_FILE <- c(PALM = "Summary_stat_SMESH.Rdata", ANCOMBC2 = "Summary_stat_ancombc2.Rdata",
                MaAsLin3 = "Summary_stat_MaAsLin3.Rdata", LinDA = "Summary_stat_LinDA.Rdata",
                Melody = "Summary_stat_Melody.Rdata")

  ## ---- inputs -----------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2)
    stop("need <application> <method>, e.g.\n",
         "  Rscript Analysis/3_summary_statistics.R CRC PALM")
  app <- args[1]; method <- args[2]
  if (!app %in% names(APPS))
    stop("unknown application '", app, "'; choose ", paste(names(APPS), collapse = ", "))
  if (!method %in% names(OUT_FILE))
    stop("unknown method '", method, "'; choose ", paste(names(OUT_FILE), collapse = ", "))
  if (!dir.exists("utility"))
    stop("run this from the project root, e.g. Rscript Analysis/3_summary_statistics.R CRC PALM")

  cfg <- APPS[[app]]
  out_file <- file.path(cfg$data_dir, OUT_FILE[[method]])
  maaslin_out <- file.path(cfg$data_dir, "MaAsLin3")

  message(app, " | ", method, " -> ", out_file)


# =============================================================================
#  1. Load the processed data and build the analysis feature set
# =============================================================================
#  Each application arrives in a different shape, so this section normalises all
#  three to the same four objects:
#    otu_filter    per dataset count matrix, samples x features
#    meta_filter   per dataset metadata, aligned to otu_filter
#    cov_adjust    covariates to adjust for (NULL where none)
#    context_of    dataset key -> context label

  if (app == "CRC") {

    e <- new.env(); load(file.path(cfg$data_dir, cfg$source), envir = e)
    otu_filter  <- e$otu_filter
    meta_filter <- e$meta_filter
    cov_adjust  <- NULL
    ## disease indicator: anything that is not a Control
    grp_of     <- function(d) as.numeric(d$Disease != "Control")
    context_of <- setNames(names(otu_filter), names(otu_filter))   # identity
    feature_ID <- unique(unlist(lapply(otu_filter, colnames)))

  } else if (app == "HGMT") {

    ## The filtered object 0_preprocessing.R writes: already reduced to the 159
    ## genera the analysis uses, so the feature set is read, not re-derived.
    e <- new.env(); load(file.path(cfg$data_dir, cfg$source), envir = e)
    otu_filter <- e$otu_filter
    feature_ID <- e$feature_ID
    cov_adjust <- e$covariate.adjust_filter

    ## meta_filter holds tibbles, which carry no rownames; index them by Run ID
    ## so they line up with otu_filter.
    meta_filter <- lapply(names(otu_filter), function(d) {
      m <- as.data.frame(e$meta_filter[[d]], stringsAsFactors = FALSE)
      rownames(m) <- as.character(m$`Run ID`)
      m[rownames(otu_filter[[d]]), , drop = FALSE]
    })
    names(meta_filter) <- names(otu_filter)

    ## context = tumour type, recovered by stripping the project prefix
    context_of <- setNames(
      str_remove(names(otu_filter), "^PRJ[A-Z0-9-]+(?:_(?:16S|WGS|PE|SE))?_"),
      names(otu_filter))

    grp_of <- function(d) as.numeric(d$`Phenotype name` != "Health")

  } else {   ## GMrepo

    rel.abd    <- readRDS(file.path(cfg$data_dir, "rel.abd.rds"))
    cov.int    <- readRDS(file.path(cfg$data_dir, "covariate.interest.rds"))
    cov_all    <- readRDS(file.path(cfg$data_dir, "covariate.adjust.rds"))
    clust      <- readRDS(file.path(cfg$data_dir, "cluster.rds"))
    psum       <- openxlsx::read.xlsx(file.path(cfg$data_dir, "GMrepo_project_summary.xlsx"))

    ## context = disease name, looked up from the project-disease key
    mesh <- unique(psum[, c("disease_mesh", "disease_name")])
    context_of <- setNames(
      mesh$disease_name[match(sub("^.*_", "", names(rel.abd)), mesh$disease_mesh)],
      names(rel.abd))
    stopifnot(!anyNA(context_of))

    ## Unlike the other two, the pan-disease feature set is decided AFTER the
    ## meta-analysis (see section 4), so every species enters the back-end here.
    feature_ID <- sort(unique(unlist(lapply(rel.abd, colnames))))
    otu_filter <- rel.abd
    meta_filter <- list(); cov_adjust <- list()
    for (d in names(otu_filter)) {
      keep <- rownames(otu_filter[[d]])
      meta_filter[[d]] <- data.frame(disease = as.integer(cov.int[[d]][keep, 1]),
                                     row.names = keep)
      if (!is.null(cov_all[[d]]))
        cov_adjust[[d]] <- cov_all[[d]][keep, , drop = FALSE]
    }
    grp_of <- function(d) d$disease
  }

  message("  contexts: ", length(unique(context_of)),
          " | datasets: ", length(otu_filter),
          " | features: ", length(feature_ID))


# =============================================================================
#  2. Helpers
# =============================================================================

  ## An empty (est, stderr) pair on the analysis feature set, so every dataset
  ## returns the same shape whatever it could fit.
  empty_sums <- function(n) {
    m <- matrix(NA_real_, length(feature_ID), 1,
                dimnames = list(feature_ID, "disease"))
    list(est = m, stderr = m, n = n)
  }

  ## The per-dataset matrix each non-PALM back-end is given: the analysis feature
  ## set, restricted to features this dataset can actually fit.
  dataset_of <- function(d) {
    Y <- otu_filter[[d]]
    keep <- intersect(feature_ID,
                      colnames(Y)[colMeans(Y != 0) >= cfg$prev_cut & apply(Y, 2, var) > 0])
    Y[, keep, drop = FALSE]
  }

  ## Pool the per-dataset summaries of one context with a fixed-effect PALM
  ## meta-analysis.  Applications whose datasets already are contexts skip this.
  meta_by_context <- function(sums) {
    if (!cfg$meta) return(sums)
    out <- list()
    for (ctx in unique(context_of)) {
      keep <- names(sums)[context_of[names(sums)] == ctx]
      res  <- PALM::palm.meta.summary(summary.stats = sums[keep])
      out[[ctx]] <- list(
        est    = matrix(res$disease$coef, ncol = 1,
                        dimnames = list(res$disease$feature, "disease")),
        stderr = matrix(res$disease$stderr, ncol = 1,
                        dimnames = list(res$disease$feature, "disease")),
        n      = sum(vapply(sums[keep], function(d) d$n, numeric(1))))
    }
    out
  }

  ## Pan-disease only: a species is kept if it is estimable in MORE than
  ## `min_ctx` diseases, counted on the UNMASKED context-level summaries.  The
  ## set is defined once, on PALM, and every other back-end is masked to it, so
  ## all methods are fitted on the same species.  PALM must therefore be run
  ## first; the set is stored alongside its output.
  dense_features_of <- function(meta_sums, X) {
    feats <- unique(unlist(lapply(meta_sums, function(s) rownames(s$est))))
    n_dis <- vapply(feats, function(f)
      sum(vapply(meta_sums, function(s)
        f %in% rownames(s$est) && !is.na(s$est[f, "disease"]), logical(1))), numeric(1))
    names(n_dis)[n_dis > X]
  }

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

  ## Applied to every back-end's context-level summaries.  A no-op unless the
  ## application decides its feature set after the meta step.
  finalize <- function(sums) {
    if (!isTRUE(cfg$mask_post)) return(sums)
    ref <- file.path(cfg$data_dir, OUT_FILE[["PALM"]])
    if (method == "PALM") {
      ## Prefer the recorded set: species_lst.Rdata holds the `dense_feature`
      ## the published analysis used.  Re-deriving it is the documented fallback
      ## for a fresh dataset, but on this one it must not be allowed to drift.
      rec <- file.path(cfg$data_dir, "species_lst.Rdata")
      if (file.exists(rec)) {
        er <- new.env(); load(rec, envir = er)
        feats <- er$dense_feature
        message("  dense feature set: ", length(feats), " (recorded in species_lst.Rdata)")
      } else {
        feats <- dense_features_of(sums, cfg$min_ctx)
        message("  dense feature set: ", length(feats), " of ",
                length(unique(unlist(lapply(sums, function(s) rownames(s$est))))),
                "  (re-derived; species_lst.Rdata not found)")
      }
      assign("dense_feature_ID", feats, envir = globalenv())
      return(mask_to(sums, feats))
    }
    if (!file.exists(ref))
      stop("run PALM first - ", basename(ref), " defines the feature set every\n",
           "  other method is masked to:  Rscript Analysis/3_summary_statistics.R ",
           app, " PALM")
    e <- new.env(); load(ref, envir = e)
    if (!"dense_feature_ID" %in% ls(e))
      stop(basename(ref), " has no `dense_feature_ID`; regenerate it with this script")
    mask_to(sums, get("dense_feature_ID", envir = e))
  }

  ## The per-dataset summaries, before pooling.  Table S1 reports how many taxa
  ## each individual study could estimate, which the pooled object cannot answer.
  ## Only meaningful where a meta step happens - for colorectal the strata ARE
  ## the contexts, so the pooled file already is the per-study object.
  save_perstudy <- function(sums) {
    if (!isTRUE(cfg$meta)) return(invisible(NULL))
    f <- file.path(cfg$data_dir,
                   if (method == "PALM") "summary_stat.rds"
                   else sprintf("summary_stat_%s.rds", method))
    saveRDS(sums, f)
    message("  saved: ", f, "  (", length(sums), " datasets, pre-meta)")
  }

  save_as <- function(obj) {
    assign(cfg$obj[[method]], obj)
    keep <- cfg$obj[[method]]
    if (isTRUE(cfg$mask_post) && method == "PALM") keep <- c(keep, "dense_feature_ID")
    save(list = keep, file = out_file)
    message("  saved: ", out_file, "  (", length(obj), " contexts, object `",
            cfg$obj[[method]], "`)")
  }


# =============================================================================
#  3. The back-ends
# =============================================================================

if (method == "PALM") {

  ## PALM fits all datasets jointly: a null model, then one summary per dataset.
  null_obj <- PALM::palm.null.model(
    rel.abd          = otu_filter,
    covariate.adjust = if (is.null(cov_adjust) || !length(cov_adjust)) NULL else cov_adjust,
    prev.filter      = cfg$null_prev)

  ci <- lapply(meta_filter, function(d)
    matrix(grp_of(d), ncol = 1, dimnames = list(rownames(d), "disease")))

  sums <- if (cfg$use_clust)
    PALM::palm.get.summary(null.obj = null_obj, covariate.interest = ci,
                           cluster = readRDS(file.path(cfg$data_dir, "cluster.rds")),
                           correct = cfg$correct)
  else
    PALM::palm.get.summary(null.obj = null_obj, covariate.interest = ci,
                           correct = cfg$correct)

  save_perstudy(sums)
  save_as(finalize(meta_by_context(sums)))


} else if (method == "Melody") {

  ## Melody has its own null model and needs no meta step - it pools internally.
  null_obj <- miMeta::melody.null.model(
    rel.abd          = otu_filter,
    covariate.adjust = if (is.null(cov_adjust) || !length(cov_adjust)) NULL else cov_adjust,
    prev.filter      = cfg$null_prev)

  ci <- lapply(meta_filter, function(d)
    matrix(grp_of(d), ncol = 1, dimnames = list(rownames(d), "disease")))

  sums <- if (cfg$use_clust)
    miMeta::melody.get.summary(null.obj = null_obj, covariate.interest = ci,
                               cluster = readRDS(file.path(cfg$data_dir, "cluster.rds")))
  else
    miMeta::melody.get.summary(null.obj = null_obj, covariate.interest = ci)

  save_as(sums)


} else if (method == "ANCOMBC2") {

  source("./utility/ancombc.R")
  sums <- list()
  for (d in names(otu_filter)) {
    Y <- dataset_of(d)
    s <- empty_sums(ncol(Y))
    fit <- try(ancombc.fun(
      feature.table = t(Y),
      meta          = meta_filter[[d]] %>% dplyr::transmute(disease = grp_of(.)),
      formula       = "disease",
      adjust.method = "fdr",
      group         = NULL,
      subject       = NULL,
      method        = "ancombc2"), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      ok <- fit$res$passed_ss_disease
      s$est[fit$res$taxon[ok], 1]    <- fit$res$lfc_disease[ok]
      s$stderr[fit$res$taxon[ok], 1] <- fit$res$se_disease[ok]
    }
    sums[[d]] <- s
    message("    ", d, " (", which(names(otu_filter) == d), "/", length(otu_filter), ")")
  }
  save_perstudy(sums)
  save_as(finalize(meta_by_context(sums)))


} else if (method == "MaAsLin3") {

  dir.create(maaslin_out, showWarnings = FALSE, recursive = TRUE)
  sums <- list()
  for (d in names(otu_filter)) {
    Y <- dataset_of(d)
    s <- empty_sums(ncol(Y))
    fit <- try(maaslin3::maaslin3(
      input_data     = Y,
      input_metadata = meta_filter[[d]] %>% dplyr::transmute(disease = grp_of(.)),
      output         = file.path(maaslin_out, d),
      formula        = "~ disease",
      normalization  = "TSS",
      transform      = "LOG",
      augment        = TRUE,
      standardize    = FALSE,
      max_significance = 0.1,
      median_comparison_abundance  = TRUE,
      median_comparison_prevalence = FALSE,
      max_pngs    = 100,
      save_models = FALSE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      r <- fit$fit_data_abundance$results %>% dplyr::filter(!is.na(coef), !is.na(stderr))
      s$est[r$feature, 1]    <- r$coef
      s$stderr[r$feature, 1] <- r$stderr
    }
    sums[[d]] <- s
    message("    ", d, " (", which(names(otu_filter) == d), "/", length(otu_filter), ")")
  }
  save_perstudy(sums)
  save_as(finalize(meta_by_context(sums)))


} else {   ## LinDA

  sums <- list()
  for (d in names(otu_filter)) {
    Y <- t(dataset_of(d))
    Y <- Y[, colSums(Y) > 0, drop = FALSE]
    s <- empty_sums(ncol(Y))
    fit <- try(MicrobiomeStat::linda(
      feature.dat      = Y,
      meta.dat         = meta_filter[[d]] %>% dplyr::transmute(disease = grp_of(.)),
      formula          = "~ disease",
      feature.dat.type = "count",
      prev.filter      = 0,
      adaptive         = TRUE,
      alpha            = 0.05), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      o <- fit$output$disease
      s$est[rownames(o), 1]    <- o$log2FoldChange
      s$stderr[rownames(o), 1] <- o$lfcSE
    }
    sums[[d]] <- s
    message("    ", d, " (", which(names(otu_filter) == d), "/", length(otu_filter), ")")
  }
  save_perstudy(sums)
  save_as(finalize(meta_by_context(sums)))
}

  message("done.")
