# =============================================================================
#   4_real_data_analysis.R  --  fit one model to one real data application
# =============================================================================
#  Takes the context-level association summaries produced upstream and fits one
#  clustering model to them, either on every context or with one context held
#  out.  One invocation = one (application, G, method, s) cell; the full set of
#  results is produced by looping over them, which is how they were run on the
#  cluster.
#
#  This merges CRC_meta_loso.R, HGMT_meta_loso.R and GMrepo_meta_loso.R.  Those
#  three were structurally identical but each carried its own settings, so the
#  differences are held in the APPS table below and nothing else varies.  No
#  method parameter has been changed.
#
#  Usage, from the project root:
#      Rscript Analysis/4_real_data_analysis.R <application> <G> <method> <s>
#      Rscript Analysis/4_real_data_analysis.R HGMT 4 SMESH 100
#
#    application  GMrepo | HGMT | CRC
#    G            number of clusters to fit
#    method       SMESH | ANCOMBC2 | MaAsLin3 | LinDA | SKM | SHC
#    s            context to hold out, 1..L for leave-one-context-out.
#                 USE s = 100 TO FIT EVERY CONTEXT: the code subsets with
#                 `summary.stats[-s]`, and a negative index past the end of the
#                 list drops nothing, so s = 100 keeps all of them.  That is the
#                 convention the saved file names use (`..._s100.Rdata`).
#
#  Output:  <app>_loso/Model<G>_<method>_s<s>.Rdata
#           SMESH / ANCOMBC2 / MaAsLin3 / LinDA save a fitted model object;
#           SKM and SHC save `tab`, a cluster/study data frame.
# =============================================================================

  rm(list = ls())

  library("phyloseq")
  library("coin")
  library("purrr")
  library("tidyverse")
  library("abess")
  library("dplyr")
  library("ggplot2")
  library("stringr")
  library("UpSetR")
  library("sparcl")
  library("MicrobiomeStat")
  library("ANCOMBC")
  library("MIDASim")
  library("PALM")
  library("cluster")
  library("mclust")

  ## The SMESH implementation, shared by all three applications.
  ## NOTE: the saved results in the *_loso/ folders were produced by earlier
  ## implementations (see git history / Backup/utility/), so re-running will not
  ## reproduce those files bit for bit.
  SMESH_SRC <- "./utility/SMESH.R"

  ## ---- per-application settings ----------------------------------------------
  ## These are the only differences between the three original scripts.  The
  ## permutation counts are NOT uniform across applications - they are carried
  ## over exactly as each application was run.
  APPS <- list(
    GMrepo = list(
      data_dir  = "./GMrepo_analysis/Data",
      out_dir   = "./GMrepo_analysis/GMrepo_loso",
      objects   = c(SMESH    = "summary_stat_meta_filter",
                    ANCOMBC2 = "summary_stat_meta_filter_ancombc2",
                    MaAsLin3 = "summary_stat_meta_filter_MaAsLin3",
                    LinDA    = "summary_stat_meta_filter_Linda"),
      nperm     = c(SMESH = 20, ANCOMBC2 = 20, MaAsLin3 = 20, LinDA = 20),
      nperms_sparse = 50),
    HGMT = list(
      data_dir  = "./HGMT_analysis/Data",
      out_dir   = "./HGMT_analysis/HGMT_loso",
      objects   = c(SMESH    = "summary_stat_meta_filter",
                    ANCOMBC2 = "summary_stat_meta_filter_ancombc2",
                    MaAsLin3 = "summary_stat_meta_filter_MaAsLin3",
                    LinDA    = "summary_stat_meta_filter_Linda"),
      nperm     = c(SMESH = 50, ANCOMBC2 = 50, MaAsLin3 = 50, LinDA = 50),
      nperms_sparse = 30),
    CRC = list(
      data_dir  = "./CRC_analysis/Data",
      out_dir   = "./CRC_analysis/CRC_loso",
      objects   = c(SMESH    = "summary_stat_filter",
                    ANCOMBC2 = "summary_stat_filter_ancombc2",
                    MaAsLin3 = "summary_stat_filter_MaAsLin3",
                    LinDA    = "summary_stat_filter_LinDA"),
      nperm     = c(SMESH = 50, ANCOMBC2 = 50, MaAsLin3 = 50, LinDA = 50),
      nperms_sparse = 50)
  )

  ## The summary statistics of each method live in their own file, and the model
  ## object each method saves keeps the name the figure scripts expect.
  SUMMARY_FILE <- c(SMESH    = "Summary_stat_SMESH.Rdata",
                    ANCOMBC2 = "Summary_stat_ancombc2.Rdata",
                    MaAsLin3 = "Summary_stat_MaAsLin3.Rdata",
                    LinDA    = "Summary_stat_LinDA.Rdata")
  MODEL_NAME   <- c(SMESH    = "SMESH_model",
                    ANCOMBC2 = "ANCOMBC2_model",
                    MaAsLin3 = "MaAsLin2_model",   # legacy name, kept for the saved files
                    LinDA    = "LinDA_model")

  ## ---- inputs -----------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4)
    stop("need <application> <G> <method> <s>, e.g.\n",
         "  Rscript Analysis/4_real_data_analysis.R HGMT 4 SMESH 100")
  print(args)

  app    <- args[1]
  G      <- as.numeric(args[2])
  method <- args[3]
  s      <- as.numeric(args[4])

  if (!app %in% names(APPS))
    stop("unknown application '", app, "'; choose ", paste(names(APPS), collapse = ", "))
  if (!method %in% c(names(SUMMARY_FILE), "SKM", "SHC"))
    stop("unknown method '", method, "'; choose ",
         paste(c(names(SUMMARY_FILE), "SKM", "SHC"), collapse = ", "))

  cfg <- APPS[[app]]
  if (!dir.exists("utility"))
    stop("run this from the project root, e.g. Rscript Analysis/4_real_data_analysis.R ...")
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(cfg$out_dir, paste0("Model", G, "_", method, "_s", s, ".Rdata"))

  #' Load one method's summary statistics and return them under a common name.
  #' SKM and SHC cluster the PALM z-scores, so they read the SMESH file too.
  load_summaries <- function(which_method) {
    f <- file.path(cfg$data_dir, SUMMARY_FILE[[which_method]])
    if (!file.exists(f)) stop("missing summary statistics: ", f)
    e <- new.env(); load(f, envir = e)
    get(cfg$objects[[which_method]], envir = e)
  }

  message(app, " | G = ", G, " | method = ", method, " | s = ", s,
          if (s >= 100) "  (all contexts)" else "  (leave-one-context-out)")


# =============================================================================
#  SMESH and the three alternative input summaries
# =============================================================================
#  Identical call in every case; only the summaries, the permutation count and
#  the name the object is saved under differ.

  if (method %in% names(SUMMARY_FILE)) {

    source(SMESH_SRC)

    summary_stats <- load_summaries(method)

    fit <- smesh.meta.summary(summary.stats = summary_stats[-s],
                             G       = G,
                             nperm   = cfg$nperm[[method]],
                             verbose = TRUE)

    ## Save under the object name the downstream figure scripts expect.
    assign(MODEL_NAME[[method]], fit)
    save(list = MODEL_NAME[[method]], file = out_file)


# =============================================================================
#  SKM + FE and SHC + FE
# =============================================================================
#  Two-step alternatives: cluster the contexts on their PALM association
#  z-scores first, then a fixed-effect meta-analysis runs within each cluster
#  (that second step lives in the figure pipeline).  Only the clustering is done
#  here, and only `tab` is saved.

  } else {

    summary_stats <- load_summaries("SMESH")

    ## beta_summary: taxa x contexts (z-scores).  sparcl cannot take missing
    ## values, so a feature a context could not estimate is entered as 0.
    feature_ID <- unique(unlist(sapply(summary_stats[-s],
                                       function(d) rownames(d$est)[!is.na(d$est)])))
    beta_summary <- matrix(NA, nrow = length(feature_ID), ncol = length(summary_stats[-s]),
                           dimnames = list(feature_ID, names(summary_stats[-s])))
    for (d in names(summary_stats)[-s]) {
      beta_summary[feature_ID, d] <-
        summary_stats[[d]]$est[feature_ID, 1] / summary_stats[[d]]$stderr[feature_ID, 1]
    }
    beta_summary[is.na(beta_summary)] <- 0
    X <- t(beta_summary)   # contexts x taxa

    set.seed(123)
    Kmax_user <- 8
    Kmax_safe <- min(Kmax_user, floor(nrow(X) / 2))
    Kmax_safe <- max(2, Kmax_safe)

    if (method == "SKM") {
      K_for_w <- min(4, Kmax_safe)

      ## The permutation step only picks the L1 bound; the clustering itself is
      ## done at the requested G.
      perm_skm <- KMeansSparseCluster.permute(
        X,
        K      = G,
        nperms = cfg$nperms_sparse,
        silent = TRUE
      )
      w_skm <- perm_skm$bestw

      if (G <= 1) {
        fit_skm <- NULL
        cl_skm  <- rep(1L, nrow(X))
        ws_skm  <- NULL
      } else {
        fit_skm <- KMeansSparseCluster(X, K = G, wbounds = w_skm)
        cl_skm  <- fit_skm[[1]]$Cs
        ws_skm  <- fit_skm[[1]]$ws
      }
      cl <- cl_skm

    } else {   ## SHC
      perm_shc <- HierarchicalSparseCluster.permute(X, nperms = cfg$nperms_sparse)
      w_shc    <- perm_shc$bestw
      fit_shc  <- HierarchicalSparseCluster(X, wbound = w_shc)

      if (G == 1) {
        cl_shc <- rep(1L, nrow(X))
        ws_shc <- NULL
      } else {
        cl_shc <- cutree(fit_shc$hc, k = G)
        ws_shc <- fit_shc$ws
      }
      names(cl_shc) <- colnames(beta_summary)
      cl <- cl_shc
    }

    tab <- tibble::tibble(cluster = cl, study = names(cl))
    save(tab, file = out_file)
  }

  message("written: ", out_file)
