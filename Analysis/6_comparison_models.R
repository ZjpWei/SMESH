# =============================================================================
#   6_comparison_models.R  --  comparison-method models at the selected G
# =============================================================================
#  Run this AFTER 4_real_data_analysis.R, once the number of clusters G has been
#  chosen for each application.  It produces the three comparison-method objects
#  the figure scripts load from <app>/Model/, for all three applications in one
#  invocation:
#
#    SKM + FE   sparse k-means clustering of the contexts, then a fixed-effect
#               PALM model within each cluster       -> SKMean_FE_G<G>.Rdata
#    SHC + FE   sparse hierarchical clustering, same second step
#                                                    -> SHC_FE_G<G>.Rdata
#    Melody     fixed-effect meta-analysis benchmark -> Melody_model.Rdata
#
#  SKM and SHC are two-step methods: 4_real_data_analysis.R already did the
#  clustering half and saved it as `tab` in <app>_loso/Model<G>_{SKM,SHC}_s100.Rdata
#  (s = 100 is the fit on every context).  This script does the second half -
#  the within-cluster fixed-effect fit - which is what the figures actually plot.
#
#  Usage, from the project root:
#      Rscript Analysis/6_comparison_models.R              # all three applications
#      Rscript Analysis/6_comparison_models.R HGMT         # just one
#      Rscript Analysis/6_comparison_models.R HGMT 3       # override the G to use
#
#  Output objects keep the names the figure scripts expect:
#      SKMean_FE_G<G>.Rdata / SHC_FE_G<G>.Rdata  ->  `detect.signal`
#                                                    (features x G coefficients)
#      Melody_model.Rdata                        ->  `Melody_mod`
# =============================================================================

  rm(list = ls())

  library(dplyr)
  library(tibble)
  ## utility/PALM_tune.R calls abess() and glmtlp() unqualified.
  library(abess)
  library(glmtlp)

  set.seed(2026)

  ## ---- per-application settings ----------------------------------------------
  ## `G` is the number of clusters selected for that application; change it here
  ## if the selected model changes, or pass a second command-line argument.
  APPS <- list(
    GMrepo = list(
      G          = 2,
      data_dir   = "./GMrepo_analysis/Data",
      loso_dir   = "./GMrepo_analysis/GMrepo_loso",
      model_dir  = "./GMrepo_analysis/Model",
      smesh_obj  = "summary_stat_meta_filter"),
    HGMT = list(
      G          = 4,
      data_dir   = "./HGMT_analysis/Data",
      loso_dir   = "./HGMT_analysis/HGMT_loso",
      model_dir  = "./HGMT_analysis/Model",
      smesh_obj  = "summary_stat_meta_filter"),
    CRC = list(
      G          = 2,
      data_dir   = "./CRC_analysis/Data",
      loso_dir   = "./CRC_analysis/CRC_loso",
      model_dir  = "./CRC_analysis/Model",
      smesh_obj  = "summary_stat_filter")
  )

  ## ---- inputs -----------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  apps <- if (length(args) >= 1) args[1] else names(APPS)
  if (!all(apps %in% names(APPS)))
    stop("unknown application '", setdiff(apps, names(APPS))[1], "'; choose ",
         paste(names(APPS), collapse = ", "))
  G_override <- if (length(args) >= 2) as.numeric(args[2]) else NA
  if (length(args) >= 2 && length(apps) > 1)
    stop("a G override applies to one application; name it, e.g.\n",
         "  Rscript Analysis/6_comparison_models.R HGMT 3")

  if (!dir.exists("utility"))
    stop("run this from the project root, e.g. Rscript Analysis/6_comparison_models.R")

  source("./utility/PALM_tune.R")


# =============================================================================
#  Helpers
# =============================================================================

  #' Load a single named object out of an .Rdata file.
  load_obj <- function(file, name) {
    if (!file.exists(file)) return(NULL)
    e <- new.env(); load(file, envir = e)
    if (!name %in% ls(e)) return(NULL)
    get(name, envir = e)
  }

  #' Fixed-effect PALM within each cluster; one column of coefficients per
  #' cluster.  palm_tune() returns coefficients over whatever feature set the
  #' supplied summaries carry, so align by name rather than assuming a length.
  #' This is the same construction 5_gmrepo_sensitivity.R uses.
  fe_signal <- function(tab, sums) {
    cols <- lapply(sort(unique(tab$cluster)), function(g) {
      cf <- palm_tune(summary.stats = sums[tab$study[tab$cluster == g]])$disease$coef
      if (is.null(names(cf))) names(cf) <- rownames(sums[[1]]$est)
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


# =============================================================================
#  Per application
# =============================================================================

for (app in apps) {

  cfg <- APPS[[app]]
  G   <- if (is.na(G_override)) cfg$G else G_override
  dir.create(cfg$model_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n== ", app, " | G = ", G, " ==")

  ## Context-level summaries: the same object SMESH was fitted on, so the
  ## comparison methods see exactly the same data.
  smesh_file <- file.path(cfg$data_dir, "Summary_stat_SMESH.Rdata")
  sums <- load_obj(smesh_file, cfg$smesh_obj)
  if (is.null(sums))
    stop("missing ", cfg$smesh_obj, " in ", smesh_file,
         " - run Analysis/4_real_data_analysis.R first")

  ## ---- SKM + FE and SHC + FE -------------------------------------------------
  ## The clustering half comes from 4_real_data_analysis.R (`tab`, at s = 100);
  ## only the within-cluster fixed-effect fit happens here.
  for (m in c("SKM", "SHC")) {
    clust_file <- file.path(cfg$loso_dir, sprintf("Model%d_%s_s100.Rdata", G, m))
    tab <- load_obj(clust_file, "tab")
    if (is.null(tab)) {
      warning(app, ": no clustering at G = ", G, " for ", m, " (", clust_file,
              ").  Run:  Rscript Analysis/4_real_data_analysis.R ", app, " ", G,
              " ", m, " 100")
      next
    }

    message("  ", m, " + FE: ", nrow(tab), " contexts, ",
            length(unique(tab$cluster)), " clusters")
    detect.signal <- fe_signal(tab, sums)

    out <- file.path(cfg$model_dir,
                     sprintf("%s_G%d.Rdata", if (m == "SKM") "SKMean_FE" else "SHC_FE", G))
    save(detect.signal, file = out)
    message("    saved: ", out, "  (", nrow(detect.signal), " x ",
            ncol(detect.signal), ", ", sum(rowSums(detect.signal != 0) > 0),
            " features selected)")
  }

  ## ---- Melody ----------------------------------------------------------------
  ## Melody is not a clustering method: it fits one fixed-effect meta-analysis
  ## across all contexts, so it needs no G.  It reads its own summary statistics,
  ## which are produced by miMeta::melody.get.summary() from the per-sample data.
  melody_file <- file.path(cfg$data_dir, "Summary_stat_Melody.Rdata")
  melody_sums <- load_obj(melody_file, "summary_stat_melody")
  if (is.null(melody_sums)) {
    message("  Melody: skipped - ", melody_file, " not found.")
    message("          Melody summary statistics are generated by ",
            "miMeta::melody.get.summary()")
    message("          from rel.abd / covariate.interest / covariate.adjust / cluster; ",
            "that step")
    message("          is not yet part of this folder.  ",
            "Model/Melody_model.Rdata is kept as is.")
  } else {
    message("  Melody: ", length(melody_sums), " contexts")
    Melody_mod <- miMeta::melody.meta.summary(summary.stats = melody_sums,
                                              output.best.one = TRUE)
    out <- file.path(cfg$model_dir, "Melody_model.Rdata")
    save(Melody_mod, file = out)
    message("    saved: ", out)
  }
}

  message("\ndone.")
