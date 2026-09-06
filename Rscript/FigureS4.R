# =============================================================================
#   Figure S4 - pan-tumor application, sensitivity analyses
#
#  Run from the project root: Rscript Rscript/FigureS4.R
# =============================================================================
#  Panel A  cluster-level effects, SMESH-PALM stacked against every competitor
#  Panel B  number of clusters selected, all data vs leave-one-context-out
#  Panel C  retention, stability and ARI under leave-one-context-out
#  Panel D  cluster assignments across methods
#  Panel E  signature counts and Jaccard overlap across methods
#
#  Display order matches the main figure (see Rscript/Figure5.R): clusters run
#  most parsimonious first, studies are ordered inside a cluster by
#  hierarchical clustering of their summary statistics.
# =============================================================================
  
  rm(list = ls())
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(clue)
  
  source("./utility/heatmap_util.R")
  
  ## ---- configuration ---------------------------------------------------------
  loso_dir  <- "./HGMT_analysis/HGMT_loso"
  model_dir <- "./HGMT_analysis/Model"
  data_dir  <- "./HGMT_analysis/Data"
  fig_dir   <- "./HGMT_analysis/Figure"
  
  n_studies <- 15        # leave-one-context-out runs; "s100" = all contexts
  G_grid    <- 2:6       # candidate cluster counts
  cut_off   <- 0.015      # relative GIC gain below which we stop adding clusters
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## Each method's fitting script saves its model under a different object name.
  model_object <- c(
    SMESH    = "SMESH_model",
    ANCOMBC2 = "ANCOMBC2_model",
    LinDA    = "LinDA_model",
    MaAsLin3 = "MaAsLin2_model"   # legacy object name, kept for the saved files
  )
  
  ## Label used in the figures -> object name inside the .Rdata files.
  gic_methods <- c(
    "SMESH-PALM"     = "SMESH",
    "SMESH-ANCOMBC2" = "ANCOMBC2",
    "SMESH-LinDA"    = "LinDA",
    "SMESH-MaAsLin3" = "MaAsLin3"
  )
  
  #' Load one fit: `Model<G>_<method>_s<s>.Rdata`.
  load_fit <- function(method, G, s) {
    f <- file.path(loso_dir, paste0("Model", G, "_", method, "_s", s, ".Rdata"))
    e <- new.env()
    load(f, envir = e)
    get(model_object[[method]], envir = e)
  }
  
  #' Load the `tab` (cluster / study) saved by the SKM and SHC scripts.
  load_tab <- function(method, G, s) {
    f <- file.path(loso_dir, paste0("Model", G, "_", method, "_s", s, ".Rdata"))
    e <- new.env()
    load(f, envir = e)
    setNames(e$tab$cluster, e$tab$study)
  }
  
  #' Hard cluster assignment, named by study.
  cluster_vec <- function(model) {
    setNames(apply(model$disease$W, 1, which.max), rownames(model$disease$W))
  }
  
  #' Smallest G whose relative GIC gain over the next model drops below cut_off.
  #'
  #' Falls back to the largest candidate when the gain never flattens out.
  select_G <- function(method, s, G_grid. = G_grid, cut_off. = cut_off) {
    gic <- vapply(G_grid., function(g) load_fit(method, g, s)$disease$gic, numeric(1))
    rel <- (head(gic, -1) - tail(gic, -1)) / head(gic, -1)
    k <- which(rel <= cut_off.)[1]
    if (is.na(k)) max(G_grid.) else G_grid.[k]
  }
  
  
  # =============================================================================
  #  Panel A - cluster-level effects, method by method
  # =============================================================================
  
  load(file.path(data_dir, "Summary_stat_SMESH.Rdata"))   # summary_stat_meta_filter
  
  ## The panels compare methods at a common G, the one reported in the main text.
  G <- 4
  
  SMESH_model    <- load_fit("SMESH",    G, 100)
  ANCOMBC2_model <- load_fit("ANCOMBC2", G, 100)
  LinDA_model    <- load_fit("LinDA",    G, 100)
  MaAsLin2_model <- load_fit("MaAsLin3", G, 100)
  
  SHC_tab <- load_tab("SHC", G, 100)
  SKM_tab <- load_tab("SKM", G, 100)
  
  load(file.path(model_dir, "SHC_FE_G4.Rdata"));    SHC_detect <- detect.signal
  load(file.path(model_dir, "SKMean_FE_G4.Rdata")); SKM_detect <- detect.signal
  load(file.path(model_dir, "Melody_model.Rdata"))  # Melody_mod
  
  W          <- round(SMESH_model$disease$W)
  feature.ID <- rownames(SMESH_model$disease$mu)
  ref_cluster <- cluster_vec(SMESH_model)
  ancombc2_cluster <- cluster_vec(ANCOMBC2_model)
  ancombc2_cluster[ancombc2_cluster == 2 | ancombc2_cluster == 4] <- 
    ancombc2_cluster[ancombc2_cluster == 2 | ancombc2_cluster == 4] %% 4 + 2
  
  stats <- summary_stat_matrices(summary_stat_meta_filter, feature.ID)
  beta  <- mask_unobserved_clusters(SMESH_model$disease$mu, W, stats$est)
  
  ## ---- display order (same rule as the main figure) --------------------------
  cluster_order <- order_clusters_by_selection(beta)
  ## Same feature subset as the main figure: only the selected signatures drive
  ## the within-cluster ordering, otherwise the two figures disagree.
  sel_rows      <- function(b) names(which(rowSums(b != 0) > 0))
  study_order   <- order_studies_within_clusters(stats$est[sel_rows(beta), ],
                                                 W, cluster_order)
  species_lst   <- build_signature_groups(beta, cluster_order)
  
  ## Each row block: the effect matrix, and its columns in reference display
  ## order so the blocks line up cluster-for-cluster with SMESH-PALM.
  meta_panels <- list(
    list(label = "SMESH-PALM",     AA = beta,                      cl = ref_cluster),
    list(label = "SMESH-ANCOMBC2", AA = ANCOMBC2_model$disease$mu, cl = ancombc2_cluster),
    list(label = "SMESH-MaAsLin3", AA = MaAsLin2_model$disease$mu, cl = cluster_vec(MaAsLin2_model)),
    list(label = "SMESH-LinDA",    AA = LinDA_model$disease$mu,    cl = cluster_vec(LinDA_model)),
    list(label = "SKM+FE",         AA = SKM_detect,                cl = SKM_tab),
    list(label = "SHC+FE",         AA = SHC_detect,                cl = SHC_tab)
  )
  
  g_A_lst <- lapply(meta_panels, function(p) {
    cols <- tag_cols_in_ref_order(ref_cluster, p$cl, cluster_order)
    plot.single.study.heatmap.meta(
      x            = rev(cols),          # y is drawn bottom-up
      species_lst  = species_lst,
      G            = G,
      cluster_cols = cluster_palette(G, prefix = "C"),
      AA           = p$AA,
      mx           = 2,
      step         = 0.01,
      text         = FALSE
    )
  })
  
  ## Melody fits one pooled model, so it contributes a single row.
  melody_mat <- matrix(Melody_mod$disease$coef, ncol = 1,
                       dimnames = list(names(Melody_mod$disease$coef), "1"))
  g_A_lst <- c(g_A_lst, list(
    plot.single.study.heatmap.meta(
      x            = 1,
      species_lst  = species_lst,
      G            = G,
      cluster_cols = cluster_palette(G, prefix = "C"),
      AA           = melody_mat,
      mx           = 2,
      step         = 0.01,
      text         = FALSE
    )
  ))
  
  g_A <- ggpubr::ggarrange(plotlist = g_A_lst, ncol = 1)
  
  ggsave(
    filename = file.path(fig_dir, "SFig6_A.png"),
    plot = g_A,
    width = 700, height = 320, units = "mm", dpi = 300
  )
  
  
  # =============================================================================
  #  Panel B - number of clusters selected
  # =============================================================================
  
  runs <- c(seq_len(n_studies), 100)
  
  G_tab <- bind_rows(lapply(names(gic_methods), function(label) {
    tibble(
      Method = label,
      G      = vapply(runs, function(d) select_G(gic_methods[[label]], d), numeric(1)),
      type   = ifelse(runs == 100, "ALL", "LOCO")
    )
  }))
  
  ## SKM and SHC take G as given rather than selecting it.
  G_tab <- bind_rows(G_tab, expand_grid(
    Method = c("SKM+FE", "SHC+FE"),
    tibble(G = 1, type = ifelse(runs == 100, "ALL", "LOCO"))
  ))
  
  method_order_B <- c(names(gic_methods), "SKM+FE", "SHC+FE")
  
  g_B <- plot_G_box(df = G_tab, method_order = method_order_B)
  
  ggsave(
    filename = file.path(fig_dir, "SFig6_B.png"),
    plot = g_B,
    width = 110, height = 110, units = "mm", dpi = 300
  )
  
  
  # =============================================================================
  #  Panel C - leave-one-context-out retention, stability and ARI
  # =============================================================================
  
  study_id <- names(summary_stat_meta_filter)
  
  ## G chosen on the full data, per method; the leave-one-context-out refits are read at
  ## that same G so the comparison is like for like.
  best_G <- vapply(gic_methods, function(m) select_G(m, 100), numeric(1))
  
  res <- bind_rows(lapply(names(gic_methods), function(label) {
    m <- gic_methods[[label]]
    cluster_global <- cluster_vec(load_fit(m, best_G[[label]], 100))
  
    ## One leave-one-context-out refit per held-out study.
    runs_lst <- lapply(seq_along(study_id), function(s) {
      cl_sub <- cluster_vec(load_fit(m, best_G[[label]], s))
      cl_sub <- cl_sub[intersect(names(cl_sub), study_id)]
      list(
        retain = retain_membership(cl_sub, cluster_global[names(cl_sub)]),
        ari    = mclust::adjustedRandIndex(cl_sub, cluster_global[names(cl_sub)])
      )
    })
    names(runs_lst) <- study_id
  
    ## Stability of study d: how often it keeps its membership across the refits
    ## it took part in.
    stability <- vapply(study_id, function(d) {
      mean(vapply(runs_lst, function(r) {
        if (d %in% names(r$retain$keep_status)) {
          as.numeric(r$retain$keep_status[[d]])
        } else {
          NA_real_
        }
      }, numeric(1)), na.rm = TRUE)
    }, numeric(1))
  
    tibble(
      Method    = label,
      retain    = vapply(runs_lst, function(r) r$retain$retain_rate, numeric(1)),
      stability = stability,
      ari       = vapply(runs_lst, function(r) r$ari, numeric(1)),
      study     = study_id
    )
  }))
  
  res$Method <- factor(res$Method, levels = names(gic_methods))
  
  ## The three metrics share everything but the y variable and its label.
  ## `coord_cartesian` zooms rather than filtering, so nothing is dropped
  ## silently if a value falls outside the nominal range.
  metric_box <- function(y, y_lab, y_lim = c(0, 1)) {
    ggplot(res, aes(x = Method, y = .data[[y]])) +
      geom_boxplot(outlier.shape = 16) +
      coord_cartesian(ylim = y_lim) +
      theme_bw() +
      labs(x = NULL, y = y_lab) +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.margin = margin(5.5, 10, 5.5, 10, unit = "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
  }
  
  g_C <- patchwork::wrap_plots(
    list(
      metric_box("retain",    "Retention"),
      metric_box("stability", "Stability"),
      ## ARI is bounded by [-1, 1], not [0, 1].
      metric_box("ari",       "ARI", y_lim = range(c(0, 1, res$ari)))
    ),
    nrow = 1
  ) & theme(plot.margin = margin(5.5, 20, 5.5, 5.5))
  
  ggsave(
    filename = file.path(fig_dir, "SFig6_C.png"),
    plot = g_C,
    width = 250, height = 100, units = "mm", dpi = 300
  )
  
  
  # # =============================================================================
  # #  Panel D - cluster assignments across methods
  # # =============================================================================
  # #  Panels D and E use each method at its own selected G, so a method that
  # #  prefers more clusters than the reference gets its own extra colours.
  # 
  # SMESH_best    <- load_fit("SMESH",    best_G[["SMESH-PALM"]],     100)
  # ANCOMBC2_best <- load_fit("ANCOMBC2", best_G[["SMESH-ANCOMBC2"]], 100)
  # LinDA_best    <- load_fit("LinDA",    best_G[["SMESH-LinDA"]],    100)
  # MaAsLin2_best <- load_fit("MaAsLin3", best_G[["SMESH-MaAsLin3"]], 100)
  # 
  # ## Re-derive the display order for the selected-G reference fit.
  # W_best     <- round(SMESH_best$disease$W)
  # G_best     <- ncol(W_best)
  # beta_best  <- mask_unobserved_clusters(SMESH_best$disease$mu, W_best, stats$est)
  # 
  # cluster_order_best <- order_clusters_by_selection(beta_best)
  # study_order_best   <- order_studies_within_clusters(stats$est[sel_rows(beta_best), ],
  #                                                     W_best, cluster_order_best)
  # 
  # ## Head-room for competitors that split further than the reference.
  # n_extra <- max(0, max(vapply(
  #   list(ANCOMBC2_best, LinDA_best, MaAsLin2_best),
  #   function(m) ncol(m$disease$W), numeric(1)
  # )) - G_best)
  # 
  # g_D <- plot.cluster(
  #   x            = rev(cluster_order_best),
  #   cluster_cols = cluster_palette_reversed(G_best, n_extra = n_extra),
  #   study_order  = rev(study_order_best),
  #   ref   = list("SMESH-PALM" = cluster_vec(SMESH_best)),
  #   other = list(
  #     "SMESH-ANCOMBC2" = cluster_vec(ANCOMBC2_best),
  #     "SMESH-LinDA"    = cluster_vec(LinDA_best),
  #     "SMESH-MaAsLin3" = cluster_vec(MaAsLin2_best)
  #   ),
  #   height_gap = 0.50
  # )
  # 
  # ggsave(
  #   filename = file.path(fig_dir, "SFig6_D.png"),
  #   plot = g_D,
  #   width = 170, height = 60, units = "mm", dpi = 300
  # )
  # 
  # 
  # # =============================================================================
  # #  Panel E - shared vs context-dependent signatures across methods
  # # =============================================================================
  # 
  # selected_in_all <- function(mu) names(which(apply(mu, 1, function(d) all(d != 0))))
  # selected_in_any <- function(mu) names(which(rowSums(mu != 0) > 0))
  # 
  # method_mu <- list(
  #   "SMESH-PALM"     = SMESH_best$disease$mu,
  #   "SMESH-ANCOMBC2" = ANCOMBC2_best$disease$mu,
  #   "SMESH-LinDA"    = LinDA_best$disease$mu,
  #   "SMESH-MaAsLin3" = MaAsLin2_best$disease$mu
  # )
  # 
  # shared_lst   <- lapply(method_mu, selected_in_all)
  # specific_lst <- Map(function(mu, shared) setdiff(selected_in_any(mu), shared),
  #                     method_mu, shared_lst)
  # 
  # g_E <- plot_panel_E_bar(
  #   shared_lst   = shared_lst,
  #   specific_lst = specific_lst,
  #   ref_method   = "SMESH-PALM",
  #   method_order = names(method_mu),
  #   col_shared   = "#222222",
  #   col_specific = "#8FA3B0"
  # )
  # 
  # ggsave(
  #   filename = file.path(fig_dir, "SFig6_E.png"),
  #   plot = g_E,
  #   width = 350, height = 100, units = "mm", dpi = 300
  # )
