# =============================================================================
#   Figure 6 - colorectal neoplasia application (26 stage-stratified studies)
# =============================================================================
#  Panel A  per-study effect heat map (A1) over the cluster-level effects (A2)
#  Panel B  leave-one-context-out co-clustering consensus
#  Panel C  forest plots for representative signatures
#  Panel D  cluster assignments across methods
#  Panel E  signature counts and Jaccard overlap across methods
#
#  Display order is decided once, in Panel A, and reused by every later panel:
#    * clusters run most parsimonious first, which here puts the smaller
#      cluster (9 studies, adenoma + early) above the larger one
#      (17 studies, early + late);
#    * studies inside a cluster are ordered by disease stage,
#      Adenoma -> early -> late.
#  Cluster colours are keyed by display position, so the panels stay consistent
#  with each other no matter how the ordering comes out.
#
#  Run from the project root: Rscript Rscript/Figure6.R
# =============================================================================

  rm(list = ls())

  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)

  source("./utility/heatmap_util.R")

  ## ---- configuration ---------------------------------------------------------
  loso_dir  <- "./CRC_analysis/CRC_loso"
  model_dir <- "./CRC_analysis/Model"
  data_dir  <- "./CRC_analysis/Data"
  fig_dir   <- "./CRC_analysis/Figure"

  fit_tag   <- "Model2"   # the G = 2 fit reported in the paper
  n_studies <- 26         # leave-one-context-out runs; "s100" = all contexts

  ## Studies are named "<study>:<stage>"; this is the order stages appear in.
  stage_levels <- c("Adenoma", "early", "late")

  ## ---- study display labels --------------------------------------------------
  ## The source repository ships three identifiers that carry the wrong year or a
  ## typo, and two that are too long for an axis.  Piccinno et al. (2025) use:
  ##   WirbelJ_2018 -> WirbelJ_2019 (Nat Med 2019),  YuJ_2015 -> YuJ_2017 (Gut 2017),
  ##   NSHII -> NHSII (Nurses' Health Study II),
  ##   ONCOBIOME_IIGM_CZ / _IT -> Cohorts 3 and 4.
  ##
  ## The names are lookup keys as well as labels (Panel B indexes the consensus
  ## matrix by them, Panels C and D join on them), so the relabelling is applied
  ## to every loaded object rather than to the axis text - otherwise the joins
  ## would silently drop rows.  Only Panel A1 prints them.
  ##
  ## MUST match `study_fig` in Make_TableS3.R.
  rename_study <- c(
    "WirbelJ_2018"      = "WirbelJ_2019",
    "YuJ_2015"          = "YuJ_2017",
    "NSHII"             = "NHSII",
    "ONCOBIOME_IIGM_CZ" = "IIGM_CZ",
    "ONCOBIOME_IIGM_IT" = "IIGM_IT"
  )

  #' Relabel "<study>:<stage>" keys, leaving the stage untouched.
  relabel_study <- function(x) {
    s <- sub(":.*$", "", x)
    g <- sub("^.*:", "", x)
    paste0(ifelse(s %in% names(rename_study), rename_study[s], s), ":", g)
  }

  #' Relabel the study names of a fitted model's membership matrix.
  relabel_model <- function(m) {
    rownames(m$disease$W) <- relabel_study(rownames(m$disease$W))
    m
  }

  ## `Model2_SMESH_s100.Rdata` and friends.
  loso_file <- function(method, s, tag = fit_tag) {
    file.path(loso_dir, paste0(tag, "_", method, "_s", s, ".Rdata"))
  }

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
#  Panel A - effect heat maps
# =============================================================================

  load(file.path(data_dir, "Summary_stat_SMESH.Rdata"))   # summary_stat_filter
  load(loso_file("SMESH", 100))                           # SMESH_model

  names(summary_stat_filter) <- relabel_study(names(summary_stat_filter))
  SMESH_model <- relabel_model(SMESH_model)

  G          <- ncol(SMESH_model$disease$W)
  W          <- round(SMESH_model$disease$W)              # study x cluster
  feature.ID <- rownames(SMESH_model$disease$mu)

  ## Per-study estimates, standard errors and FDR q-values, aligned on feature.ID.
  stats <- summary_stat_matrices(summary_stat_filter, feature.ID)

  ## A cluster whose member studies never observed a taxon cannot speak to it, so
  ## drop those effects rather than let them read as a selected zero.
  beta <- mask_unobserved_clusters(SMESH_model$disease$mu, W, stats$est)

  ## ---- display order ---------------------------------------------------------
  ## Clusters: fewest selected features first.  Here that is the 9-study
  ## adenoma/early cluster (30 signatures) above the 17-study early/late one
  ## (89 signatures).
  cluster_order <- order_clusters_by_selection(beta)

  ## Studies: by disease stage within each cluster.
  stage <- setNames(sub(".*:", "", rownames(W)), rownames(W))
  study_order <- order_studies_by_key(
    key           = stage,
    W             = W,
    cluster_order = cluster_order,
    group_levels  = stage_levels
  )

  ## Features grouped by exactly which clusters select them, named by display
  ## position ("1|2" = all-cluster shared, "2" = specific to the second).
  species_lst <- build_signature_groups(beta, cluster_order)

  ## Membership matrix with the clusters in display order.
  W_mat <- W[, cluster_order, drop = FALSE]
  colnames(W_mat) <- paste0("C", seq_len(G))

  ## ---- A1: one row per study -------------------------------------------------
  g_A1 <- plot.single.study.heatmap.ref(
    x            = rev(study_order),   # y is drawn bottom-up
    species_lst  = species_lst,
    cluster_cols = cluster_palette(G),
    AA           = stats$est,
    AA.test      = stats$pval,
    AA.test.q    = stats$qval,
    W            = W_mat,
    mx           = 2,
    step         = 0.001
  )

  ggsave(
    filename = file.path(fig_dir, "FigureA1.png"),
    plot = g_A1,
    width = 700, height = 220, units = "mm", dpi = 300
  )

  ## ---- A2: one row per cluster (meta effects) --------------------------------
  g_A2 <- plot.single.study.heatmap.meta(
    x            = rev(cluster_order),
    species_lst  = species_lst,
    G            = G,
    cluster_cols = cluster_palette(G, prefix = "C"),
    AA           = beta,
    mx           = 2,
    step         = 0.01
  )

  ggsave(
    filename = file.path(fig_dir, "FigureA2.png"),
    plot = g_A2,
    width = 700, height = 150, units = "mm", dpi = 300
  )


# =============================================================================
#  Panel B - leave-one-context-out consensus
# =============================================================================

  ## Cluster assignment from each leave-one-context-out refit.
  cluster_list <- lapply(seq_len(n_studies), function(s) {
    f <- loso_file("SMESH", s)
    if (!file.exists(f)) return(NULL)
    load(f)   # SMESH_model, local to this call
    setNames(apply(SMESH_model$disease$W, 1, which.max),
             relabel_study(rownames(SMESH_model$disease$W)))
  })
  names(cluster_list) <- paste0("s", seq_len(n_studies))
  cluster_list <- Filter(Negate(is.null), cluster_list)

  ## Co-clustering frequency, restricted to the runs where both studies were in.
  C <- build_consensus_lodo(cluster_list)$C

  df_consensus <- as.data.frame(C[study_order, study_order, drop = FALSE]) %>%
    rownames_to_column("i") %>%
    pivot_longer(-i, names_to = "j", values_to = "p") %>%
    mutate(
      i = factor(i, levels = rev(study_order)),
      j = factor(j, levels = study_order)
    )

  ## Outline the clusters on the diagonal; sizes follow the display order.
  rect_df <- consensus_block_rects(colSums(W_mat))

  g_B <- ggplot(df_consensus, aes(x = j, y = i, fill = p)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
      fill = NA,
      linewidth = 3,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = setNames(CLUSTER_COLS[seq_len(G)], rect_df$group)) +
    scale_fill_gradientn(
      colours = c("white", "#FDFFE5", "#7AD151", "#22A884", "#2A788E", "#1E5A6A", "#174552"),
      limits = c(0, 1),
      oob = scales::squish,
      name = "Co-clustering\nfrequency"
    ) +
    coord_fixed() +
    labs(title = "Leave-one-context-out consensus", x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(fig_dir, "FigureB.png"),
    plot = g_B,
    width = 135, height = 100, units = "mm", dpi = 300
  )


# =============================================================================
#  Panel C - forest plots for representative signatures
# =============================================================================

  ## Taxa to highlight.  Their signature type ("All-cluster shared (+)",
  ## "C1-specific", ...) is derived from `beta` rather than hard-coded, so the
  ## subtitles stay correct if the fit or the cluster order changes.
  ##
  ## To pick a taxon of a given type, read it off `species_lst`, whose names are
  ## the display positions sharing the signature:
  ##   species_lst[["1|2"]]  all-cluster shared
  ##   species_lst[["1"]]    specific to C1 (the adenoma/early cluster)
  panel_c_taxa <- c(
    "s_Bifidobacterium_longum",
    "s_Fusobacterium_nucleatum",
    "s_Phascolarctobacterium_succinatutens",
    "s_Clostridium_symbiosum"
  )

  ## signature_label() lives in utility/heatmap_util.R; direction = "none"
  ## suppresses the sign suffix, since the intervals already show direction.

  ## Strip the "s_" prefix and underscores for the italic subtitle.
  pretty_taxon <- function(x) gsub("_", " ", sub("^s_", "", x))

  ## Display position of each study's cluster.
  study_cluster <- setNames(match(apply(W, 1, which.max), cluster_order),
                            rownames(W))

  g_C_lst <- lapply(panel_c_taxa, function(taxon) {
    df_forest <- tibble(
      id      = factor(colnames(stats$est), levels = rev(study_order)),
      est     = as.numeric(stats$est[taxon, ]),
      se      = as.numeric(stats$stderr[taxon, ]),
      cluster = factor(study_cluster[colnames(stats$est)], levels = seq_len(G))
    ) %>%
      mutate(
        lower = est - 1.96 * se,
        upper = est + 1.96 * se
      )

    ggplot(df_forest, aes(y = id, x = est, color = cluster)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
      geom_errorbar(aes(xmin = lower, xmax = upper),
                    orientation = "y", width = 0, linewidth = 1) +
      geom_point(size = 2) +
      ggtitle(bquote(atop(bold(.(signature_label(taxon, beta, cluster_order, direction = "none"))),
                          italic(.(pretty_taxon(taxon)))))) +
      facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
      scale_color_manual(values = cluster_palette(G)) +
      labs(x = "Effect estimate (95% CI)", y = NULL, color = NULL) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
  })

  g_C <- wrap_plots(g_C_lst, nrow = 1) &
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

  ggsave(
    filename = file.path(fig_dir, "FigureC.png"),
    plot = g_C,
    width = 360, height = 120, units = "mm", dpi = 300
  )


# =============================================================================
#  Panel D - cluster assignments across methods
# =============================================================================

  load(loso_file("SHC", 100));  shier_tab <- tab   # cluster / study data frame
  load(loso_file("SKM", 100));  skm_tab   <- tab
  load(loso_file("ANCOMBC2", 100))                 # ANCOMBC2_model
  load(loso_file("MaAsLin3", 100))                 # MaAsLin2_model (legacy name)
  load(loso_file("LinDA", 100))                    # LinDA_model

  shier_tab$study <- relabel_study(shier_tab$study)
  skm_tab$study   <- relabel_study(skm_tab$study)
  ANCOMBC2_model  <- relabel_model(ANCOMBC2_model)
  MaAsLin2_model  <- relabel_model(MaAsLin2_model)
  LinDA_model     <- relabel_model(LinDA_model)

  g_D <- plot.cluster(
    x            = rev(cluster_order),   # plot.cluster relabels G..1 along this
    ## Every method here lands on G clusters; raise n_extra if one ever splits
    ## further, otherwise plot.cluster() stops with "Not enough colors".
    cluster_cols = cluster_palette_reversed(G, n_extra = 0),
    study_order  = rev(study_order),
    ref   = list("SMESH-PALM" = apply(SMESH_model$disease$W, 1, which.max)),
    other = list(
      "SKM"            = skm_tab   %>% pull(cluster, name = study),
      "SHC"            = shier_tab %>% pull(cluster, name = study),
      "SMESH-ANCOMBC2" = apply(ANCOMBC2_model$disease$W, 1, which.max),
      "SMESH-LinDA"    = apply(LinDA_model$disease$W, 1, which.max),
      "SMESH-MaAsLin3" = apply(MaAsLin2_model$disease$W, 1, which.max)
    ),
    height_gap = 0.5
  )

  ggsave(
    filename = file.path(fig_dir, "FigureD.png"),
    plot = g_D,
    width = 160, height = 60, units = "mm", dpi = 300
  )


# =============================================================================
#  Panel E - shared vs context-dependent signatures across methods
# =============================================================================
#  Reads the same fits as Panel D.  The previous version of this script pulled
#  ANCOMBC2 / LinDA / MaAsLin3 from CRC_loso_v2, which is a different feature
#  space (379 vs 274 taxa) from the assignments drawn in Panel D; point
#  `loso_dir` at CRC_loso_v2 above to go back to that.

  load(file.path(model_dir, "SHC_FE_G2.Rdata"));    SHC_detect <- detect.signal
  load(file.path(model_dir, "SKMean_FE_G2.Rdata")); SKM_detect <- detect.signal
  load(file.path(model_dir, "Melody_model.Rdata"))  # Melody_mod

  ## A signature is "shared" when every cluster selects it, and
  ## "context-dependent" when at least one - but not all - cluster does.
  selected_in_all <- function(mu) names(which(apply(mu, 1, function(d) all(d != 0))))
  selected_in_any <- function(mu) names(which(rowSums(mu != 0) > 0))

  method_mu <- list(
    "SMESH-PALM"     = SMESH_model$disease$mu,
    "SMESH-ANCOMBC2" = ANCOMBC2_model$disease$mu,
    "SMESH-MaAsLin3" = MaAsLin2_model$disease$mu,
    "SMESH-LinDA"    = LinDA_model$disease$mu,
    "SKM+FE"         = SKM_detect,
    "SHC+FE"         = SHC_detect
  )

  shared_lst   <- lapply(method_mu, selected_in_all)
  specific_lst <- Map(function(mu, shared) setdiff(selected_in_any(mu), shared),
                      method_mu, shared_lst)

  ## Melody fits a single pooled model, so every signature it finds is shared and
  ## it has no context-dependent set to report.
  shared_lst[["Melody"]]   <- names(which(Melody_mod$disease$coef != 0))
  specific_lst[["Melody"]] <- character(0)

  g_E <- plot_panel_E_bar(
    shared_lst   = shared_lst,
    specific_lst = specific_lst,
    ref_method   = "SMESH-PALM",
    method_order = c(
      "SMESH-PALM",
      "SMESH-ANCOMBC2",
      "SMESH-LinDA",
      "SMESH-MaAsLin3",
      "SKM+FE",
      "SHC+FE",
      "Melody"
    ),
    col_shared   = "#222222",
    col_specific = "#8FA3B0"
  )

  ggsave(
    filename = file.path(fig_dir, "FigureE.png"),
    plot = g_E,
    width = 350, height = 100, units = "mm", dpi = 300
  )
