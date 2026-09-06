# =============================================================================
#   Figure S3 - pan-disease application, sensitivity analyses
#
#  Run from the project root: Rscript Rscript/FigureS3.R
# =============================================================================
#  Panel A  cluster-level effects, SMESH-PALM stacked against every competitor
#  Panel B  number of clusters selected, all data vs leave-one-context-out
#  Panel C  retention, stability and ARI under leave-one-context-out
#  Panel D  clustering agreement under feature-filter and taxonomy perturbation
#
#  Display order matches the main figure (see Rscript/Figure4.R): clusters run
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
  loso_dir  <- "./GMrepo_analysis/GMrepo_loso"
  model_dir <- "./GMrepo_analysis/Model"
  data_dir  <- "./GMrepo_analysis/Data"
  fig_dir   <- "./GMrepo_analysis/Figure"
  
  G_grid  <- 2:5      # candidate cluster counts
  cut_off <- 0.015     # relative GIC gain below which we stop adding clusters
  
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
  #' Falls back to the largest candidate when the gain never flattens out.  Panel
  #' B and Panels C-E share this function, so the G they report and the G they
  #' plot cannot drift apart.
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
  G <- 2
  
  SMESH_model    <- load_fit("SMESH",    G, 100)
  ANCOMBC2_model <- load_fit("ANCOMBC2", G, 100)
  LinDA_model    <- load_fit("LinDA",    G, 100)
  MaAsLin2_model <- load_fit("MaAsLin3", G, 100)
  
  SHC_tab <- load_tab("SHC", G, 100)
  SKM_tab <- load_tab("SKM", G, 100)
  
  load(file.path(model_dir, "SHC_FE_G2.Rdata"));    SHC_detect <- detect.signal
  load(file.path(model_dir, "SKMean_FE_G2.Rdata")); SKM_detect <- detect.signal
  load(file.path(model_dir, "Melody_model.Rdata"))  # Melody_mod
  
  W           <- round(SMESH_model$disease$W)
  feature.ID  <- rownames(SMESH_model$disease$mu)
  ref_cluster <- cluster_vec(SMESH_model)
  
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
    list(label = "SMESH-ANCOMBC2", AA = ANCOMBC2_model$disease$mu, cl = cluster_vec(ANCOMBC2_model)),
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
    ## ggarrange splits the height evenly, so the 7th row needs proportionally
    ## more canvas to keep the blocks the same size as before.
    width = 700, height = 292, units = "mm", dpi = 300
  )
  
  
  # =============================================================================
  #  Panel B - number of clusters selected
  # =============================================================================
  
  study_id  <- names(summary_stat_meta_filter)
  n_studies <- length(study_id)
  runs      <- c(seq_len(n_studies), 100)
  
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
        axis.text.y = element_text(size = 12),
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
  
  
  # =============================================================================
  #  Panel D - sensitivity to the feature filter and to taxonomic level
  # =============================================================================
  #  Both panels re-run the whole pipeline on a perturbed version of the input
  #  data and ask how much of the primary result survives:
  #
  #    left   feature filter   species estimable in > 8 diseases (primary, 244)
  #                            vs > 4 diseases (relaxed, 339)
  #    right  taxonomic level  species (primary) vs genus, aggregating species
  #                            counts by the NCBI lineage of their taxon id
  #
  #  The panel draws the adjusted Rand index between the perturbed and primary
  #  context clusterings.  The Jaccard overlap of the selected signatures is
  #  computed alongside and written to CSV; it is restricted to features
  #  estimable under BOTH settings, otherwise it would partly re-measure the
  #  perturbation itself, since a feature one setting cannot estimate can never
  #  be recovered there.

  sens_dir <- "./GMrepo_analysis/Sensitivity"

  #' Features a set of summaries could actually speak to.  The matrices carry the
  #' full universe with filtered-out features as NA, so the estimable set is read
  #' off the summaries rather than off rownames().
  tested_of <- function(sums) {
    rn <- rownames(sums[[1]]$est)
    rn[Reduce(`|`, lapply(sums, function(s) !is.na(s$est[, 1])))]
  }
  selected_of <- function(mu) rownames(mu)[rowSums(mu != 0, na.rm = TRUE) > 0]

  load_all <- function(file) {
    e <- new.env(); load(file, envir = e); e
  }

  ## The primary species-level fits at G = 2, already loaded for Panel A.
  primary_mu <- list(
    "SMESH-PALM"     = SMESH_model$disease$mu,
    "SMESH-ANCOMBC2" = ANCOMBC2_model$disease$mu,
    "SMESH-LinDA"    = LinDA_model$disease$mu,
    "SMESH-MaAsLin3" = MaAsLin2_model$disease$mu
  )
  primary_cl <- list(
    "SMESH-PALM"     = ref_cluster,
    "SMESH-ANCOMBC2" = cluster_vec(ANCOMBC2_model),
    "SMESH-LinDA"    = cluster_vec(LinDA_model),
    "SMESH-MaAsLin3" = cluster_vec(MaAsLin2_model)
  )
  sens_methods <- names(primary_mu)

  #' One sensitivity comparison.  `map_primary` lifts the primary selections into
  #' the perturbed feature space (identity for the filter, species -> genus for
  #' the taxonomy comparison).
  sens_table <- function(file, universe, map_primary = identity) {
    e <- load_all(file)
    alt_mu <- list(
      "SMESH-PALM"     = e$SMESH_model$disease$mu,
      "SMESH-ANCOMBC2" = e$ANCOMBC2_model$disease$mu,
      "SMESH-LinDA"    = e$LinDA_model$disease$mu,
      "SMESH-MaAsLin3" = e$MaAsLin2_model$disease$mu
    )
    alt_cl <- lapply(list(e$SMESH_model, e$ANCOMBC2_model, e$LinDA_model, e$MaAsLin2_model),
                     cluster_vec)
    names(alt_cl) <- sens_methods

    bind_rows(lapply(sens_methods, function(m) {
      sel_a <- intersect(selected_of(alt_mu[[m]]), universe)
      sel_p <- intersect(map_primary(selected_of(primary_mu[[m]])), universe)
      common <- intersect(names(alt_cl[[m]]), names(primary_cl[[m]]))
      tibble(
        Method  = m,
        jaccard = jaccard_index(sel_a, sel_p),
        ari     = mclust::adjustedRandIndex(alt_cl[[m]][common], primary_cl[[m]][common])
      )
    }))
  }

  #' One bar per method, value printed at the bar end.
  sens_bar <- function(df, value, title, x_lab) {
    df <- df %>% mutate(Method = factor(Method, levels = rev(sens_methods)),
                        label  = sprintf("%.2f", .data[[value]]))
    ggplot(df, aes(x = .data[[value]], y = Method)) +
      geom_col(fill = "#6F678F", width = 0.66) +
      geom_text(aes(label = label), hjust = -0.15, size = 3.8) +
      scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
      labs(title = title, x = x_lab, y = NULL) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )
  }

  ## ---- inputs for both halves of the panel ------------------------------------
  filter_file <- file.path(sens_dir, "res_species_X4_G2.Rdata")
  genus_file  <- file.path(sens_dir, "res_genus_G2.Rdata")
  stopifnot(file.exists(filter_file), file.exists(genus_file))

  tested_primary <- tested_of(summary_stat_meta_filter)

  ## Left half: relaxed feature filter, species level throughout.
  tested_relaxed  <- tested_of(load_all(filter_file)$summary_stat_meta_filter)
  universe_filter <- intersect(tested_primary, tested_relaxed)
  message("Panel D left : ", length(universe_filter), " species estimable under both filters")
  tab_filter <- sens_table(filter_file, universe_filter)

  ## Right half: genus level.  Species -> genus by the NCBI lineage of the taxon
  ## id carried in the feature name, not by parsing the name: "[Clostridium]
  ## symbiosum" is genus Otoolea.
  genus_lookup   <- read.csv(file.path("./Data/GMrepo", "taxon_id_genus_lookup.csv"),
                             stringsAsFactors = FALSE, colClasses = "character")
  taxid_to_genus <- setNames(genus_lookup$genus_name, genus_lookup$species_taxon_id)
  to_genus <- function(species) {
    g <- taxid_to_genus[sub(".*\\[(\\d+)\\]$", "\\1", species)]
    g[is.na(g) | g == ""] <- "Unclassified"
    setdiff(unique(unname(g)), "Unclassified")
  }

  tested_genus   <- tested_of(load_all(genus_file)$summary_stat_meta_filter)
  universe_genus <- intersect(setdiff(tested_genus, "Unclassified"),
                              to_genus(tested_primary))
  message("Panel D right: ", length(universe_genus), " genera estimable at both resolutions")
  tab_genus <- sens_table(genus_file, universe_genus, map_primary = to_genus)

  print(bind_rows(mutate(tab_filter, comparison = "feature filter"),
                  mutate(tab_genus,  comparison = "taxonomy")))

  ## Only the clustering-agreement metric is drawn; the signature overlap is kept
  ## in the CSV for the text.
  g_D <- patchwork::wrap_plots(
    list(
      ## Two-line titles; `atop` stacks them and the theme's bold face applies to
      ## both, so K comes out bold-italic as in the requested layout.
      sens_bar(tab_filter, "ari",
               bquote(atop("Feature-filtering sensitivity",
                           italic(K) ~ "= 244 versus" ~ italic(K) ~ "= 339 retained species")),
               "ARI between clusterings"),
      sens_bar(tab_genus, "ari",
               bquote(atop("Taxonomic-resolution sensitivity",
                           "Species-level versus genus-level analysis")),
               "ARI between clusterings")
    ),
    nrow = 1
  ) & theme(plot.margin = margin(5.5, 12, 5.5, 5.5))

  ggsave(
    filename = file.path(fig_dir, "SFig6_D.png"),
    plot = g_D,
    width = 300, height = 75, units = "mm", dpi = 300
  )

  write.csv(bind_rows(mutate(tab_filter, comparison = "feature filter (X>8 vs X>4)"),
                      mutate(tab_genus,  comparison = "taxonomy (species vs genus)")),
            file.path(fig_dir, "SFig6_D_sensitivity.csv"), row.names = FALSE)

