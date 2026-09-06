# =============================================================================
#  Shared helpers for the MOSIAC figure scripts
# =============================================================================
#
#  Sections
#    1. Shared constants (cluster palette, method ordering)
#    2. Summary-statistic assembly
#    3. Display ordering (clusters, studies, signature groups)
#    4. Heat maps
#    5. Method-comparison panels
#    6. Cluster-stability helpers
#
#  Vocabulary used throughout:
#    * cluster label   - the raw column index of `W` / `mu` coming out of a model.
#    * display position - 1..G, the order clusters are drawn in (1 = first = top
#                        of a heat map).  `cluster_order[k]` is the raw label
#                        shown at display position `k`.
#  Panels stay colour-consistent because every colour is keyed by *display
#  position*, never by the raw cluster label.
# =============================================================================

library(patchwork)


# =============================================================================
#  1. Shared constants
# =============================================================================

## Cluster colours, indexed by display position.
CLUSTER_COLS <- c(
  "#C8C400",  # 1 olive-yellow
  "#E64B9A",  # 2 pink
  "#00BFC4",  # 3 cyan
  "#F28E00",  # 4 orange
  "#39B600",  # 5 green
  "#F91A00",  # 6 red
  "#A6417A"   # 7 purple
)

## Row order for the method-comparison panels.  Methods not listed here keep the
## order in which the caller supplied them and are appended at the end.
METHOD_DISPLAY_ORDER <- c(
  "SMESH-PALM", "SMESH-ANCOMBC2", "SMESH-LinDA", "SMESH-MaAsLin3",
  "SKM", "SHC", "SKM+FE", "SHC+FE", "Melody"
)

#' Cluster colours keyed by display position.
#'
#' @param n      number of clusters to colour.
#' @param prefix name prefix, e.g. "C" for the "C1".."CG" keys used by the
#'               meta-effect heat map.
cluster_palette <- function(n, prefix = "") {
  if (n > length(CLUSTER_COLS)) {
    stop("cluster_palette(): only ", length(CLUSTER_COLS),
         " cluster colours are defined, ", n, " requested.")
  }
  setNames(CLUSTER_COLS[seq_len(n)], paste0(prefix, seq_len(n)))
}

#' Cluster colours for `plot.cluster()`.
#'
#' `plot.cluster()` relabels clusters as G..1 along the reversed display order,
#' so its palette has to be keyed the other way round.  Extra slots are for
#' competing methods that split the data into more clusters than the reference.
#'
#' @param G       number of reference clusters.
#' @param n_extra number of additional colours for unmatched method clusters.
cluster_palette_reversed <- function(G, n_extra = 0) {
  if (G + n_extra > length(CLUSTER_COLS)) {
    stop("cluster_palette_reversed(): only ", length(CLUSTER_COLS),
         " cluster colours are defined, ", G + n_extra, " requested.")
  }
  setNames(
    c(rev(CLUSTER_COLS[seq_len(G)]), CLUSTER_COLS[G + seq_len(n_extra)]),
    as.character(seq_len(G + n_extra))
  )
}

## Preferred order first, anything unknown appended in the order given.
order_methods <- function(methods, preferred = METHOD_DISPLAY_ORDER) {
  c(intersect(preferred, methods), setdiff(methods, preferred))
}

## plotmath labels for a method axis: the reference method is bold, the rest
## plain.  Returned as strings to be fed through parse().
method_label_expr <- function(methods, ref_method = NULL) {
  setNames(
    ifelse(methods %in% ref_method,
           paste0("bold('", methods, "')"),
           paste0("'", methods, "'")),
    methods
  )
}


# =============================================================================
#  2. Summary-statistic assembly
# =============================================================================

#' Stack per-study summary statistics into feature x study matrices.
#'
#' Features a study never observed stay `NA` in `est` (so that
#' `mask_unobserved_clusters()` can find them) and are given a p/q-value of 1.1,
#' which is the sentinel the heat maps use to suppress a point.
#'
#' @param summary_stats named list; each element has `$est` and `$stderr`
#'                      matrices with features in the rows.
#' @param feature_id    features to keep, in the desired row order.
#'
#' @return list of feature x study matrices: `est`, `stderr`, `pval`, `qval`.
summary_stat_matrices <- function(summary_stats, feature_id) {
  studies <- names(summary_stats)
  est  <- matrix(NA_real_, length(feature_id), length(studies),
                 dimnames = list(feature_id, studies))
  stderr <- pval <- qval <- est

  for (l in studies) {
    s  <- summary_stats[[l]]
    id <- intersect(feature_id, rownames(s$est))

    est[id, l]    <- s$est[id, 1]
    stderr[id, l] <- s$stderr[id, 1]
    pval[id, l]   <- 1 - pchisq((s$est[id, 1] / s$stderr[id, 1])^2, df = 1)
    ## FDR over the full feature set, so studies stay comparable even when a
    ## few features are missing.
    qval[id, l]   <- p.adjust(pval[id, l], method = "fdr", n = length(feature_id))
  }

  pval[is.na(pval)] <- 1.1
  qval[is.na(qval)] <- 1.1

  list(est = est, stderr = stderr, pval = pval, qval = qval)
}

#' Zero out cluster-level effects that no study in the cluster could observe.
#'
#' @param beta     feature x cluster effect matrix.
#' @param W        study x cluster membership matrix (rows named by study).
#' @param est_mat  feature x study estimates, rows aligned with `beta`.
mask_unobserved_clusters <- function(beta, W, est_mat) {
  observed <- apply(W, 2, function(w) {
    members <- rownames(W)[w != 0]
    apply(est_mat[, members, drop = FALSE], 1, function(x) any(!is.na(x)))
  })
  beta[!observed] <- 0
  beta
}


# =============================================================================
#  3. Display ordering
# =============================================================================

#' Order clusters by how many features they select.
#'
#' @param beta       feature x cluster effect matrix.
#' @param decreasing `FALSE` (default) puts the most parsimonious cluster first.
#'
#' @return raw cluster labels, in display order.
order_clusters_by_selection <- function(beta, decreasing = FALSE) {
  n_selected <- colSums(beta != 0)
  order(n_selected, decreasing = decreasing)
}

#' Order studies inside each cluster by summary-statistic similarity.
#'
#' Studies of a cluster are ordered by average-linkage hierarchical clustering
#' of their effect-size profiles.  Effects a study is missing are imputed with
#' that study's own mean, which keeps them from driving the distance.
#'
#' @param effect_mat    feature x study effect estimates.
#' @param W             study x cluster membership matrix (rows named by study).
#' @param cluster_order raw cluster labels in display order.
#' @param hclust_method linkage passed to `hclust()`.
#'
#' @return study IDs, in display order (grouped by cluster).
order_studies_within_clusters <- function(effect_mat, W, cluster_order,
                                          hclust_method = "average") {
  unlist(lapply(cluster_order, function(g) {
    ids <- rownames(W)[W[, g] != 0]
    if (length(ids) < 2) return(ids)

    mat <- effect_mat[, ids, drop = FALSE]
    mat <- apply(mat, 2, function(v) {
      v[is.na(v)] <- mean(v, na.rm = TRUE)
      v
    })

    hc <- hclust(dist(t(mat)), method = hclust_method)
    hc$labels[hc$order]
  }), use.names = FALSE)
}

#' Describe how a signature is shared across clusters.
#'
#' Reads the sharing pattern off `beta` rather than being hard-coded, so the
#' forest-plot subtitles stay correct if the fit or the cluster order changes.
#'
#' @param taxon         row name of `beta`.
#' @param beta          feature x cluster effect matrix.
#' @param cluster_order raw cluster labels in display order.
#' @param direction     when to append a "(+)" / "(-)" sign suffix:
#'                      "shared" (default) only for all-cluster-shared
#'                      signatures, "all" for every label, "none" never.
#'                      The suffix is only added when the non-zero effects
#'                      agree in sign.
#'
#' @return e.g. "All-cluster shared (-)", "C2-specific (+)", "Shared by C1 and C3".
signature_label <- function(taxon, beta, cluster_order,
                            direction = c("shared", "all", "none")) {
  direction <- match.arg(direction)

  eff <- beta[taxon, cluster_order]
  pos <- which(eff != 0)
  if (length(pos) == 0) return("Not selected")

  sign_suffix <- function() {
    if (all(eff[pos] > 0)) " (+)" else if (all(eff[pos] < 0)) " (-)" else ""
  }

  if (length(pos) == ncol(beta)) {
    sfx <- if (direction %in% c("shared", "all")) sign_suffix() else ""
    return(paste0("All-cluster shared", sfx))
  }

  sfx <- if (direction == "all") sign_suffix() else ""

  if (length(pos) == 1) return(paste0("C", pos, "-specific", sfx))

  labs <- paste0("C", pos)
  paste0("Shared by ", paste(paste(head(labs, -1), collapse = ", "),
                             tail(labs, 1), sep = " and "), sfx)
}

#' Order studies inside each cluster by a categorical key.
#'
#' The counterpart to `order_studies_within_clusters()`, for when a meaningful
#' ordering already exists (disease stage, timepoint, ...) and hierarchical
#' clustering of the summary statistics is not what you want.
#'
#' @param key           named vector, study -> group label.
#' @param W             study x cluster membership matrix (rows named by study).
#' @param cluster_order raw cluster labels in display order.
#' @param group_levels  group labels in the order they should appear; defaults
#'                      to `sort(unique(key))`.  Studies tied on the key keep
#'                      the order they appear in `W`.
#'
#' @return study IDs, in display order (grouped by cluster, then by key).
order_studies_by_key <- function(key, W, cluster_order, group_levels = NULL) {
  if (is.null(group_levels)) group_levels <- sort(unique(key))

  unknown <- setdiff(unique(key), group_levels)
  if (length(unknown)) {
    stop("order_studies_by_key(): key values missing from group_levels: ",
         paste(unknown, collapse = ", "))
  }

  unlist(lapply(cluster_order, function(g) {
    ids <- rownames(W)[W[, g] != 0]
    ids[order(match(key[ids], group_levels), seq_along(ids))]
  }), use.names = FALSE)
}

#' Features selected in exactly the clusters `ID`, ordered by total effect.
get_tax <- function(beta, ID) {
  tmp <- apply(beta, 1, function(d) {
    all(d[ID] != 0) && all(d[setdiff(seq_len(ncol(beta)), ID)] == 0)
  })
  names(sort(rowSums(beta)[tmp]))
}

#' Split features into signature groups, one per non-empty subset of clusters.
#'
#' Groups run from the all-cluster-shared signature down to the
#' cluster-specific ones, and are named by the *display positions* they cover
#' ("1|2|3|4", "1|3", "2", ...) so downstream code never has to re-derive them.
#'
#' @param beta          feature x cluster effect matrix.
#' @param cluster_order raw cluster labels in display order.
build_signature_groups <- function(beta, cluster_order) {
  G <- ncol(beta)
  subsets <- unlist(
    lapply(G:1, function(k) combn(seq_len(G), k, simplify = FALSE)),
    recursive = FALSE
  )

  groups <- lapply(subsets, function(pos) get_tax(beta, cluster_order[pos]))
  names(groups) <- vapply(subsets, paste, character(1), collapse = "|")
  groups
}

## Collapse signature groups into the three facet columns the heat maps draw:
##   1 = all-cluster shared, 2 = partially shared, 3 = cluster specific.
signature_block_index <- function(species_lst, G = NULL) {
  n_grp <- length(species_lst)
  nm    <- names(species_lst)

  if (!is.null(nm) && all(nzchar(nm))) {
    ## Named by build_signature_groups(): read the subset size off the name.
    k <- lengths(strsplit(nm, "|", fixed = TRUE))
    G_data <- max(k)
  } else {
    ## Unnamed legacy input: assume the full power set in the same order.
    G_data <- as.integer(round(log2(n_grp + 1)))
    if (n_grp != 2^G_data - 1) {
      stop("signature_block_index(): species_lst must hold every cluster subset ",
           "(2^G - 1 groups), or carry names from build_signature_groups(); ",
           "got ", n_grp, " groups.")
    }
    k <- rep(G_data:1, times = choose(G_data, G_data:1))
  }

  ## `species_lst` is the authority on how many clusters there are; a `G` that
  ## disagrees used to be truncated silently, mislabelling the facet blocks.
  if (!is.null(G) && G != G_data) {
    warning("signature_block_index(): species_lst describes G = ", G_data,
            ", ignoring the supplied G = ", G, ".")
  }

  ifelse(k == G_data, 1L, ifelse(k == 1L, 3L, 2L))
}

## Facet-column label ("C1"/"C2"/"C3") repeated once per feature.
signature_class_id <- function(species_lst, G = NULL) {
  block <- signature_block_index(species_lst, G)
  rep(paste0("C", block), times = lengths(species_lst))
}

#' Outline coordinates for the cluster blocks of a consensus heat map.
#'
#' Assumes the heat map runs left-to-right along the study display order and
#' top-to-bottom along its reverse, i.e. blocks sit on the main diagonal.
#'
#' @param sizes number of studies per cluster, in display order.
consensus_block_rects <- function(sizes) {
  n     <- sum(sizes)
  end   <- cumsum(sizes)
  start <- end - sizes + 1

  data.frame(
    xmin  = start - 0.5,
    xmax  = end + 0.5,
    ymin  = n + 0.5 - end,
    ymax  = n + 1.5 - start,
    group = paste0("block", seq_along(sizes)),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
#  4. Heat maps
# =============================================================================

#' Per-study effect heat map, one row per study, grouped by cluster.
#'
#' @param x          study IDs, bottom-to-top (i.e. reversed display order).
#' @param species_lst signature groups from `build_signature_groups()`.
#' @param cluster_cols colours keyed by display position, `cluster_palette(G)`.
#' @param AA,AA.test,AA.test.q feature x study estimate / p-value / q-value.
#' @param W          study x cluster membership, columns in display order.
#' @param mx,step    colour-scale limit and resolution.
plot.single.study.heatmap.ref <- function(x, species_lst, cluster_cols,
                                          AA, AA.test, AA.test.q, W, mx, step) {

  G <- ncol(W)
  AA[AA >  mx] <-  mx
  AA[AA < -mx] <- -mx

  num.col.steps <- length(seq(-mx, mx, by = step)) - 1
  n <- floor(0.45 * num.col.steps)
  col.hm <- c(
    colorRampPalette(c("#3B6FB6", "#F7F7F7"))(n),
    colorRampPalette(c("#F7F7F7", "#B35836"))(n + 1)[-1]
  )

  taxa     <- unlist(species_lst, use.names = FALSE)
  class_id <- signature_class_id(species_lst, G)

  df.plot <- do.call(rbind, lapply(x, function(l) {
    tibble(
      species     = factor(taxa, levels = taxa),
      cmpd        = factor(l, levels = l),
      AA          = AA[taxa, l],
      pval        = AA.test[taxa, l],
      qval        = AA.test.q[taxa, l],
      cluster     = factor(which.max(W[l, ]), levels = as.character(1:G)),
      HMDB.Class2 = factor(class_id, levels = paste0("C", 1:3))
    )
  }))

  df.plot2 <- tibble(
    cluster = factor(apply(W[x, , drop = FALSE], 1, which.max),
                     levels = as.character(1:G)),
    disease = factor(x, levels = x)
  )

  df.plot <- df.plot %>%
    mutate(
      lp = -log10(pmax(qval, .Machine$double.xmin)),
      ## qval == 1.1 (the "never observed" sentinel) gives lp < 0, which drops
      ## the point entirely rather than drawing a misleading zero effect.
      size_cat = case_when(
        lp < 0            ~ NA_character_,
        lp < -log10(0.05) ~ "q-value >= 0.05",
        TRUE              ~ "q-value < 0.05"
      )
    )

  g1 <- df.plot %>%
    ggplot(aes(x = species, y = cmpd)) +
    geom_point(
      aes(color = AA, size = size_cat),
      shape = 16,
      na.rm = FALSE
    ) +
    facet_grid(
      rows   = vars(cluster),
      cols   = vars(HMDB.Class2),
      scales = "free",
      space  = "free",
      switch = "y"
    ) +
    scale_color_gradientn(
      colours = col.hm,
      values  = scales::rescale(c(
        -mx, -0.6 * mx, -0.3 * mx, -0.1 * mx, 0,
        0.1 * mx, 0.3 * mx, 0.6 * mx, mx
      )),
      limits = c(-mx, mx),
      breaks = c(-mx, 0, mx),
      labels = c(paste0("-", mx), "0", as.character(mx)),
      name   = "Association effect",
      na.value = "white"
    ) +
    scale_size_manual(
      values = c(
        "q-value >= 0.05" = 5,
        "q-value < 0.05"  = 8
      ),
      breaks = c("q-value >= 0.05", "q-value < 0.05"),
      labels = c(
        expression(italic(q) >= 0.05),
        expression(italic(q) < 0.05)
      ),
      name = NULL,
      na.translate = FALSE
    ) +
    guides(
      color = guide_colorbar(
        order = 1,
        direction = "horizontal",
        title.position = "top",
        barwidth = unit(4, "cm"),
        barheight = unit(0.4, "cm"),
        title.theme = element_text(
          size = 20,
          margin = margin(b = 12)
        )
      ),
      size = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = "black", linewidth = 2),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 20)),
      legend.text = element_text(size = 20),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.height = unit(1.2, "cm"),
      legend.spacing.x = unit(1.2, "cm"),
      legend.box.spacing = unit(0.8, "cm"),
      strip.text = element_blank()
    )

  ## Narrow colour strip carrying the study labels and their cluster.
  p_cluster <- ggplot(df.plot2, aes(x = "Cluster", y = disease, fill = cluster)) +
    geom_tile() +
    scale_x_discrete(position = "bottom") +
    scale_fill_manual(values = cluster_cols) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(
      rows   = vars(cluster),
      scales = "free_y",
      space  = "free_y",
      switch = "y"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(size = 21),
      strip.placement = "outside",
      strip.background = element_blank()
    )

  (p_cluster | g1) + plot_layout(widths = c(0.05, 3.5))
}

#' Cluster-level (meta) effect heat map, one row per cluster.
#'
#' @param x           cluster columns of `AA`, bottom-to-top.  Row `k` from the
#'                    bottom is labelled "Ck".
#' @param species_lst signature groups from `build_signature_groups()`.
#' @param G           number of clusters; derived from `species_lst` when NULL.
#' @param cluster_cols colours keyed "C1".."CG", `cluster_palette(G, "C")`.
#' @param AA          feature x cluster effect matrix.
#' @param mx,step     colour-scale limit and resolution.
#' @param text        draw the taxa names along the x axis.
plot.single.study.heatmap.meta <- function(x, species_lst, G = NULL,
                                           cluster_cols, AA, mx, step,
                                           text = TRUE) {

  AA[AA >  mx] <-  mx
  AA[AA < -mx] <- -mx

  taxa <- unlist(species_lst, use.names = FALSE)

  ## Features a model never selected are absent from AA; draw them as zero.
  taxa_miss <- setdiff(taxa, rownames(AA))
  if (length(taxa_miss) > 0) {
    AA <- rbind(AA, matrix(0, nrow = length(taxa_miss), ncol = ncol(AA),
                           dimnames = list(taxa_miss, colnames(AA))))
  }

  num.col.steps <- length(seq(-mx, mx, by = step)) - 1
  n <- floor(0.49 * num.col.steps)
  col.hm <- c(
    colorRampPalette(c("#3B6FB6", "#FFFFFF"))(n),
    colorRampPalette(c("#FFFFFF", "#B35836"))(n + 1)[-1]
  )

  class_id     <- signature_class_id(species_lst, G)
  row_labels   <- paste0("C", length(x):1)

  df.plot <- do.call(rbind, Map(function(l, lab) {
    tibble(
      species     = factor(taxa, levels = taxa),
      cmpd        = factor(lab, levels = lab),
      AA          = AA[taxa, l],
      HMDB.Class2 = factor(class_id, levels = paste0("C", 1:3))
    )
  }, x, row_labels))

  ## Levels stay in drawing order, so "C1" ends up at the top of the strip.
  df.plot2 <- tibble(cluster = factor(row_labels, levels = row_labels))

  g1 <- df.plot %>%
    ggplot(aes(x = species, y = cmpd, fill = AA)) +
    geom_tile() +
    facet_grid(cols = vars(HMDB.Class2),
               scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(
      colours = col.hm,
      values  = scales::rescale(c(
        -mx, -0.6 * mx, -0.3 * mx, -0.1 * mx, 0,
        0.1 * mx, 0.3 * mx, 0.6 * mx, mx
      )),
      limits = c(-mx, mx),
      breaks = c(-mx, 0, mx),
      labels = c(paste0("-", mx), "0", as.character(mx))
    ) +
    scale_x_discrete(position = "bottom") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, colour = "black", linewidth = 2),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.y = element_blank(),
      axis.text.x = if (text) {
        element_text(angle = 45, hjust = 1, vjust = 1, size = 18)
      } else {
        element_blank()
      },
      strip.text = element_blank(),
      legend.position = "none"
    )

  p_cluster <- ggplot(df.plot2, aes(x = "Cluster", y = cluster, fill = cluster)) +
    geom_tile() +
    scale_x_discrete(position = "bottom") +
    scale_fill_manual(values = cluster_cols) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 6), "cm"),
      axis.text.y = element_text(size = 24, face = "bold"),
      strip.placement = "outside",
      strip.background = element_blank()
    )

  (p_cluster | g1) + plot_layout(widths = c(0.05, 3.5))
}


# =============================================================================
#  5. Method-comparison panels
# =============================================================================

#' Compare cluster assignments across methods, one row per method.
#'
#' Competing methods are relabelled onto the reference clusters by greedy
#' best-overlap matching; clusters that find no partner (methods with more
#' clusters than the reference) take the leftover colours.
#'
#' @param x            reference cluster labels in *reversed* display order, so
#'                     the first display cluster gets the highest new label.
#' @param cluster_cols colours from `cluster_palette_reversed()`; must hold one
#'                     slot per reference cluster plus one per extra cluster.
#' @param ref          one-element named list: method name -> named cluster vector.
#' @param other        named list of the same shape for the competing methods.
#' @param study_order  study IDs, right-to-left along the x axis.
#' @param method_order row order top-to-bottom; defaults to `METHOD_DISPLAY_ORDER`.
#' @param height_gap   tile height, i.e. the gap between method rows.
plot.cluster <- function(x,
                         cluster_cols,
                         ref,
                         other = list(),
                         height_gap = 0.8,
                         study_order = NULL,
                         method_order = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required.")
  }

  if (!is.list(ref) || length(ref) != 1) {
    stop("ref must be a named list with exactly one element.")
  }

  ref_name <- names(ref)[1]
  ref_vec  <- ref[[1]]

  if (is.null(ref_name) || ref_name == "") {
    stop("ref must have a method name.")
  }
  if (is.null(names(ref_vec))) {
    stop("Reference cluster vector must be named by study.")
  }

  study_ids <- names(ref_vec)
  ref_vec <- as.integer(ref_vec)
  names(ref_vec) <- study_ids

  ref_labels <- sort(unique(ref_vec))
  G_ref <- length(ref_labels)

  if (length(x) != G_ref || !setequal(x, ref_labels)) {
    stop("x must contain exactly the reference cluster labels, in the desired order.")
  }

  if (is.null(study_order)) {
    study_order <- study_ids
  } else {
    if (!all(study_ids %in% study_order)) {
      stop("study_order must contain all studies in ref.")
    }
    study_order <- unique(study_order)
  }
  study_rank <- setNames(seq_along(study_order), study_order)

  # old ref cluster -> reordered display cluster
  ref_old_to_new <- setNames(seq_len(G_ref), x)
  ref_reordered <- unname(ref_old_to_new[as.character(ref_vec)])
  names(ref_reordered) <- study_ids

  # make sure cluster_cols has enough colors
  if (is.null(names(cluster_cols))) {
    names(cluster_cols) <- as.character(seq_along(cluster_cols))
  }

  map_to_ref_with_extra <- function(method_vec, ref_vec, ref_old_to_new, cluster_cols) {
    if (is.null(names(method_vec))) {
      stop("Each method cluster vector must be named by study.")
    }

    method_ids <- names(method_vec)
    method_vec <- as.integer(method_vec)
    names(method_vec) <- method_ids

    if (!setequal(names(method_vec), names(ref_vec))) {
      stop("Each method must contain exactly the same study IDs as ref.")
    }

    method_vec <- method_vec[names(ref_vec)]
    method_labels <- sort(unique(method_vec))
    overlap_tab <- table(method = method_vec, ref = ref_vec)

    available_new_labels <- as.integer(names(cluster_cols))
    used_new_labels <- integer(0)
    method_to_new <- setNames(rep(NA_integer_, length(method_labels)), method_labels)

    # process method clusters by strongest overlap first
    best_overlap <- sapply(method_labels, function(k) max(overlap_tab[as.character(k), ]))
    method_labels_ord <- method_labels[order(best_overlap, decreasing = TRUE)]

    for (k in method_labels_ord) {
      overlaps <- overlap_tab[as.character(k), ]
      ref_old_ranked <- as.integer(names(sort(overlaps, decreasing = TRUE)))
      ref_new_ranked <- unname(ref_old_to_new[as.character(ref_old_ranked)])

      pick <- ref_new_ranked[ref_new_ranked %in% available_new_labels & !(ref_new_ranked %in% used_new_labels)]

      if (length(pick) > 0) {
        method_to_new[as.character(k)] <- pick[1]
        used_new_labels <- c(used_new_labels, pick[1])
      }
    }

    # assign extra unmatched clusters to remaining colors
    remaining_method <- names(method_to_new)[is.na(method_to_new)]
    remaining_colors <- setdiff(available_new_labels, used_new_labels)

    if (length(remaining_method) > length(remaining_colors)) {
      stop("Not enough colors in cluster_cols to assign all extra clusters.")
    }

    if (length(remaining_method) > 0) {
      method_to_new[remaining_method] <- remaining_colors[seq_along(remaining_method)]
    }

    mapped <- unname(method_to_new[as.character(method_vec)])
    names(mapped) <- names(method_vec)

    list(
      mapped_cluster = mapped,
      mapping = method_to_new,
      overlap = overlap_tab
    )
  }

  plot_df <- data.frame(
    study = study_ids,
    stringsAsFactors = FALSE
  )
  plot_df[[ref_name]] <- ref_reordered

  if (length(other) > 0) {
    if (is.null(names(other)) || any(names(other) == "")) {
      stop("other must be a named list.")
    }

    for (nm in names(other)) {
      tmp <- map_to_ref_with_extra(
        method_vec = other[[nm]],
        ref_vec = ref_vec,
        ref_old_to_new = ref_old_to_new,
        cluster_cols = cluster_cols
      )
      plot_df[[nm]] <- tmp$mapped_cluster[plot_df$study]
    }
  }

  plot_df$study_rank <- study_rank[plot_df$study]
  plot_df <- plot_df[order(plot_df[[ref_name]], plot_df$study_rank), , drop = FALSE]
  plot_df$study <- factor(plot_df$study, levels = rev(plot_df$study))

  method_names <- c(ref_name, names(other))
  if (is.null(method_order)) {
    method_order <- order_methods(method_names)
  } else {
    method_order <- c(intersect(method_order, method_names),
                      setdiff(method_names, method_order))
  }

  long_df <- tidyr::pivot_longer(
    plot_df,
    cols = all_of(method_names),
    names_to = "method",
    values_to = "cluster"
  )

  ## y is drawn bottom-up, so reverse to get method_order top-to-bottom.
  long_df$method  <- factor(long_df$method, levels = rev(method_order))
  long_df$cluster <- factor(long_df$cluster, levels = names(cluster_cols))

  lab_y <- method_label_expr(method_order, ref_name)

  ggplot2::ggplot(long_df, ggplot2::aes(y = method, x = study, fill = cluster)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.2, width = 1, height = height_gap) +
    ggplot2::scale_y_discrete(
      labels = function(lbl) {
        sapply(lbl, function(z) parse(text = lab_y[[as.character(z)]])[[1]])
      }
    ) +
    labs(title = "Comparison of clustering assignments") +
    ggplot2::scale_fill_manual(values = cluster_cols, drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10)
    )
}


plot_shared_specific <- function(shared_lst,
                                 specific_lst,
                                 method_order = names(shared_lst),
                                 col_shared = "#222222",
                                 col_specific = "#9AA7B1") {

  # count features
  df <- data.frame(
    method = method_order,
    Shared = sapply(method_order, function(m) length(shared_lst[[m]])),
    `Cluster-specific` = sapply(method_order, function(m) length(specific_lst[[m]])),
    check.names = FALSE
  )

  # keep chosen order
  df$method <- factor(df$method, levels = rev(method_order))

  # left panel: shared
  p1 <- ggplot(df, aes(x = Shared, y = method)) +
    geom_col(fill = col_shared, width = 0.62) +
    labs(title = "Shared signatures") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.grid = element_blank(),
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
      plot.margin = margin(5.5, 10, 5.5, 5.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  # right panel: cluster-specific
  p2 <- ggplot(df, aes(x = `Cluster-specific`, y = method)) +
    geom_col(fill = col_specific, width = 0.62) +
    labs(title = "Cluster-specific signatures") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  p1 + p2 + plot_layout(widths = c(1, 1))
}


plot_jaccard_shared_specific <- function(shared_lst,
                                         specific_lst,
                                         ref_method = "SMESH-PALM",
                                         method_order = NULL,
                                         cols = c("#F2F2F2", "#9ECAE1", "#3182BD")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (!(ref_method %in% names(shared_lst)) || !(ref_method %in% names(specific_lst))) {
    stop("ref_method must exist in both shared_lst and specific_lst.")
  }

  other_methods <- intersect(names(shared_lst), names(specific_lst))

  if (is.null(method_order)) {
    method_order <- other_methods
  } else {
    method_order <- method_order[method_order %in% other_methods]
  }

  df <- data.frame(
    method = rep(method_order, 2),
    type = rep(c("Shared signatures", "Cluster-specific signatures"), each = length(method_order)),
    jaccard = c(
      sapply(method_order, function(m) {
        jaccard_index(shared_lst[[ref_method]], shared_lst[[m]])
      }),
      sapply(method_order, function(m) {
        jaccard_index(specific_lst[[ref_method]], specific_lst[[m]])
      })
    ),
    stringsAsFactors = FALSE
  )

  df$method <- factor(df$method, levels = rev(method_order))
  df$type <- factor(df$type, levels = c("All-cluster shared", "Cluster-specific"))
  df$label <- ifelse(is.na(df$jaccard), "NA", sprintf("%.2f", df$jaccard))

  ggplot2::ggplot(df, ggplot2::aes(x = type, y = method, fill = jaccard)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.8, width = 1, height = 1) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 6, fontface = "bold") +
    ggplot2::scale_fill_gradientn(
      colours = cols,
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Jaccard overlap",
      na.value = "white"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

#' Signature counts and Jaccard overlap against a reference method.
plot_panel_E_bar <- function(shared_lst,
                             specific_lst,
                             ref_method = "SMESH-PALM",
                             method_order = c(
                               "SMESH-PALM",
                               "SMESH-ANCOMBC2",
                               "SMESH-LinDA",
                               "SMESH-MaAsLin3",
                               "SKM+FE",
                               "SHC+FE",
                               "Melody"
                             ),
                             col_shared = "#222222",
                             col_specific = "#8FA3B0") {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)

  df <- data.frame(
    method = method_order,
    shared_count = sapply(method_order, function(m) length(shared_lst[[m]])),
    specific_count = sapply(method_order, function(m) length(specific_lst[[m]])),
    shared_jaccard = sapply(method_order, function(m) jaccard_index(shared_lst[[ref_method]], shared_lst[[m]])),
    specific_jaccard = sapply(method_order, function(m) jaccard_index(specific_lst[[ref_method]], specific_lst[[m]])),
    stringsAsFactors = FALSE
  )

  df$method <- factor(df$method, levels = rev(method_order))

  lab_y <- method_label_expr(method_order, ref_method)

  df_count <- df %>%
    dplyr::select(method, shared_count, specific_count) %>%
    pivot_longer(
      cols = c(shared_count, specific_count),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(
      type = factor(type,
                    levels = c("specific_count", "shared_count"),
                    labels = c("Context-dependent", "All-cluster shared"))
    )

  df_jaccard <- df %>%
    dplyr::select(method, shared_jaccard, specific_jaccard) %>%
    pivot_longer(
      cols = c(shared_jaccard, specific_jaccard),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(
      type = factor(type,
                    levels = c("specific_jaccard", "shared_jaccard"),
                    labels = c("Context-dependent", "All-cluster shared")),
      label = sprintf("%.2f", value)
    )

  xmax_count <- max(df_count$value, na.rm = TRUE)

  p_count <- ggplot(df_count, aes(x = value, y = method, fill = type)) +
    geom_col(width = 0.62, position = position_dodge(width = 0.75)) +
    geom_text(
      aes(label = ifelse(value == 0, "", value)),
      position = position_dodge(width = 0.75),
      hjust = -0.15,
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "All-cluster shared" = col_shared,
        "Context-dependent" = col_specific
      ),
      name = NULL
    ) +
    scale_y_discrete(
      labels = function(lbl) {
        sapply(lbl, function(z) parse(text = lab_y[[as.character(z)]])[[1]])
      }
    ) +
    scale_x_continuous(
      limits = c(0, xmax_count * 1.18),
      breaks = pretty(c(0, xmax_count)),
      labels = function(x) abs(x),
      expand = c(0, 0)
    ) +
    labs(title = "Number of selected signatures") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 19),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5)
    )

  p_jaccard <- ggplot(df_jaccard, aes(x = value, y = method, fill = type)) +
    geom_col(width = 0.62, position = position_dodge(width = 0.75)) +
    geom_text(
      aes(label = ifelse(value == 0, "", label)),
      position = position_dodge(width = 0.75),
      hjust = -0.15,
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "All-cluster shared" = col_shared,
        "Context-dependent" = col_specific
      ),
      breaks = c("All-cluster shared", "Context-dependent"),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = c(0, 0)
    ) +
    labs(title = paste0("Jaccard overlap (vs ", ref_method, ")")) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.text = element_text(size = 18),
      plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 12)
    )

  p_count + p_jaccard + plot_layout(widths = c(1.55, 1.55))
}

#' Selected number of clusters, all-data value against the leave-one-context-out
#' spread.
#'
#' `df$type` is "ALL" for the all-data fit; every other value is treated as a
#' leave-one-context-out run and labelled "LOCO" (older scripts wrote "LODO").
plot_G_box <- function(df, method_order = NULL) {
  library(dplyr)
  library(ggplot2)

  df_loco <- df %>%
    filter(type != "ALL") %>%
    mutate(
      Method = as.character(Method),
      G = as.numeric(G),
      point_type = "LOCO"
    )

  df_all <- df %>%
    filter(type == "ALL") %>%
    mutate(
      Method = as.character(Method),
      G = as.numeric(G),
      point_type = "ALL"
    ) %>%
    distinct(Method, .keep_all = TRUE)

  if (is.null(method_order)) {
    method_order <- unique(df$Method)
  }

  df_loco$Method <- factor(df_loco$Method, levels = method_order)
  df_all$Method  <- factor(df_all$Method, levels = method_order)

  ggplot(df_loco, aes(x = Method, y = G)) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.35,
      outlier.shape = NA
    ) +
    geom_jitter(
      aes(shape = point_type, fill = point_type),
      width = 0.12,
      height = 0,
      size = 2.4,
      alpha = 0.8,
      color = "black"
    ) +
    geom_point(
      data = df_all,
      aes(x = Method, y = G, shape = point_type, fill = point_type),
      inherit.aes = FALSE,
      size = 4,
      stroke = 1.1,
      color = "black"
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("ALL" = 23, "LOCO" = 16)
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("ALL" = "gold", "LOCO" = "black")
    ) +
    guides(
      fill = "none",
      shape = guide_legend(
        override.aes = list(
          shape = c(23, 16),
          fill  = c("gold", "black"),
          color = c("black", "black"),
          size  = c(4, 2.4)
        )
      )
    ) +
    scale_y_continuous(breaks = sort(unique(df$G))) +
    labs(
      x = NULL,
      y = NULL
    ) +
    ggtitle("Selected number of clusters (G)") +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 13),
      plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.margin = margin(5.5, 10, 5.5, 10, unit = "mm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}


# =============================================================================
#  6. Cluster-stability helpers
# =============================================================================

#' Relabel clusters C1, C2, ... from largest to smallest.
relabel_cluster_by_size <- function(df) {
  cl_map <- df %>%
    dplyr::count(cluster_raw, name = "n") %>%
    dplyr::arrange(desc(n), cluster_raw) %>%
    dplyr::mutate(cluster_label = paste0("C", dplyr::row_number())) %>%
    dplyr::select(cluster_raw, cluster_label)

  df %>% dplyr::left_join(cl_map, by = "cluster_raw")
}

#' Co-clustering frequency across leave-one-out runs.
#'
#' @param cluster_list list of named cluster vectors, one per run.  Studies held
#'                     out of a run are simply absent from its vector and are
#'                     excluded from that run's denominator.
build_consensus_lodo <- function(cluster_list) {
  ids <- sort(unique(unlist(lapply(cluster_list, names))))
  n <- length(ids)

  agree_mat <- matrix(0, n, n, dimnames = list(ids, ids))
  count_mat <- matrix(0, n, n, dimnames = list(ids, ids))

  for (cl in cluster_list) {
    cl_full <- setNames(rep(NA_integer_, n), ids)
    cl_full[names(cl)] <- cl

    obs <- !is.na(cl_full)

    ## pairs observed in this run
    obs_pair <- outer(obs, obs, FUN = "&")

    ## same-cluster indicator only for observed pairs
    same_pair <- outer(cl_full, cl_full, FUN = "==")
    same_pair[!obs_pair] <- FALSE

    agree_mat <- agree_mat + same_pair
    count_mat <- count_mat + obs_pair
  }

  C <- agree_mat / count_mat
  diag(C) <- 1

  list(C = C, count_mat = count_mat, agree_mat = agree_mat)
}

jaccard_index <- function(x, y) {
  x <- unique(x)
  y <- unique(y)
  u <- union(x, y)
  if (length(u) == 0) return(NA_real_)
  length(intersect(x, y)) / length(u)
}

#' Share of studies keeping their cluster membership between two clusterings.
#'
#' Cluster labels are matched one-to-one by the Hungarian algorithm before
#' comparing, so a pure relabelling counts as full retention.
retain_membership <- function(vec1, vec2) {
  stopifnot(length(vec1) == length(vec2))

  # if named, require same names (same studies)
  if (!is.null(names(vec1)) || !is.null(names(vec2))) {
    stopifnot(identical(names(vec1), names(vec2)))
  }

  f1 <- factor(vec1)
  f2 <- factor(vec2)

  tab <- table(f1, f2)  # overlap counts (can be rectangular)

  nr <- nrow(tab); nc <- ncol(tab)
  n  <- max(nr, nc)

  # Pad to square with zeros
  tab_sq <- matrix(0, n, n)
  tab_sq[1:nr, 1:nc] <- as.matrix(tab)

  # nonnegative cost
  cost <- max(tab_sq) - tab_sq

  # Hungarian assignment: one-to-one between rows and cols (including dummy padded ones)
  assignment <- clue::solve_LSAP(cost)  # assignment[i] = chosen column for row i

  # Build mapping only for REAL columns (original vec2 clusters)
  real_row_names <- rownames(tab)
  real_col_names <- colnames(tab)

  map <- rep(NA_character_, length(real_col_names))
  names(map) <- real_col_names

  for (i in seq_len(nr)) {
    j <- as.integer(assignment[i])
    if (j <= nc) {
      # vec1 cluster (row i) matched to vec2 cluster (col j)
      map[j] <- real_row_names[i]
    }
  }

  # Relabel vec2 using mapping; unmatched vec2 clusters become NA
  vec2_aligned <- map[as.character(f2)]
  names(vec2_aligned) <- names(vec2)

  # Per-study 0/1 (NA mapping => 0 by default)
  keep <- as.integer(as.character(f1) == vec2_aligned)
  keep[is.na(keep)] <- 0L
  names(keep) <- names(vec1)

  list(
    retain_rate = mean(keep),
    keep_status = keep,
    mapping_vec2_to_vec1 = map,
    overlap_table = tab
  )
}

#' Map another method's clusters onto the reference clusters, in display order.
map_cluster_to_ref_order <- function(ref_cluster,
                                     tag_cluster,
                                     ref_order = sort(unique(ref_cluster))) {
  if (is.null(names(ref_cluster)) || is.null(names(tag_cluster))) {
    stop("Both 'ref_cluster' and 'tag_cluster' must be named vectors, with study IDs as names.")
  }

  common_ids <- intersect(names(ref_cluster), names(tag_cluster))
  if (length(common_ids) == 0) {
    stop("No overlapping study IDs found between 'ref_cluster' and 'tag_cluster'.")
  }

  ref_cluster <- ref_cluster[common_ids]
  tag_cluster <- tag_cluster[common_ids]

  smesh_labels <- sort(unique(ref_cluster))

  if (length(ref_order) != length(smesh_labels)) {
    stop("'ref_order' must have the same length as the number of unique reference clusters.")
  }

  if (!setequal(ref_order, smesh_labels)) {
    stop("'ref_order' must contain exactly the reference cluster labels, just in the desired order.")
  }

  overlap_tab <- table(ref = ref_cluster, tag = tag_cluster)
  ordered_tab <- overlap_tab[as.character(ref_order), , drop = FALSE]

  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required. Please install it with install.packages('clue').")
  }

  ## assignment[i] is the tag column matched to reference row i, so the mapping
  ## has to be filled by scattering ref_order into those columns - indexing
  ## `ref_order[assignment]` instead silently applies the inverse permutation.
  assignment <- clue::solve_LSAP(max(ordered_tab) - ordered_tab)

  # tag cluster j -> reference cluster
  mapping <- setNames(rep(NA, ncol(ordered_tab)), colnames(ordered_tab))
  mapping[as.integer(assignment)] <- ref_order

  mapped_tag <- unname(mapping[as.character(tag_cluster)])
  names(mapped_tag) <- names(tag_cluster)

  list(
    overlap_table = overlap_tab,
    ordered_overlap_table = ordered_tab,
    mapping = mapping,
    mapped_tag_cluster = mapped_tag,
    mapped_overlap_table = table(ref = ref_cluster, tag_mapped = mapped_tag)
  )
}

#' Columns of a competing model's effect matrix, in reference display order.
#'
#' Answers "which column of the other method's `mu` plays the role of the
#' reference cluster shown at display position k?", so its heat map can be
#' stacked row-for-row against the reference one.
#'
#' @param ref_cluster,tag_cluster named cluster vectors over the same studies.
#' @param cluster_order reference cluster labels, in display order.
#'
#' @return tag column indices, one per entry of `cluster_order`.
tag_cols_in_ref_order <- function(ref_cluster, tag_cluster, cluster_order) {
  mp <- map_cluster_to_ref_order(
    ref_cluster = ref_cluster,
    tag_cluster = tag_cluster,
    ref_order   = sort(unique(ref_cluster))
  )$mapping

  inv <- setNames(names(mp), mp)          # reference cluster -> tag cluster
  as.integer(inv[as.character(cluster_order)])
}
