# =============================================================================
#   Make_TableS4_S6.R
#   SMESH-PALM cluster-level effect estimates for the three real data
#   applications, as Supplementary Tables S4, S5 and S6.
# =============================================================================
#  Writes ONE workbook with exactly three tabs:
#
#    "Supplementary table 4"   pan-disease  (GMrepo, 12 diseases,     G = 2)
#    "Supplementary table 5"   pan-tumor    (HGMT,   15 tumour types, G = 4)
#    "Supplementary table 6"   colorectal   (26 study-stage strata,   G = 2)
#
#  Each tab has the same three column groups:
#
#    Feature ID    the modelled taxon, in the order it is drawn along the
#                  x axis of panel A of the corresponding main figure
#                  (Figure 4, 5 and 6)
#    Shared type   which clusters select the feature: all of them, a proper
#                  subset, or exactly one
#    C1 ... CG     the cluster-level effect estimate, one column per cluster
#
#  CONVENTIONS
#   * The reported fit is the one used in the main figures: the G selected by
#     the generalized information criterion on the full data (G = 2, 4 and 2).
#   * Clusters are numbered in DISPLAY order (C1, C2, ...), the order used in
#     the figures, not by the raw mixture-component index.
#   * Row order reproduces panel A exactly.  Features are grouped by which
#     clusters select them - all-cluster shared first, then the partially
#     shared subsets, then the cluster-specific ones - and ordered inside a
#     group by the sum of their cluster effects.  This is
#     `build_signature_groups()`, the same call the figure scripts make, so the
#     n-th row of a tab is the n-th column of the heat map.
#   * Only selected signatures are listed: a feature with no non-zero cluster
#     effect belongs to no signature group and is drawn in no panel.
#   * A cluster in which no context could estimate a feature is reported as the
#     literal string NA, not as a blank and not as 0, so "no data in this
#     cluster" stays distinguishable from "tested, effect zero".  These are the
#     cells the figures leave white.
#   * SMESH returns point estimates for the cluster-level effects; the fitted
#     object carries no standard errors, so none are reported.
# =============================================================================

  rm(list = ls())
  suppressMessages({ library(dplyr); library(tibble) })
  source("./utility/heatmap_util.R")

  OUT <- "./Figure/TableS4_S6_SMESH_effects.xlsx"
  dir.create("./Figure", showWarnings = FALSE, recursive = TRUE)

  lo <- function(f, nm) { e <- new.env(); load(f, envir = e); get(nm, envir = e) }

  APPS <- list(
    list(sheet = "Supplementary table 4", app = "Pan-disease (GMrepo)",
         loso = "./GMrepo_analysis/GMrepo_loso", data = "./GMrepo_analysis/Data",
         ss = "summary_stat_meta_filter", G = 2),
    list(sheet = "Supplementary table 5", app = "Pan-tumor (HGMT)",
         loso = "./HGMT_analysis/HGMT_loso", data = "./HGMT_analysis/Data",
         ss = "summary_stat_meta_filter", G = 4),
    list(sheet = "Supplementary table 6", app = "Colorectal neoplasia",
         loso = "./CRC_analysis/CRC_loso", data = "./CRC_analysis/Data",
         ss = "summary_stat_filter", G = 2)
  )

  ## Group name from build_signature_groups() ("1|3") -> reader-facing label.
  shared_type <- function(nm, G) {
    pos <- strsplit(nm, "|", fixed = TRUE)[[1]]
    if (length(pos) == G) "All-cluster shared"
    else if (length(pos) == 1) sprintf("Cluster-specific (C%s)", pos)
    else sprintf("Partially shared (%s)", paste0("C", pos, collapse = "|"))
  }

  sheets <- list()

  for (a in APPS) {
    ss  <- lo(file.path(a$data, "Summary_stat_SMESH.Rdata"), a$ss)
    mdl <- lo(file.path(a$loso, sprintf("Model%d_SMESH_s100.Rdata", a$G)), "SMESH_model")

    W  <- round(mdl$disease$W)                      # context x cluster
    mu <- mdl$disease$mu                            # feature x cluster
    st <- summary_stat_matrices(ss, rownames(mu))
    G  <- ncol(W)

    ## Same two calls the figure scripts make, in the same order.
    beta          <- mask_unobserved_clusters(mu, W, st$est)
    cluster_order <- order_clusters_by_selection(beta)
    species_lst   <- build_signature_groups(beta, cluster_order)

    ## Panel A draws `unlist(species_lst)` left to right.
    taxa <- unlist(species_lst, use.names = FALSE)
    type <- rep(vapply(names(species_lst), shared_type, character(1), G = G),
                times = lengths(species_lst))

    ## Which (feature, cluster) pairs no member context could estimate.
    estimable <- vapply(seq_len(G), function(g) {
      members <- rownames(W)[W[, g] != 0]
      rowSums(!is.na(st$est[, members, drop = FALSE])) > 0
    }, logical(nrow(mu)))

    est <- beta[taxa, cluster_order, drop = FALSE]
    est[!estimable[taxa, cluster_order, drop = FALSE]] <- NA_real_
    colnames(est) <- paste0("C", seq_len(G))

    sheets[[a$sheet]] <- tibble(`Feature ID` = taxa, `Shared type` = type) %>%
      bind_cols(as_tibble(round(est, 4)))

    cat(sprintf("%-22s %2d contexts, G = %d, %3d signatures, %4d NA effects\n",
                a$app, nrow(W), G, length(taxa), sum(is.na(est))))
  }

  ## keepNA/na.string write the string "NA" rather than an empty cell.
  openxlsx::write.xlsx(sheets, file = OUT, overwrite = TRUE,
                       keepNA = TRUE, na.string = "NA",
                       headerStyle = openxlsx::createStyle(textDecoration = "bold"),
                       colWidths = "auto")

  cat("\nwritten:", OUT, "\n")
  for (s in names(sheets))
    cat(sprintf("  %-24s %3d rows x %d cols\n", s, nrow(sheets[[s]]), ncol(sheets[[s]])))
