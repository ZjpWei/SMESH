# =============================================================================
#   Check_cluster_covariates.R
#   Are the learned SMESH clusters explained by context-level demographic,
#   study-design or technical characteristics?
# =============================================================================
#  This is a supporting check, NOT a supplementary table.  It backs the claim
#  made three times in the main text that the learned clusters are not a
#  restatement of context-level characteristics:
#     p.9  L200-202  pan-disease   (12 contexts, 2 clusters)
#     p.11 L265-266  pan-tumor     (15 contexts, 4 clusters)
#     p.13 L315-317  colorectal    (26 contexts, 2 clusters)
#  Nothing it writes is shipped with the paper.
#
#  DESIGN
#    unit      : the CONTEXT, not the sample
#    outcome   : SMESH cluster label (display order)
#    statistic : numeric     -> Kruskal-Wallis
#                categorical -> Pearson chi-square
#    null      : permutation of the cluster labels, which preserves the cluster
#                sizes.  With two clusters the number of distinct labellings is
#                small enough to ENUMERATE EXHAUSTIVELY (exact p-value); with
#                more clusters it is not, so a Monte Carlo permutation p-value
#                is used instead and the number of draws is reported.
#    multiple  : Benjamini-Hochberg within each application
#
#  Sample-level fields are aggregated to the context on DISTINCT samples first,
#  so a control shared between two contrasts is not counted twice.
# =============================================================================

  rm(list = ls())
  suppressMessages({ library(dplyr); library(tidyr) })
  source("./utility/heatmap_util.R")

  lo <- function(f, nm) { e <- new.env(); load(f, envir = e); get(nm, envir = e) }

  N_PERM  <- 20000     # Monte Carlo draws when exhaustive enumeration is infeasible
  OUT_DIR <- "./Figure"
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  blank_to_na <- function(v) {
    v <- as.character(v)
    v[v %in% c("", " ", "NA", "na", "N/A", "unknown", "Unknown", "not available", "-")] <- NA
    v
  }


# -----------------------------------------------------------------------------
#  cluster labels in DISPLAY order, for one application
# -----------------------------------------------------------------------------

  display_clusters <- function(loso_dir, data_dir, ss_name, G) {
    ss   <- lo(file.path(data_dir, "Summary_stat_SMESH.Rdata"), ss_name)
    mdl  <- lo(file.path(loso_dir, sprintf("Model%d_SMESH_s100.Rdata", G)), "SMESH_model")
    W    <- round(mdl$disease$W)
    st   <- summary_stat_matrices(ss, rownames(mdl$disease$mu))
    beta <- mask_unobserved_clusters(mdl$disease$mu, W, st$est)
    co   <- order_clusters_by_selection(beta)
    raw  <- setNames(apply(mdl$disease$W, 1, which.max), rownames(mdl$disease$W))
    list(cluster = setNames(paste0("C", match(raw, co)), names(raw)),
         n_feature = vapply(ss, function(s) sum(!is.na(s$est[, 1])), integer(1)))
  }


# -----------------------------------------------------------------------------
#  context-level characteristic tables, one builder per application
# -----------------------------------------------------------------------------

  build_gmrepo <- function() {
    cl   <- display_clusters("./GMrepo_analysis/GMrepo_loso", "./GMrepo_analysis/Data",
                             "summary_stat_meta_filter", 2)
    meta <- readRDS("GMrepo_analysis/Data/metadata_list_final.rds")
    ci   <- readRDS("GMrepo_analysis/Data/covariate.interest.rds")
    psum <- openxlsx::read.xlsx("GMrepo_analysis/Data/GMrepo_project_summary.xlsx")
    m2n  <- setNames(psum$disease_name, psum$disease_mesh)

    d <- bind_rows(lapply(names(meta), function(k) {
      z <- as.data.frame(meta[[k]], stringsAsFactors = FALSE)
      v <- ci[[k]]
      data.frame(context = unname(m2n[sub("^.+_", "", k)]),
                 study   = sub("_[^_]+$", "", k),
                 id      = as.character(z$run_id),
                 is_case = v[match(z$run_id, rownames(v)), 1],
                 sex     = blank_to_na(z$sex),
                 age     = suppressWarnings(as.numeric(blank_to_na(z$host_age))),
                 country = blank_to_na(z$country),
                 platform = blank_to_na(z$instrument_model),
                 depth   = suppressWarnings(as.numeric(blank_to_na(z$nr_reads_sequenced))),
                 stringsAsFactors = FALSE)
    })) %>% distinct(context, id, .keep_all = TRUE)
    list(d = d, cl = cl, app = "pan-disease (GMrepo)")
  }

  build_hgmt <- function() {
    cl <- display_clusters("./HGMT_analysis/HGMT_loso", "./HGMT_analysis/Data",
                           "summary_stat_meta_filter", 4)
    e  <- new.env(); load("HGMT_analysis/Data/HGMT_16S.Rdata", envir = e)
    otu <- e$otu_final; mt <- e$meta_final

    d <- bind_rows(lapply(names(otu), function(k) {
      md <- mt[[k]]
      ctx <- sub("^PRJ[A-Z0-9-]+(?:_(?:16S|WGS|PE|SE))?_", "", k)
      lib <- rowSums(as.matrix(otu[[k]]), na.rm = TRUE)
      data.frame(context = ctx,
                 study   = as.character(md$`Project ID`),
                 id      = as.character(md$`Run ID`),
                 is_case = as.integer(md$`Phenotype name` != "Health"),
                 sex     = blank_to_na(md$sex),
                 age     = suppressWarnings(as.numeric(blank_to_na(md$age))),
                 country = blank_to_na(md$country),
                 platform = blank_to_na(md$`Sequencing method`),
                 depth   = unname(lib[as.character(md$`Run ID`)]),
                 stringsAsFactors = FALSE)
    })) %>% distinct(context, id, .keep_all = TRUE)
    list(d = d, cl = cl, app = "pan-tumor (HGMT)")
  }

  build_crc <- function() {
    cl <- display_clusters("./CRC_analysis/CRC_loso", "./CRC_analysis/Data",
                           "summary_stat_filter", 2)
    md <- lo("CRC_analysis/Data/CRC_strata.Rdata", "meta_data_sp")
    ctxs <- names(cl$cluster)
    d <- bind_rows(lapply(ctxs, function(k) {
      z <- as.data.frame(md[[k]], stringsAsFactors = FALSE)
      data.frame(context = k,
                 study   = sub(":.*$", "", k),
                 id      = as.character(z$SampleID),
                 is_case = as.integer(z$Group),
                 sex     = blank_to_na(z$Sex),
                 age     = suppressWarnings(as.numeric(blank_to_na(z$Age))),
                 country = NA_character_,
                 platform = NA_character_,
                 depth   = NA_real_,
                 stringsAsFactors = FALSE)
    })) %>% distinct(context, id, .keep_all = TRUE)
    list(d = d, cl = cl, app = "colorectal neoplasia")
  }


# -----------------------------------------------------------------------------
#  aggregate to the context
# -----------------------------------------------------------------------------

  summarise_contexts <- function(b) {
    female <- function(s) {
      s <- tolower(s); s[!s %in% c("female", "male", "f", "m")] <- NA
      100 * mean(s %in% c("female", "f"), na.rm = TRUE)
    }
    ctx <- b$d %>%
      group_by(context) %>%
      summarise(
        n_studies    = n_distinct(study),
        n_samples    = n(),
        n_case       = sum(is_case == 1, na.rm = TRUE),
        n_control    = sum(is_case == 0, na.rm = TRUE),
        case_frac    = n_case / (n_case + n_control),
        pct_female   = female(sex),
        median_age   = median(age, na.rm = TRUE),
        sex_reported = 100 * mean(!is.na(sex)),
        age_reported = 100 * mean(!is.na(age)),
        median_depth = median(depth, na.rm = TRUE),
        n_countries  = ifelse(all(is.na(country)), NA_integer_, n_distinct(na.omit(country))),
        top_country  = ifelse(all(is.na(country)), NA_character_,
                              names(sort(table(country), decreasing = TRUE))[1]),
        n_platforms  = ifelse(all(is.na(platform)), NA_integer_, n_distinct(na.omit(platform))),
        .groups = "drop") %>%
      mutate(n_feature = unname(b$cl$n_feature[context]),
             cluster   = unname(b$cl$cluster[context])) %>%
      arrange(cluster, context)
    ctx[!is.nan(ctx$pct_female) | TRUE, ]
  }


# -----------------------------------------------------------------------------
#  permutation test
# -----------------------------------------------------------------------------

  MAX_EXACT <- 200000   # enumerate exhaustively only below this many labellings

  #' Permutation test of one characteristic against the cluster label.
  #'
  #' Permuting the labels is equivalent to permuting the values, so the null is
  #' built by shuffling the values against a FIXED cluster-indicator matrix -
  #' which turns the whole null distribution into one matrix product for numeric
  #' variables.  Two clusters with few enough labellings are enumerated exactly;
  #' otherwise a Monte Carlo p-value is reported and `mode` says so.
  #'
  #' Numeric     : Kruskal-Wallis H.  The tie correction is a constant factor
  #'               that does not depend on the labels, so it is omitted - it
  #'               cannot change a permutation p-value.
  #' Categorical : Pearson chi-square.  Row and column margins are fixed under
  #'               permutation, so the expected counts are computed once.
  perm_test <- function(x, grp, n_perm = N_PERM, seed = 1) {
    keep <- !is.na(x) & !(is.numeric(x) & is.nan(suppressWarnings(as.numeric(x))))
    x <- x[keep]; grp <- droplevels(factor(grp[keep]))
    n <- length(x); G <- nlevels(grp)
    if (G < 2 || n < 4 || length(unique(x)) < 2)
      return(list(p = NA_real_, n = n, mode = "n/a", n_null = NA_integer_))

    M  <- model.matrix(~ grp - 1)          # n x G cluster indicators (fixed)
    ns <- colSums(M)

    ## how the null is generated: matrix of permutation index columns
    if (G == 2 && choose(n, min(ns)) <= MAX_EXACT) {
      combos <- combn(n, ns[1])
      idx <- apply(combos, 2, function(k) c(k, setdiff(seq_len(n), k)))
      mode <- "exact"; n_null <- ncol(idx)
    } else {
      set.seed(seed)
      idx <- replicate(n_perm, sample.int(n))
      mode <- "MC"; n_null <- n_perm
    }

    if (is.numeric(x)) {
      r    <- rank(x)
      obs  <- 12 / (n * (n + 1)) * sum(colSums(M * r)^2 / ns) - 3 * (n + 1)
      S    <- crossprod(M, matrix(r[idx], nrow = n))          # G x n_null
      null <- 12 / (n * (n + 1)) * colSums(S^2 / ns) - 3 * (n + 1)
    } else {
      xf <- factor(x); L <- nlevels(xf); xi <- as.integer(xf); gi <- as.integer(grp)
      E  <- outer(as.vector(table(xf)), ns) / n
      chi2 <- function(v) sum(matrix(tabulate((v - 1L) * G + gi, L * G),
                                     nrow = L, byrow = TRUE)^2 / E) - n
      obs  <- chi2(xi)
      null <- apply(idx, 2, function(k) chi2(xi[k]))
    }

    p <- if (mode == "exact") mean(null >= obs - 1e-9)
         else (1 + sum(null >= obs - 1e-9)) / (n_null + 1)
    list(p = p, n = n, mode = mode, n_null = n_null)
  }

  VARS <- c("n_studies", "n_samples", "case_frac",
            "pct_female", "median_age", "sex_reported", "age_reported",
            "median_depth", "n_countries", "n_platforms", "top_country", "n_feature")

  ## Human-readable labels and the group each variable belongs to, used for the
  ## published table.
  VAR_INFO <- tibble::tribble(
    ~variable,       ~group,         ~label,
    "n_studies",     "Study design", "Number of contributing studies",
    "n_samples",     "Study design", "Number of samples",
    "case_frac",     "Study design", "Case fraction (cases / total)",
    "pct_female",    "Demographic",  "Percentage female",
    "median_age",    "Demographic",  "Median age (years)",
    "sex_reported",  "Technical",    "Percentage of samples with sex recorded",
    "age_reported",  "Technical",    "Percentage of samples with age recorded",
    "median_depth",  "Technical",    "Median library size",
    "n_countries",   "Technical",    "Number of countries represented",
    "n_platforms",   "Technical",    "Number of sequencing platforms",
    "top_country",   "Technical",    "Most frequent country",
    "n_feature",     "Technical",    "Number of estimable features")


# -----------------------------------------------------------------------------
#  run all three applications
# -----------------------------------------------------------------------------

  all_ctx <- list(); all_res <- list()

  for (b in list(build_gmrepo(), build_hgmt(), build_crc())) {
    ctx <- summarise_contexts(b)

    cat("\n", strrep("=", 96), "\n", sep = "")
    cat("  ", b$app, "  (", nrow(ctx), " contexts, ",
        length(unique(ctx$cluster)), " clusters: ",
        paste(names(table(ctx$cluster)), table(ctx$cluster), sep = "=", collapse = " "), ")\n", sep = "")
    cat(strrep("=", 96), "\n")
    print(as.data.frame(ctx %>% select(context, cluster, n_studies, n_samples, case_frac,
                                       pct_female, median_age, n_countries, n_feature)),
          row.names = FALSE, digits = 3)

    res <- bind_rows(lapply(VARS, function(v) {
      if (all(is.na(ctx[[v]]))) return(NULL)
      r <- perm_test(ctx[[v]], ctx$cluster)
      tibble(application = b$app, variable = v,
             type = ifelse(is.numeric(ctx[[v]]), "numeric", "categorical"),
             n_ctx = r$n, mode = r$mode, n_null = r$n_null, p = r$p)
    })) %>% mutate(p_BH = p.adjust(p, method = "BH"))

    cat("\n  permutation tests:\n")
    print(as.data.frame(res %>% select(-application)), row.names = FALSE, digits = 3)
    cat(sprintf("\n  min p = %.4f   min adjusted p = %.4f   variables with p_BH < 0.1 : %d\n",
                min(res$p, na.rm = TRUE), min(res$p_BH, na.rm = TRUE),
                sum(res$p_BH < 0.1, na.rm = TRUE)))

    all_ctx[[b$app]] <- ctx %>% mutate(application = b$app, .before = 1)
    all_res[[b$app]] <- res
  }

# -----------------------------------------------------------------------------
#  write one workbook, one tab per piece of evidence
# -----------------------------------------------------------------------------

  res_all <- bind_rows(all_res) %>%
    left_join(VAR_INFO, by = "variable") %>%
    transmute(Application = application,
              Group       = group,
              Characteristic = label,
              Type        = type,
              N_contexts  = n_ctx,
              Null        = dplyr::case_when(
                              mode == "exact" ~ paste0("exact (", n_null, " labellings)"),
                              mode == "MC"    ~ paste0("Monte Carlo (", n_null, " permutations)"),
                              TRUE            ~ "not testable (constant across contexts)"),
              P           = round(p, 4),
              P_adjusted  = round(p_BH, 4))

  ctx_tab <- function(app) {
    d <- all_ctx[[app]]
    d %>% transmute(Context = context, Cluster = cluster,
                    N_studies = n_studies, N_samples = n_samples,
                    N_case = n_case, N_control = n_control,
                    Case_fraction = round(case_frac, 3),
                    Pct_female = round(pct_female, 1),
                    Median_age = median_age,
                    Pct_sex_recorded = round(sex_reported, 1),
                    Pct_age_recorded = round(age_reported, 1),
                    Median_library_size = round(median_depth),
                    N_countries = n_countries, Top_country = top_country,
                    N_platforms = n_platforms, N_estimable_features = n_feature)
  }

  readme <- tibble::tibble(Item = c(
    "Title",
    "Question",
    "Unit of analysis",
    "Outcome",
    "Statistic (numeric)",
    "Statistic (categorical)",
    "Null distribution",
    "Multiple testing",
    "Note on power",
    "Note on missing data",
    "Generated by"),
    Description = c(
    "Association between context-level characteristics and the SMESH cluster assignment",
    "Do the learned clusters reflect demographic, study-design or technical differences between contexts rather than microbiome signal?",
    "One row per context (disease, tumour type, or stage stratum) - not per sample",
    "SMESH cluster label, in the display order used in the main figures",
    "Kruskal-Wallis H (the tie correction is a label-independent constant and is omitted; it cannot change a permutation p-value)",
    "Pearson chi-square; row and column margins are fixed under permutation",
    "Cluster labels permuted, preserving cluster sizes. Enumerated exhaustively where feasible, otherwise 20,000 random permutations; the 'Null' column states which was used for each test.",
    "Benjamini-Hochberg, applied within each application",
    "Contexts are few (12, 15 and 26), so these tests have limited power; a null result does not exclude an association.",
    "Sample-level fields are aggregated on distinct samples, so a control shared between two contexts is counted once. Contexts with a missing value are dropped from that test and 'N_contexts' reports how many remained.",
    "Check_cluster_covariates.R"))

  sheets <- list(
    "README"            = readme,
    "Association tests" = res_all,
    "Contexts pan-disease" = ctx_tab("pan-disease (GMrepo)"),
    "Contexts pan-tumor"   = ctx_tab("pan-tumor (HGMT)"),
    "Contexts CRC"         = ctx_tab("colorectal neoplasia"))

  out_xlsx <- file.path(OUT_DIR, "cluster_covariates.xlsx")
  openxlsx::write.xlsx(sheets, file = out_xlsx,
                       asTable = FALSE, overwrite = TRUE,
                       headerStyle = openxlsx::createStyle(textDecoration = "bold"),
                       colWidths = "auto")

  cat("\n", strrep("=", 96), "\n", sep = "")
  cat("  workbook written: ", out_xlsx, "\n", sep = "")
  for (s in names(sheets)) cat(sprintf("    %-22s %3d rows\n", s, nrow(sheets[[s]])))
