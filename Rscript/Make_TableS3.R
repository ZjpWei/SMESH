# =============================================================================
#   Make_TableS3.R  --  regenerate Supplementary Table S3 (colorectal neoplasia)
# =============================================================================
#  One row per study-stage stratum used in the colorectal neoplasia analysis
#  (26 strata from 12 studies).
#
#  STUDY NAMING.  The SegataLab repository ships files named with
#  curatedMetagenomicData-style dataset identifiers (WirbelJ_2018, YuJ_2015,
#  NSHII, ...).  Piccinno et al. (2025) label the same datasets by citation year
#  in their figures, and two of them differ from the dataset identifier:
#
#      repository file        Piccinno et al. (2025)
#      ------------------     ----------------------
#      WirbelJ_2018       ->  Wirbel, J. (2019)
#      YuJ_2015           ->  Yu, J. (2017)
#      NSHII              ->  NHSII  (Cohort 5; Nurses' Health Study II)
#      IIGM_TU            ->  IIGM TUCRC (Cohort 6)
#      ONCOBIOME_IIGM_CZ  ->  Cohort 3
#      ONCOBIOME_IIGM_IT  ->  Cohort 4
#
#  The dataset identifiers are not wrong - WirbelJ_2018 and YuJ_2015 are the
#  conventional identifiers, and the year refers to the data rather than the
#  publication - but the table should match the source publication, so the
#  published labels are used here and the identifier is kept alongside.
#
#  "No. taxa" is the number of the 260 analysed species present in that stratum.
# =============================================================================

  rm(list = ls())
  suppressMessages(library(dplyr))

  out_tex <- "./Figure/TableS3_colorectal.tex"
  out_csv <- "./Figure/TableS3_colorectal.csv"
  dir.create("./Figure", showWarnings = FALSE, recursive = TRUE)

  lo <- function(f, nm) { e <- new.env(); load(f, envir = e); get(nm, envir = e) }
  md  <- lo("CRC_analysis/Data/CRC_strata.Rdata", "meta_data_sp")
  ot  <- lo("CRC_analysis/Data/CRC_strata.Rdata", "otu_data_sp")
  mdl <- lo("CRC_analysis/CRC_loso/Model2_SMESH_s100.Rdata", "SMESH_model")
  keep <- rownames(mdl$disease$W)

  ## The 260 species carried into the analysis.
  feature_ID <- names(which(table(unlist(lapply(ot, function(d)
    colnames(d)[colMeans(d != 0) >= 0.1]))) >= 23))

  ## Two names per study.
  ##
  ##  study_paper : the designation used by Piccinno et al. (2025).  Three of the
  ##                repository identifiers carry the wrong year or a typo:
  ##                  WirbelJ_2018 -> WirbelJ_2019   (Nat Med 2019)
  ##                  YuJ_2015     -> YuJ_2017       (Gut 2017)
  ##                  NSHII        -> NHSII          (Nurses' Health Study II)
  ##                The four newly sequenced cohorts are Cohorts 3-6 there.
  ##
  ##  study_fig   : the short label drawn in Figure 6.  MUST match the
  ##                `rename_study` map in CRC_analysis/Figure_CRC.R - if one
  ##                changes, change the other.
  study_paper <- c(
    "FengQ_2015"        = "Feng, Q. (2015)",
    "GuptaA_2019"       = "Gupta, A. (2019)",
    "VogtmannE_2016"    = "Vogtmann, E. (2016)",
    "WirbelJ_2018"      = "Wirbel, J. (2019)",
    "YachidaS_2019"     = "Yachida, S. (2019)",
    "YangJ_2020"        = "Yang, J. (2020)",
    "YuJ_2015"          = "Yu, J. (2017)",
    "ZellerG_2014"      = "Zeller, G. (2014)",
    "ONCOBIOME_IIGM_CZ" = "Cohort 3",
    "ONCOBIOME_IIGM_IT" = "Cohort 4",
    "NSHII"             = "Cohort 5",
    "IIGM_TU"           = "Cohort 6")

  study_fig <- c(
    "FengQ_2015"        = "FengQ_2015",
    "GuptaA_2019"       = "GuptaA_2019",
    "VogtmannE_2016"    = "VogtmannE_2016",
    "WirbelJ_2018"      = "WirbelJ_2019",
    "YachidaS_2019"     = "YachidaS_2019",
    "YangJ_2020"        = "YangJ_2020",
    "YuJ_2015"          = "YuJ_2017",
    "ZellerG_2014"      = "ZellerG_2014",
    "ONCOBIOME_IIGM_CZ" = "IIGM_CZ",
    "ONCOBIOME_IIGM_IT" = "IIGM_IT",
    "NSHII"             = "NHSII",
    "IIGM_TU"           = "IIGM_TU")

  stage_label <- c(Adenoma = "Adenoma", early = "Early", late = "Late")

  tab <- bind_rows(lapply(keep, function(k) {
    z <- md[[k]]; d <- ot[[k]]
    sid <- sub(":.*$", "", k)
    tibble(
      `CRC stage`     = unname(stage_label[sub("^.*:", "", k)]),
      `Study`         = unname(study_paper[sid]),
      `Label in Fig. 6` = unname(study_fig[sid]),
      `No. cases`     = sum(z$Group == 1),
      `No. controls`  = sum(z$Group == 0),
      `No. taxa`      = length(intersect(feature_ID, colnames(d)[colSums(d) > 0]))
    )
  })) %>%
    mutate(`CRC stage` = factor(`CRC stage`, levels = c("Adenoma", "Early", "Late"))) %>%
    arrange(`CRC stage`, Study)

  if (anyNA(tab$Study) || anyNA(tab$`Label in Fig. 6`))
    stop("unmapped study identifier")

  print(as.data.frame(tab), row.names = FALSE)
  ids <- unlist(lapply(md[keep], function(z) as.character(z$SampleID)))
  cat("\nstrata:", nrow(tab), "| studies:", length(unique(tab$Study)),
      "| samples:", length(unique(ids)),
      "| column sum:", sum(tab$`No. cases`) + sum(tab$`No. controls`), "\n")

  write.csv(tab, out_csv, row.names = FALSE)

  esc  <- function(x) gsub("_", "\\\\_", as.character(x))
  body <- apply(tab, 1, function(r) paste0(paste(esc(r), collapse = " & "), " \\\\"))

  tex <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{0.82}",
    "\\caption*{{\\bf{Table S3}}: List of studies for the colorectal neoplasia application.}",
        "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{lllrrr}",
    "\\toprule",
        "CRC stage & Study & Fig.~6 label & No. cases & No. controls & No. taxa \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    paste0("\\caption*{\\footnotesize Study designations follow Piccinno et al. (2025); ",
           "``Label in Fig.~6'' is the short name drawn in the main figure. ",
           "The analysis comprises ", nrow(tab), " study--stage strata from ",
           length(unique(tab$Study)), " studies and ", length(unique(ids)),
           " samples. Controls were partitioned across the stage strata of a study, so no ",
           "sample contributes to more than one stratum. No. taxa is the number of the 260 ",
           "analysed species present in that stratum.}"),
    "\\end{table}")

  writeLines(tex, out_tex)
  cat("\nwritten:", out_tex, "and", out_csv, "\n")
