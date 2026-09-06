# =============================================================================
#   Make_TableS2.R  --  regenerate Supplementary Table S2 (pan-tumor studies)
# =============================================================================
#  One row per study-tumor case-control comparison actually used in the
#  pan-tumor analysis (36 rows over 28 projects, 15 tumor types).
#
#  Read directly from HGMT_analysis/Data/HGMT_16S.Rdata, the object the
#  analysis was run on.  HGMT_info.xlsx is NOT used: it still lists 19 tumor
#  types (including four no longer modelled) and its N.taxa column matches
#  neither the raw nor the filtered genus count of any current study.
#
#  "No. taxa" is the number of genera in the analysis feature set reaching 10%
#  prevalence in that study, i.e. the genera estimable there.
#
#  NOTE the assay column reports HGMT's per-sample `Assay type`, which is WGS
#  for PRJEB38625 and PRJEB6070 even though those batches are named *_16S.
# =============================================================================

  rm(list = ls())
  suppressMessages({ library(dplyr); library(stringr) })

  out_tex <- "./Figure/TableS2_pan_tumor.tex"
  out_csv <- "./Figure/TableS2_pan_tumor.csv"
  dir.create("./Figure", showWarnings = FALSE, recursive = TRUE)

  e <- new.env(); load("HGMT_analysis/Data/HGMT_16S.Rdata", envir = e)
  otu <- e$otu_final; mt <- e$meta_final; ca <- e$covariate.adjust

  ctx_of <- function(k) str_remove(k, "^PRJ[A-Z0-9-]+(?:_(?:16S|WGS|PE|SE))?_")

  ## Display labels, identical to the ones drawn in the figures.
  ## SOURCE OF TRUTH: HGMT_analysis/Figure_HGMT.R (`rename_map`).
  ## Types absent from the map keep their pipeline name (Melanoma, Glioma,
  ## Meningioma), exactly as dplyr::recode leaves them in the figure.
  rename_map <- c(
    "Esophageal Squamous Cell Carcinoma" = "Esophageal ESCC",
    "Acute Lymphoblastic Leukemia"       = "Hematologic ALL",
    "Non-Small-Cell Lung Cancer"         = "Lung NSCLC",
    "Hepatocellular Carcinoma"           = "Liver HCC",
    "Pancreatic Ductal Adenocarcinoma"   = "Pancreatic PDAC",
    "Colorectal Cancer"                  = "Colorectal cancer",
    "Colorectal Polyps"                  = "Colorectal polyps",
    "Breast Cancer"                      = "Breast cancer",
    "Brain Metastasis"                   = "Brain metastases",
    "Thyroid Cancer"                     = "Thyroid cancer",
    "Endometrial Cancer"                 = "Endometrial cancer",
    "Gastric Cancer"                     = "Gastric cancer")

  ## PRJNA626591 is carried as "Brain Metastasis" throughout (the tumour-type
  ## recode in Analysis/0_preprocessing.R, Figure5.R, and here), displayed as
  ## "Brain metastases".  HGMT's own
  ## phenotype for this project is "Brain Neoplasms" (MeSH D001932); the name
  ## used here is a project decision, kept for consistency with the figures.
  fix_label <- function(x) dplyr::recode(x, !!!rename_map)

  ## Feature set: genus present in >=10% of samples in >= 6 tumour types.
  lab <- vapply(names(otu), ctx_of, character(1))
  dl  <- lapply(unique(lab), function(l)
    unique(unlist(lapply(otu[lab == l], function(d) colnames(d)[colMeans(d != 0) >= 0.1]))))
  feature_ID <- names(which(table(unlist(dl)) >= 6))

  pretty_cov <- c(age = "Age", sex = "Sex", BMI = "BMI", country = "Country")

  tab <- bind_rows(lapply(names(otu), function(k) {
    d  <- otu[[k]]; md <- mt[[k]]
    adj <- ca[[k]]
    adj_vars <- if (is.null(adj) || ncol(adj) == 0) "--" else
      paste(sort(unname(pretty_cov[colnames(adj)])), collapse = ", ")
    assay <- unique(as.character(md$`Assay type`))
    assay <- ifelse(assay == "16", "16S", assay)
    tibble(
      `Tumor type`         = unname(fix_label(ctx_of(k))),
      `Project ID`         = gsub("-", ", ", as.character(md$`Project ID`[1])),
      `Assay type`         = paste(sort(unique(assay)), collapse = ", "),
      `Adjusted variables` = adj_vars,
      `No. cases`          = sum(md$`Phenotype name` != "Health"),
      `No. controls`       = sum(md$`Phenotype name` == "Health"),
      `No. taxa`           = length(intersect(feature_ID,
                                              colnames(d)[colMeans(d != 0) >= 0.1]))
    )
  })) %>% arrange(`Tumor type`, `Project ID`)

  print(as.data.frame(tab), row.names = FALSE)

  ids <- unlist(lapply(mt, function(z) as.character(z$`Run ID`)))
  cat("\nrows:", nrow(tab),
      "| projects:", length(unique(tab$`Project ID`)),
      "| tumour types:", length(unique(tab$`Tumor type`)),
      "| distinct samples:", length(unique(ids)),
      "| column sum:", sum(tab$`No. cases`) + sum(tab$`No. controls`), "\n")
  cat("assay mix:\n"); print(table(tab$`Assay type`))

  write.csv(tab, out_csv, row.names = FALSE)

  esc  <- function(x) gsub("_", "\\\\_", as.character(x))
  body <- apply(tab, 1, function(r) paste0(paste(esc(r), collapse = " & "), " \\\\"))

  tex <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption*{\\textbf{Table S2}: List of studies for the pan-tumor application.}",
    "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{llllrrr}",
    "\\toprule",
    "Tumor type & Project ID & Assay type & Adjusted variables & No. cases & No. controls & No. taxa \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    paste0("\\caption*{\\footnotesize The analysis comprises ", nrow(tab),
           " case--control comparisons from ", length(unique(tab$`Project ID`)),
           " projects across ", length(unique(tab$`Tumor type`)), " tumor types. ",
           "Seven projects contribute to more than one tumor type and share their control ",
           "samples between comparisons, so the column totals (",
           sum(tab$`No. cases`) + sum(tab$`No. controls`),
           ") exceed the number of distinct samples (", length(unique(ids)), "). ",
           "Adjusted variables were selected per study by a within-study imbalance test; ",
           "``--'' denotes no adjustment. No. taxa is the number of genera estimable in ",
           "that study.}"),
    "\\end{table}")

  writeLines(tex, out_tex)
  cat("\nwritten:", out_tex, "and", out_csv, "\n")
