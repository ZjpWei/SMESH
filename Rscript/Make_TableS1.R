# =============================================================================
#   Make_TableS1.R  --  regenerate Supplementary Table S1 (pan-disease studies)
# =============================================================================
#  One row per study-disease case-control comparison actually used in the
#  pan-disease analysis (32 rows over 26 projects, 12 diseases).
#
#  Columns follow the existing Table S1 layout:
#    Disease type | Project ID | Assay type | Adjusted variables
#                 | No. cases | No. controls | No. taxa
#
#  "No. taxa" is the number of species with an estimable PALM association in
#  that study, i.e. species passing the within-study 10% prevalence filter -
#  not the number observed before filtering.
#
#  Writes TableS1_pan_disease.tex (a table environment, ready to \input) and
#  the same content as .csv.
# =============================================================================

  rm(list = ls())
  suppressMessages(library(dplyr))

  proc_dir <- "./GMrepo_analysis/Data"  # processed objects from 0_preprocessing.R
  out_tex  <- "./Figure/TableS1_pan_disease.tex"
  out_csv  <- "./Figure/TableS1_pan_disease.csv"
  dir.create("./Figure", showWarnings = FALSE, recursive = TRUE)

  ci   <- readRDS(file.path(proc_dir, "covariate.interest.rds"))
  ca   <- readRDS(file.path(proc_dir, "covariate.adjust.rds"))
  ss   <- readRDS(file.path(proc_dir, "summary_stat.rds"))
  psum <- openxlsx::read.xlsx(file.path(proc_dir, "GMrepo_project_summary.xlsx"))

  mesh_to_name <- setNames(psum$disease_name, psum$disease_mesh)

  ## Display labels, identical to the ones drawn in the figures.
  ## SOURCE OF TRUTH: GMrepo_analysis/Figure_GMrepo.R (`rename_map`).
  ## Keep the two in step - if a label changes there, change it here.
  rename_map <- c(
    "Arthritis, Rheumatoid"             = "Rheumatoid arthritis",
    "Colitis, Ulcerative"               = "Ulcerative colitis",
    "Crohn Disease"                     = "Crohn disease",
    "Diabetes Mellitus, Type 1"         = "Type 1 diabetes",
    "Diabetes Mellitus, Type 2"         = "Type 2 diabetes",
    "Fatigue Syndrome, Chronic"         = "Chronic fatigue syndrome",
    "Fibromyalgia"                      = "Fibromyalgia",
    "Irritable Bowel Syndrome"          = "IBS",
    "Kidney Failure, Chronic"           = "Chronic kidney failure",
    "Non-alcoholic Fatty Liver Disease" = "Liver NAFLD",
    "Obesity"                           = "Obesity",
    "Psoriasis"                         = "Psoriasis")
  rename_disease <- function(x) dplyr::recode(x, !!!rename_map)

  ## Covariate columns carry their raw metadata names; print them the way the
  ## table already does.
  pretty_cov <- c(host_age = "Age", sex = "Sex", BMI = "BMI", country = "Country")

  tab <- bind_rows(lapply(names(ci), function(k) {
    v    <- ci[[k]]
    adj  <- ca[[k]]
    adj_vars <- if (is.null(adj) || ncol(adj) == 0) "--" else
      paste(sort(unname(pretty_cov[colnames(adj)])), collapse = ", ")
    tibble(
      `Disease type`       = unname(rename_disease(mesh_to_name[sub("^.+_", "", k)])),
      `Project ID`         = sub("_[^_]+$", "", k),
      `Assay type`         = "Metagenomics",
      `Adjusted variables` = adj_vars,
      `No. cases`          = sum(v[, 1] == 1),
      `No. controls`       = sum(v[, 1] == 0),
      `No. taxa`           = sum(!is.na(ss[[k]]$est[, 1]))
    )
  })) %>%
    arrange(`Disease type`, `Project ID`)

  print(as.data.frame(tab), row.names = FALSE)
  cat("\nrows:", nrow(tab),
      "| projects:", length(unique(tab$`Project ID`)),
      "| diseases:", length(unique(tab$`Disease type`)),
      "| distinct samples:", length(unique(unlist(lapply(ci, rownames)))), "\n")

  write.csv(tab, out_csv, row.names = FALSE)

  ## ---- LaTeX ---------------------------------------------------------------
  esc <- function(x) gsub("_", "\\\\_", as.character(x))
  body <- apply(tab, 1, function(r)
    paste0(paste(esc(r), collapse = " & "), " \\\\"))

  tex <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption*{\\textbf{Table S1}: List of studies for the pan-disease application.}",
    "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{llllrrr}",
    "\\toprule",
    "Disease type & Project ID & Assay type & Adjusted variables & No. cases & No. controls & No. taxa \\\\",
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    paste0("\\caption*{\\footnotesize The analysis comprises ", nrow(tab),
           " case--control comparisons from ", length(unique(tab$`Project ID`)),
           " projects across ", length(unique(tab$`Disease type`)), " diseases. ",
           "Six projects contribute to both ulcerative colitis and Crohn disease and share ",
           "their control samples between the two comparisons, so the column totals (",
           sum(tab$`No. cases`) + sum(tab$`No. controls`),
           ") exceed the number of distinct samples (",
           length(unique(unlist(lapply(ci, rownames)))), "). ",
           "Adjusted variables were selected per study by a within-study imbalance test; ",
           "``--'' denotes no adjustment. No. taxa is the number of species with an ",
           "estimable association in that study.}"),
    "\\end{table}")

  writeLines(tex, out_tex)
  cat("\nwritten:", out_tex, "and", out_csv, "\n")
