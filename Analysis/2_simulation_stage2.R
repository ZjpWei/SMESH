# ============================================================================
#   2_simulation_stage2.R  --  stage 2 simulations, all four input summaries
# ============================================================================
#  One replicate of the stage-2 simulation: build semi-synthetic studies from
#  the colorectal template, derive context-level association summaries with the
#  requested method, fit SMESH, and score cluster recovery and signature
#  detection against the truth.
#
#  This merges Simulation_PALM.r, Simulation_ANCOMBC2.r, Simulation_LinDA.r and
#  Simulation_MaAsLin3.r.  Those four differed only in how the summary
#  statistics are produced - the data simulation and the scoring were identical
#  - so the shared parts are written once and the method-specific block is
#  chosen by the sixth argument.  No method parameter has been changed.
#
#  The shared section runs before the branch, so a given (s, Ka.per, pos.pt, u,
#  scenario) yields the SAME simulated data for every method, exactly as when
#  the four scripts were run separately.
#
#  Usage, from the project root:
#      Rscript Analysis/2_simulation_stage2.R <s> <Ka.per> <pos.pt> <u> <scenario> <method>
#      Rscript Analysis/2_simulation_stage2.R 1 0.4 0.5 0 3 PALM
#
#    s         replicate number, 1..100
#    Ka.per    signature percentage {0.1, 0.2, 0.4}
#    pos.pt    signature effect direction {0.5, 0.75, 1}
#    u         unevenness
#    scenario  scenario id, 1..5
#    method    PALM | ANCOMBC2 | LinDA | MaAsLin3
# ============================================================================

  rm(list = ls())

  library("coin")
  library("purrr")
  library("MIDASim")
  library("PALM")
  library("phyloseq")
  library("glmnet")
  library("glmtlp")
  library("dplyr")
  library("tidyr")
  library("tibble")
  library("ggplot2")
  library("mclust")
  library("sClust")
  library("sparcl")
  library("tidyverse")
  library("cluster")
  library("ANCOMBC")
  library("maaslin3")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6)
    stop("need <s> <Ka.per> <pos.pt> <u> <scenario> <method>, e.g.\n",
         "  Rscript Analysis/2_simulation_stage2.R 1 0.4 0.5 0 3 PALM")
  print(args)

  s        <- as.numeric(args[1])
  Ka.per   <- as.numeric(args[2])
  pos.pt   <- as.numeric(args[3])
  u        <- as.numeric(args[4])
  scenario <- as.numeric(args[5])
  method   <- args[6]

  if (!method %in% names(FILE_PREFIX <- c(PALM = "palm", ANCOMBC2 = "ancombc2",
                                          LinDA = "linda", MaAsLin3 = "maaslin3")))
    stop("unknown method '", method, "'; choose PALM, ANCOMBC2, LinDA or MaAsLin3")

  method_label <- c(PALM = "SMESH-PALM", ANCOMBC2 = "SMESH-ANCOMBC2",
                    LinDA = "SMESH-LinDA", MaAsLin3 = "SMESH-MaAsLin3")[[method]]

  ## ---- paths, relative to the repository root --------------------------------
  if (!dir.exists("utility"))
    stop("run this from the project root, e.g. Rscript Analysis/2_simulation_stage2.R ...")

  out_dir <- "./Simulation/Sim_CRC_stage2"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ## MaAsLin3 writes its per-study working files here.
  dir.create("./tmp_file", showWarnings = FALSE, recursive = TRUE)

  data.loc <- file.path(out_dir, paste0(FILE_PREFIX[[method]], "_res_Ka", Ka.per,
                                        "_pos", pos.pt, "_u", u,
                                        "_scenario", scenario, "_s", s, ".Rdata"))

  ## Template: the colorectal strata object, whichever copy this checkout has.
  template <- c("./Simulation/Data/CRC.Rdata", "./CRC_analysis/Data/CRC_strata.Rdata")
  template <- template[file.exists(template)]
  if (!length(template)) stop("no CRC template found under ./Data/")
  load(template[1])
  message("template: ", template[1], " | method: ", method)

#================ Melody-test one data simulation ================#
  
  # Packages ----
  
  # Simulation settings ----
  ## Sample sizes are NOT set by the simulation: each simulated study inherits
  ## the sample size of the real stratum it is built from (n <- nrow(...) below),
  ## along with that stratum's library-size distribution and taxon structure.
  ## The two arms are floor(n/2) each, so the design is balanced even where the
  ## real stratum is not.
  
  ## random seed
  set.seed(s + 2026)
  
  ## effect size and variation
  delta.gb <- c(2, 2.5)
  delta.12 <- c(3, 4)
  delta.1 = delta.2 = delta.3 <- c(5, 7)
  
  ## study proportion and homogeneity proportion
  if(scenario == 2){
    global.hom <- 1/2
    G1.prop = G2.prop <- 1/4
  }else if(scenario == 3){
    global.hom <- 1/2
    G12.prop = G1.prop = G2.prop = G3.prop <- 1/8
  }else if(scenario == 4){
    global.hom <- 1/4
    G1.prop = G2.prop <- 3/8
  }else if(scenario == 5){
    global.hom <- 1/4
    G12.prop = G1.prop = G2.prop = G3.prop <- 3/16
  }
  
  ## Case/control sequence depth unevenness {0, 1}
  mu <- 0
  
  ## Signature sparsity
  abd.pt <- 0.5
  
  ## directory for saving results
  
  # Simulate data ----
  
  ## Consider disease level feature cutoff
  feature_ID <- names(which(table(unlist(sapply(otu_data_sp, function(d){
    colnames(d)[colMeans(d!=0) >= 0.1]
  }))) >= 23))
  taxa_num <- length(feature_ID)
  
  rd_study <- sort(sample(length(otu_data_sp), size = 18))
  otu_data_filter <- sapply(otu_data_sp[rd_study], function(d){
    tmp_d <- d[,intersect(colnames(d), feature_ID)]
    tmp_d <- tmp_d[rowSums(tmp_d) > 0,]
    tmp_d <- tmp_d[,colSums(tmp_d) > 0]
    return(tmp_d)
  })
  
  re.lst <- names(otu_data_filter)
  L <- length(re.lst)
  uni.features <- unique(unlist(lapply(otu_data_filter, function(d){colnames(d)})))
  
  ## Signal add
  if(scenario == 1){
    ##  ---- Signal for global shared
    global.id <- sample(uni.features, length(uni.features) * Ka.per)
    dir.global <- rbernoulli(n = length(global.id), p = pos.pt)
    beta.global <- runif(length(global.id), min = delta.gb[1], max = delta.gb[2])
    names(beta.global) <- global.id
    
    G.lst <- rep(1, L)
    names(G.lst) <- re.lst
    signal.lst     <- list(global.id)
    signal.lst.hom <- global.id
    signal.lst.het <- NULL
  }else{
    ## ---- Signal setting
    hom.id <- sample(uni.features, length(uni.features) * Ka.per * global.hom)
    G1.id  <- sample(setdiff(uni.features, hom.id), length(uni.features) * Ka.per * G1.prop)
    G2.id  <- sample(setdiff(uni.features, c(hom.id, G1.id)), length(uni.features) * Ka.per * G2.prop)
    
    ## ---- Signal for global shared (hom)
    dir.hom   <- rbernoulli(n = length(hom.id), p = pos.pt)
    beta.hom  <- runif(length(hom.id), min = delta.gb[1], max = delta.gb[2])
    names(beta.hom) <- hom.id
    
    ## ---- Signal for G = 1
    dir.G1  <- rbernoulli(n = length(G1.id), p = pos.pt)
    beta.G1 <- runif(length(G1.id), min = delta.1[1], max = delta.1[2])
    names(beta.G1) <- G1.id
    
    ## ---- Signal for G = 2
    dir.G2  <- rbernoulli(n = length(G2.id), p = pos.pt)
    beta.G2 <- runif(length(G2.id), min = delta.2[1], max = delta.2[2])
    names(beta.G2) <- G2.id
    
    if(scenario == 3 | scenario == 5){
      ## ---- Signal setting for G = 3
      G3.id  <- sample(setdiff(uni.features, c(hom.id, G1.id, G2.id)), length(uni.features) * Ka.per * G3.prop)
      G12.id <- sample(setdiff(uni.features, c(hom.id, G1.id, G2.id, G3.id)), length(uni.features) * Ka.per * G12.prop)
      
      ## ---- Signal for G = 1&2
      dir.G12  <- rbernoulli(n = length(G12.id), p = pos.pt)
      beta.G12 <- runif(length(G12.id), min = delta.12[1], max = delta.12[2])
      names(beta.G12) <- G12.id
      
      ## ---- Signal for G = 3
      dir.G3  <- rbernoulli(n = length(G3.id), p = pos.pt)
      beta.G3 <- runif(length(G3.id), min = delta.3[1], max = delta.3[2])
      names(beta.G3) <- G3.id
      
      G.lst <- c(rep(1, round(L/2)), rep(2, round(L/4)), rep(3, L - round(L/2) - round(L/4)))
      names(G.lst)   <- re.lst
      signal.lst     <- list(c(hom.id, G12.id, G1.id), c(hom.id, G12.id, G2.id), c(hom.id, G3.id))
      signal.lst.hom <- hom.id
      signal.lst.het <- list(c(G12.id, G1.id), c(G12.id, G2.id), G3.id)
    }else{
      G.lst <- c(rep(1, round(L/2)), rep(2, L - round(L/2)))
      names(G.lst)   <- re.lst
      signal.lst     <- list(c(hom.id, G1.id), c(hom.id, G2.id))
      signal.lst.hom <- hom.id
      signal.lst.het <- list(G1.id, G2.id)
    }
  }
  G <- max(G.lst)
  
  ## Original data
  data.rel <- list()
  for(d in re.lst){
    n <- nrow(otu_data_filter[[d]])
    
    ## Add signal
    beta.neg = beta.pos <- rep(1, length(uni.features))
    names(beta.pos) = names(beta.neg) <- uni.features
    if(scenario == 1){
      beta.neg[names(beta.global[ dir.global])] <- beta.global[ dir.global]
      beta.pos[names(beta.global[!dir.global])] <- beta.global[!dir.global]
    }else if(scenario == 3 || scenario == 5){
      if(G.lst[d] == 1){
        beta.pos[names(c(beta.hom[ dir.hom], beta.G12[ dir.G12], beta.G1[ dir.G1]))] <- c(beta.hom[ dir.hom], beta.G12[ dir.G12], beta.G1[ dir.G1])
        beta.neg[names(c(beta.hom[!dir.hom], beta.G12[!dir.G12], beta.G1[!dir.G1]))] <- c(beta.hom[!dir.hom], beta.G12[!dir.G12], beta.G1[!dir.G1])
      }else if(G.lst[d] == 2){
        beta.pos[names(c(beta.hom[ dir.hom], beta.G12[ dir.G12], beta.G2[ dir.G2]))] <- c(beta.hom[ dir.hom], beta.G12[ dir.G12], beta.G2[ dir.G2])
        beta.neg[names(c(beta.hom[!dir.hom], beta.G12[!dir.G12], beta.G2[!dir.G2]))] <- c(beta.hom[!dir.hom], beta.G12[!dir.G12], beta.G2[!dir.G2])
      }else{
        beta.pos[names(c(beta.hom[ dir.hom], beta.G3[ dir.G3]))] <- c(beta.hom[ dir.hom], beta.G3[ dir.G3])
        beta.neg[names(c(beta.hom[!dir.hom], beta.G3[!dir.G3]))] <- c(beta.hom[!dir.hom], beta.G3[!dir.G3])
      }
    }else{
      if(G.lst[d] == 1){
        beta.pos[names(c(beta.hom[ dir.hom], beta.G1[ dir.G1]))] <- c(beta.hom[ dir.hom], beta.G1[ dir.G1])
        beta.neg[names(c(beta.hom[!dir.hom], beta.G1[!dir.G1]))] <- c(beta.hom[!dir.hom], beta.G1[!dir.G1])
      }else{
        beta.pos[names(c(beta.hom[ dir.hom], beta.G2[ dir.G2]))] <- c(beta.hom[ dir.hom], beta.G2[ dir.G2])
        beta.neg[names(c(beta.hom[!dir.hom], beta.G2[!dir.G2]))] <- c(beta.hom[!dir.hom], beta.G2[!dir.G2])
      }
    }
    
    ## MIDSim setup
    Y <- otu_data_filter[[d]]
    count.ibd.setup = MIDASim.setup(Y, mode = 'nonparametric', n.break.ties = 100)
    
    ## control data
    lib.size <- sample(rowSums(Y), size = floor(n/2), replace = TRUE)
    mean.avg.prop <- colMeans(Y / rowSums(Y))
    
    ## add signal  version 2 (positive on control)
    mean.avg.prop <- mean.avg.prop * beta.neg[names(mean.avg.prop)] 
    mean.avg.prop <- mean.avg.prop / sum(mean.avg.prop)
    
    ## formula from GitHub
    obs.sample.1.ct <- rowSums(Y > 0)
    xvar <- log10(rowSums(Y))
    scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
    sample.1.ct = 10^(predict(scamfit.non0, newdata = data.frame(xvar = log10(lib.size) )) )
    n.taxa = ncol(Y)
    input.sample.prop = sample.1.ct / n.taxa
    input.taxa.prop <- count.ibd.setup$taxa.1.prop * (sum(input.sample.prop) * ncol(Y) / floor(n/2)) / sum(count.ibd.setup$taxa.1.prop)
    
    ## Simulate data
    count.ibd.modified <- MIDASim.modify(count.ibd.setup,
                                         lib.size = lib.size,
                                         mean.rel.abund = mean.avg.prop,
                                         sample.1.prop = input.sample.prop,
                                         taxa.1.prop = input.taxa.prop)
    
    simulated.data <- MIDASim(count.ibd.modified)
    control.data <- simulated.data$sim_count
    
    ## case data
    lib.size <- sample(rowSums(Y), size = floor(n/2), replace = TRUE) 
    mean.avg.prop <- colMeans(Y/rowSums(Y))
    
    ## add signal  version 2 (positive on case)
    mean.avg.prop <- mean.avg.prop * beta.pos[names(mean.avg.prop)] 
    mean.avg.prop <- mean.avg.prop / sum(mean.avg.prop)
    
    ## formula from GitHub
    obs.sample.1.ct <- rowSums(Y > 0)
    xvar <- log10(rowSums(Y))
    scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
    sample.1.ct = 10^(predict(scamfit.non0, newdata = data.frame(xvar = log10(lib.size) )) )
    n.taxa = ncol(Y)
    input.sample.prop = sample.1.ct / n.taxa
    input.taxa.prop <- count.ibd.setup$taxa.1.prop * (sum(input.sample.prop) * ncol(Y) / floor(n/2)) / sum(count.ibd.setup$taxa.1.prop)
    
    ## Simulate data
    count.ibd.modified <- MIDASim.modify(count.ibd.setup,
                                         lib.size = lib.size,
                                         mean.rel.abund = mean.avg.prop,
                                         sample.1.prop = input.sample.prop,
                                         taxa.1.prop = input.taxa.prop)
    
    simulated.data <- MIDASim(count.ibd.modified)
    case.data <- simulated.data$sim_count
    
    ## summarize data
    data.tmp <- rbind(control.data, case.data)
    rownames(data.tmp) <- paste0("CS_", d, "_", 1:(2*floor(n/2)))
    data.rel[[d]] <- list(Y = data.tmp, X = rep(c(0,1), each = floor(n/2)))
  }

  # ============================================================================================== #

  ## ---- method-specific: summary statistics, then the SMESH fit ---------------
  switch(method,
    "PALM" = {
  ## output information
  Precision = Recall = F1 = Precision_hom = Recall_hom = F1_hom = Precision_het = Recall_het = F1_het = ARI = G_hat <- NULL
  rel.abd <- list()
  covariate.interest <- list()
  for(d in names(data.rel)){
    rel.abd[[d]] <- data.rel[[d]]$Y
    covariate.interest[[d]] <- matrix(data.rel[[d]]$X, ncol = 1,
                                      dimnames = list(rownames(data.rel[[d]]$Y), "disease"))
  }
  
  ## PALM summary statistics (median)
  null.obj <- PALM::palm.null.model(rel.abd = rel.abd, prev.filter = 0)
  
  summary.stat.obj <- PALM::palm.get.summary(null.obj = null.obj,
                                             covariate.interest = covariate.interest,
                                             correct = "tune")
  
  ## SMESH model
  source("./utility/SMESH.R")
  
  SMESH_fit <- smesh.meta.summary(summary.stats = summary.stat.obj, 
                                 nperm = 5,
                                 doc = "./",
                                 verbose = TRUE)
  
  ## ---- summarize the results
    },
    "ANCOMBC2" = {
  ## output information
  Precision = Recall = F1 = Precision_hom = Recall_hom = F1_hom = Precision_het = Recall_het = F1_het = ARI = G_hat <- NULL
  rel.abd <- list()
  covariate.interest <- list()
  outcome = feature.table = studys <- NULL
  for(d in names(data.rel)){
    rel.abd[[d]] <- data.rel[[d]]$Y
    tmp.count <- matrix(0, nrow = nrow(data.rel[[d]]$Y), ncol = length(uni.features),
                        dimnames = list(rownames(data.rel[[d]]$Y), uni.features))
    tmp.count[,colnames(data.rel[[d]]$Y)] <- data.rel[[d]]$Y
    feature.table <- rbind(feature.table, tmp.count)
    outcome <- c(outcome, data.rel[[d]]$X)
    studys <- c(studys, rep(d, length(data.rel[[d]]$X)))
    covariate.interest[[d]] <- matrix(data.rel[[d]]$X, ncol = 1, dimnames = list(rownames(data.rel[[d]]$Y), "disease"))
  }
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(studys))
  rownames(meta.data) <- colnames(feature.table)
  feature.ID <- rownames(feature.table)

  ## ANCOM-BC2
  source("./utility/ancombc.R")
  
  ancombc.sums <- list()
  for(d in names(rel.abd)){
    AA.est <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    AA.std <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    
    data.tmp <- feature.table[,meta.data$study == d]
    data.tmp <- data.tmp[rowMeans(data.tmp!=0)>=0.1,]
    data.tmp <- data.tmp[apply(data.tmp, 1, var) > 0, ]
    ANCOMBC2.model <- ancombc.fun(feature.table = data.tmp,
                                  meta = meta.data %>% dplyr::filter(study == d),
                                  formula = "labels",
                                  adjust.method = "fdr",
                                  group = NULL,
                                  subject = NULL,
                                  method = "ancombc2")
    
    AA.est[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],1] <-
      ANCOMBC2.model$res$lfc_labels1[ANCOMBC2.model$res$passed_ss_labels1]
    AA.std[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],1] <-
      ANCOMBC2.model$res$se_labels1[ANCOMBC2.model$res$passed_ss_labels1]
    
    ancombc.sums[[d]] <- list(est = AA.est, stderr = AA.std, n = ncol(data.tmp))
  }
  
  source("./utility/SMESH.R")
  
  SMESH_fit <-  smesh.meta.summary(summary.stats = ancombc.sums, 
                                    nperm = 5,
                                    doc = "./",
                                    verbose = FALSE)

  ## ---- summarize the results
    },
    "LinDA" = {
  ## output information
  Precision = Recall = F1 = Precision_hom = Recall_hom = F1_hom = Precision_het = Recall_het = F1_het = ARI = G_hat <- NULL
  rel.abd <- list()
  covariate.interest <- list()
  outcome = feature.table = studys <- NULL
  for(d in names(data.rel)){
    rel.abd[[d]] <- data.rel[[d]]$Y
    tmp.count <- matrix(0, nrow = nrow(data.rel[[d]]$Y), ncol = length(uni.features),
                        dimnames = list(rownames(data.rel[[d]]$Y), uni.features))
    tmp.count[,colnames(data.rel[[d]]$Y)] <- data.rel[[d]]$Y
    feature.table <- rbind(feature.table, tmp.count)
    outcome <- c(outcome, data.rel[[d]]$X)
    studys <- c(studys, rep(d, length(data.rel[[d]]$X)))
    covariate.interest[[d]] <- matrix(data.rel[[d]]$X, ncol = 1, dimnames = list(rownames(data.rel[[d]]$Y), "disease"))
  }
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(studys))
  rownames(meta.data) <- colnames(feature.table)
  feature.ID <- rownames(feature.table)

  ## LinDA
  Linda.sums <- list()
  for(d in names(rel.abd)){
    AA.est <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    AA.std <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    data.tmp <- feature.table[,meta.data$study == d]
    data.tmp <- data.tmp[rowMeans(data.tmp!=0)>=0.1,]
    data.tmp <- data.tmp[apply(data.tmp, 1, var) > 0, ]
    
    Linda.model <- MicrobiomeStat::linda(feature.dat = data.tmp,
                                         meta.dat = meta.data %>% dplyr::filter(study == d),
                                         formula = '~ labels',
                                         feature.dat.type = "count",
                                         prev.filter = 0,
                                         adaptive = TRUE,
                                         alpha = 0.05)
    
    AA.est[rownames(Linda.model$output$labels1),1] <- Linda.model$output$labels1$log2FoldChange
    AA.std[rownames(Linda.model$output$labels1),1] <- Linda.model$output$labels1$lfcSE
    Linda.sums[[d]] <- list(est = AA.est, stderr = AA.std, n = ncol(data.tmp))
  }
  
  source("./utility/SMESH.R")
  
  SMESH_fit <-  smesh.meta.summary(summary.stats = Linda.sums,  
                                     nperm = 5,
                                     doc = "./",
                                     verbose = FALSE)
  
  ## ---- summarize the results
    },
    "MaAsLin3" = {
  ## output information
  Precision = Recall = F1 = Precision_hom = Recall_hom = F1_hom = Precision_het = Recall_het = F1_het = ARI = G_hat <- NULL
  rel.abd <- list()
  covariate.interest <- list()
  outcome = feature.table = studys <- NULL
  for(d in names(data.rel)){
    rel.abd[[d]] <- data.rel[[d]]$Y
    tmp.count <- matrix(0, nrow = nrow(data.rel[[d]]$Y), ncol = length(uni.features),
                        dimnames = list(rownames(data.rel[[d]]$Y), uni.features))
    tmp.count[,colnames(data.rel[[d]]$Y)] <- data.rel[[d]]$Y
    feature.table <- rbind(feature.table, tmp.count)
    outcome <- c(outcome, data.rel[[d]]$X)
    studys <- c(studys, rep(d, length(data.rel[[d]]$X)))
    covariate.interest[[d]] <- matrix(data.rel[[d]]$X, ncol = 1, dimnames = list(rownames(data.rel[[d]]$Y), "disease"))
  }
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(studys))
  rownames(meta.data) <- colnames(feature.table)
  feature.ID <- rownames(feature.table)

  ## MaAslin3
  MaAsLin3.sums <- list()
  for(d in names(rel.abd)){
    AA.est <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    AA.std <- matrix(NA, nrow = length(feature.ID), ncol = 1, dimnames = list(feature.ID, "disease"))
    
    data.tmp <- feature.table[,meta.data$study == d]
    meta.tmp <- meta.data[meta.data$study == d, , drop = FALSE]
    
    Maaslin3_res <- maaslin3(input_data = data.tmp,
                             input_metadata = data.frame(meta.tmp),
                             output = paste0('./tmp_file/',d), #paste0('./log/',d), #
                             formula = '~ labels',
                             normalization = 'TSS',
                             transform = 'LOG',
                             augment = TRUE,
                             standardize = FALSE,
                             max_significance = 0.1,
                             median_comparison_abundance = TRUE,
                             median_comparison_prevalence = FALSE,
                             max_pngs = 100,
                             save_models = FALSE)
    
    tmp_model <- Maaslin3_res$fit_data_abundance$results
    tmp_model <- tmp_model %>% dplyr::filter(!is.na(coef), !is.na(stderr))
    
    AA.est[tmp_model$feature,1] <- tmp_model$coef
    AA.std[tmp_model$feature,1] <- tmp_model$stderr
    MaAsLin3.sums[[d]] <- list(est = AA.est, stderr = AA.std, n = ncol(data.tmp))
  }
  
  source("./utility/SMESH.R")
  
  SMESH_fit <-  smesh.meta.summary(summary.stats = MaAsLin3.sums,
                                       nperm = 5,
                                       doc = "./",
                                       verbose = FALSE)
  
  ## ---- summarize the results
    },
    stop("unreachable")
  )

  est_clust <- apply(round(SMESH_fit$disease$W), 1, function(d){which(d == 1)})
  ARI <- c(ARI, adjustedRandIndex(est_clust, G.lst))
  G_hat <- c(G_hat, sum(colSums(round(SMESH_fit$disease$W)) > 0))
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- G.lst[d]
    g2 <- which.max(SMESH_fit$disease$W[d,])
    truth.signal <- signal.lst[[G.lst[d]]]
    selected.signal <- rownames(SMESH_fit$disease$mu)[SMESH_fit$disease$mu[,g2]!=0]
    pre <- c(pre, length(intersect(truth.signal, selected.signal))/length(selected.signal))
    rec <- c(rec, length(intersect(truth.signal, selected.signal))/length(truth.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre, na.rm = TRUE))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))
  
  ## ---- homogeneous part
  detect_hom <- names(which(apply(SMESH_fit$disease$mu, 1, function(d){all(d !=0)})))
  selected.signal <- signal.lst.hom
  detected.signal <- detect_hom
  pre <- length(intersect(detected.signal, selected.signal)) / length(detected.signal)
  rec <- length(intersect(detected.signal, selected.signal)) / length(selected.signal)
  f1 <- 2 / (1/pre + 1/rec)
  Precision_hom <- c(Precision_hom, pre)
  Recall_hom <- c(Recall_hom, rec)
  F1_hom <- c(F1_hom, f1)
  
  ## ---- heterogeneous part
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- G.lst[d]
    g2 <- which.max(SMESH_fit$disease$W[d,])
    selected.signal <- signal.lst.het[[g]]
    detected.signal <- setdiff(rownames(SMESH_fit$disease$mu)[SMESH_fit$disease$mu[,g2]!=0], detect_hom)
    pre <- c(pre, length(intersect(detected.signal, selected.signal)) / length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal)) / length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, mean(pre))
  Recall_het <- c(Recall_het, mean(rec))
  F1_het <- c(F1_het, mean(f1))

  ## ---- output
  result_mat <- data.frame(
    Precision = Precision, Recall = Recall, F1 = F1,
    Precision_hom = Precision_hom, Recall_hom = Recall_hom, F1_hom = F1_hom,
    Precision_het = Precision_het, Recall_het = Recall_het, F1_het = F1_het,
    G_hat = G_hat, ARI = ARI,
    method = method_label
  )

  save(result_mat, file = data.loc)

