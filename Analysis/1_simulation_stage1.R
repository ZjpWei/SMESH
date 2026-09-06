# ============================================================================
#   1_simulation_stage1.R  --  stage 1 simulations, benchmark methods
# ============================================================================
#  One replicate of the stage-1 simulation: build semi-synthetic studies from
#  the colorectal template, then compare SMESH against True-cluster FE,
#  SKM + FE, SHC + FE and Melody on the same data.
#
#  Stage 1 and stage 2 are split because running every method on one machine
#  takes a long time; stage 2 varies the input summary statistics instead
#  (see Analysis/2_simulation_stage2.R).
#
#  Usage, from the project root:
#      Rscript Analysis/1_simulation_stage1.R <s> <Ka.per> <pos.pt> <u> <scenario>
#      Rscript Analysis/1_simulation_stage1.R 1 0.4 0.5 0 3
#
#    s         replicate number, 1..100
#    Ka.per    signature percentage {0.1, 0.2, 0.4}
#    pos.pt    signature effect direction {0.5, 0.75, 1}
#    u         unevenness
#    scenario  scenario id, 1..5
#
#  No method parameter has been changed from the original Simulation_1.r; only
#  the input and output paths are resolved relative to the repository root.
# ============================================================================

  rm(list = ls())

  library("coin")
  library("purrr")
  library("MIDASim")
  library("PALM")
  library("phyloseq")
  library("glmnet")
  library("abess")
  library("mclust")
  library("sClust")
  library("sparcl")
  library("cluster")
  library("tidyverse")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5)
    stop("need <s> <Ka.per> <pos.pt> <u> <scenario>, e.g.\n",
         "  Rscript Analysis/1_simulation_stage1.R 1 0.4 0.5 0 3")
  print(args)

  s        <- as.numeric(args[1])
  Ka.per   <- as.numeric(args[2])
  pos.pt   <- as.numeric(args[3])
  u        <- as.numeric(args[4])
  scenario <- as.numeric(args[5])

  ## ---- paths, relative to the repository root --------------------------------
  if (!dir.exists("utility"))
    stop("run this from the project root, e.g. Rscript Analysis/1_simulation_stage1.R ...")

  out_dir <- "./Simulation/Sim_CRC_stage1"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  data.loc <- file.path(out_dir, paste0("res_Ka", Ka.per, "_pos", pos.pt, "_u", u,
                                        "_scenario", scenario, "_s", s, ".Rdata"))

  ## Template: the colorectal strata object, whichever copy this checkout has.
  template <- c("./Simulation/Data/CRC.Rdata", "./CRC_analysis/Data/CRC_strata.Rdata")
  template <- template[file.exists(template)]
  if (!length(template)) stop("no CRC template found under ./Data/")
  load(template[1])
  message("template: ", template[1])

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
  
  SMESH_mod <- smesh.meta.summary(summary.stats = summary.stat.obj, 
                                 nperm = 5,
                                 doc = "./",
                                 verbose = TRUE)
  
  ## ---- summarize the results
  est_clust <- apply(round(SMESH_mod$disease$W), 1, function(d){which(d == 1)})
  ARI <- c(ARI, adjustedRandIndex(est_clust, G.lst))
  G_hat <- c(G_hat, sum(colSums(round(SMESH_mod$disease$W)) > 0))
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- G.lst[d]
    g2 <- which.max(SMESH_mod$disease$W[d,])
    truth.signal <- signal.lst[[G.lst[d]]]
    selected.signal <- rownames(SMESH_mod$disease$mu)[SMESH_mod$disease$mu[,g2]!=0]
    pre <- c(pre, length(intersect(truth.signal, selected.signal))/length(selected.signal))
    rec <- c(rec, length(intersect(truth.signal, selected.signal))/length(truth.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre, na.rm = TRUE))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))

  ## ---- homogeneous part
  detect_hom <- names(which(apply(SMESH_mod$disease$mu, 1, function(d){all(d !=0)})))
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
    g2 <- which.max(SMESH_mod$disease$W[d,])
    selected.signal <- signal.lst.het[[g]]
    detected.signal <- setdiff(rownames(SMESH_mod$disease$mu)[SMESH_mod$disease$mu[,g2]!=0], detect_hom)
    pre <- c(pre, length(intersect(detected.signal, selected.signal)) / length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal)) / length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, mean(pre))
  Recall_het <- c(Recall_het, mean(rec))
  F1_het <- c(F1_het, mean(f1))

  ## True-cluster-FE
  source("./utility/PALM_tune.R")

  ARI <- c(ARI, NA)
  detect.signal <- NULL
  for(g in 1:G){
    palm.hom.sin <- palm_tune(summary.stats = summary.stat.obj[which(G.lst == g)])
    detect.signal <- cbind(detect.signal, palm.hom.sin$disease$coef)
  }

  ## ---- summarize the results
  G_hat <- c(G_hat, NA)
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- G.lst[d]
    truth.signal <- intersect(signal.lst[[G.lst[d]]], colnames(data.rel[[d]]$Y))
    selected.signal <- rownames(detect.signal)[detect.signal[,g]!=0]
    pre <- c(pre, length(intersect(truth.signal, selected.signal))/length(selected.signal))
    rec <- c(rec, length(intersect(truth.signal, selected.signal))/length(truth.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))

  ## ---- homogeneous part
  detect_hom <- names(which(apply(detect.signal, 1, function(d){all(d != 0)})))
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
    selected.signal <- signal.lst.het[[g]]
    detected.signal <- setdiff(rownames(detect.signal)[detect.signal[,g]!=0], detect_hom)
    pre <- c(pre, length(intersect(detected.signal, selected.signal)) / length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal)) / length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, mean(pre))
  Recall_het <- c(Recall_het, mean(rec))
  F1_het <- c(F1_het, mean(f1))

  ## SKM + FE
  source("./utility/PALM_tune.R")

  feature_ID <- unique(unlist(sapply(summary.stat.obj, function(d){rownames(d$est)})))
  beta_summary <- matrix(NA, nrow = length(feature_ID), ncol = length(data.rel),
                         dimnames = list(feature_ID, names(data.rel)))
  for(d in names(data.rel)){
    beta_summary[rownames(summary.stat.obj[[d]]$est),d] <- summary.stat.obj[[d]]$est / summary.stat.obj[[d]]$stderr
  }
  beta_summary[is.na(beta_summary)] <- 0
  X <- t(beta_summary) # studies x taxa

  # choose a reasonable upper bound for K
  # at least 2-3 samples per cluster is a practical minimum
  Kmax_user <- 8
  Kmax_safe <- min(Kmax_user, floor(nrow(X) / 2))
  Kmax_safe <- max(2, Kmax_safe)
  K_for_w <- min(4, Kmax_safe)

  perm_skm <- KMeansSparseCluster.permute(
    X,
    K = K_for_w,
    nperms = 50,
    silent = TRUE
  )

  w_skm <- perm_skm$bestw
  B <- 200
  Khat_skm <- G

  # optional: visualize
  if (Khat_skm <= 1) {
    fit_skm <- NULL
    cl_skm <- rep(1L, nrow(X))
    ws_skm <- NULL
  } else {
    fit_skm <- KMeansSparseCluster(X, K = Khat_skm, wbounds = w_skm)
    cl_skm <- fit_skm[[1]]$Cs
    ws_skm <- fit_skm[[1]]$ws
  }
  KM_c <- cl_skm
  names(KM_c) <- names(data.rel)
  ARI <- c(ARI, adjustedRandIndex(KM_c, G.lst))

  # PALM-HOM model
  detect.signal <- NULL
  for(g in sort(unique(KM_c))){
    palm.hom.sin <- palm_tune(summary.stats = summary.stat.obj[which(KM_c == g)])
    detect.signal <- cbind(detect.signal, palm.hom.sin$disease$coef)
  }

  ## ---- summarize the results
  G_hat <- c(G_hat, NA)
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- KM_c[d]
    truth.signal <- intersect(signal.lst[[G.lst[d]]], colnames(data.rel[[d]]$Y))
    selected.signal <- rownames(detect.signal)[detect.signal[,g]!=0]
    pre <- c(pre, length(intersect(truth.signal, selected.signal))/length(selected.signal))
    rec <- c(rec, length(intersect(truth.signal, selected.signal))/length(truth.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))

  ## ---- homogeneous part
  detect_hom <- names(which(apply(detect.signal, 1, function(d){all(d != 0)})))
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
    g <- KM_c[d]
    selected.signal <- signal.lst.het[[G.lst[d]]]
    detected.signal <- setdiff(rownames(detect.signal)[detect.signal[,g]!=0], detect_hom)
    pre <- c(pre, length(intersect(detected.signal, selected.signal)) / length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal)) / length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, mean(pre))
  Recall_het <- c(Recall_het, mean(rec))
  F1_het <- c(F1_het, mean(f1))

  ## SHC + FE
  source("./utility/PALM_tune.R")

  feature_ID <- unique(unlist(sapply(summary.stat.obj, function(d){rownames(d$est)})))
  beta_summary <- matrix(NA, nrow = length(feature_ID), ncol = length(data.rel),
                         dimnames = list(feature_ID, names(data.rel)))
  for(d in names(data.rel)){
    beta_summary[rownames(summary.stat.obj[[d]]$est),d] <- summary.stat.obj[[d]]$est / summary.stat.obj[[d]]$stderr
  }
  beta_summary[is.na(beta_summary)] <- 0
  X <- t(beta_summary) # studies x taxa

  Kmax_user <- 10
  Kmax_safe <- min(Kmax_user, floor(nrow(X) / 2))
  Kmax_safe <- max(2, Kmax_safe)
  B <- 200
  perm_shc <- HierarchicalSparseCluster.permute(
    X,
    nperms = 50
  )
  w_shc <- perm_shc$bestw

  Khat_shc <- G

  fit_shc <- HierarchicalSparseCluster(X, wbound = w_shc)
  if (Khat_shc == 1) {
    cl_shc <- rep(1L, nrow(X))
    ws_shc <- NULL
  } else {
    cl_shc <- cutree(fit_shc$hc, k = Khat_shc)
    ws_shc <- fit_shc$ws
  }
  KM_c <- cl_shc
  names(KM_c) <- names(data.rel)
  ARI <- c(ARI, adjustedRandIndex(KM_c, G.lst))

  ## PALM-HOM model
  detect.signal <- NULL
  for(g in sort(unique(KM_c))){
    palm.hom.sin <- palm_tune(summary.stats = summary.stat.obj[which(KM_c == g)])
    detect.signal <- cbind(detect.signal, palm.hom.sin$disease$coef)
  }

  ## ---- summarize the results
  G_hat <- c(G_hat, NA)
  pre = rec <- NULL
  for(d in names(data.rel)){
    g <- KM_c[d]
    truth.signal <- intersect(signal.lst[[G.lst[d]]], colnames(data.rel[[d]]$Y))
    selected.signal <- rownames(detect.signal)[detect.signal[,g]!=0]
    pre <- c(pre, length(intersect(truth.signal, selected.signal))/length(selected.signal))
    rec <- c(rec, length(intersect(truth.signal, selected.signal))/length(truth.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))

  ## ---- homogeneous part
  detect_hom <- names(which(apply(detect.signal, 1, function(d){all(d != 0)})))
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
    g <- KM_c[d]
    selected.signal <- signal.lst.het[[G.lst[d]]]
    detected.signal <- setdiff(rownames(detect.signal)[detect.signal[,g]!=0], detect_hom)
    pre <- c(pre, length(intersect(detected.signal, selected.signal)) / length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal)) / length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, mean(pre))
  Recall_het <- c(Recall_het, mean(rec))
  F1_het <- c(F1_het, mean(f1))

  ## Melody
  ARI <- c(ARI, NA)
  null.obj2 <- miMeta::melody.null.model(rel.abd = rel.abd, prev.filter = 0, parallel.core = NULL)

  summary.stat.obj2 <- miMeta::melody.get.summary(null.obj = null.obj2,
                                                  covariate.interest = covariate.interest,
                                                  parallel.core = NULL)

  Melody_mod <- miMeta::melody.meta.summary(summary.stats = summary.stat.obj2, output.best.one = TRUE)

  ## ---- summarize the results
  G_hat <- c(G_hat, NA)
  pre = rec <- NULL
  for(d in names(rel.abd)){
    g <- G.lst[d]
    selected.signal <- signal.lst[[g]]
    detected.signal <- names(Melody_mod$disease$coef)[Melody_mod$disease$coef!=0]
    pre <- c(pre, length(intersect(detected.signal, selected.signal))/length(detected.signal))
    rec <- c(rec, length(intersect(detected.signal, selected.signal))/length(selected.signal))
  }
  f1 <- 2 / (1/pre + 1/rec)
  Precision <- c(Precision, mean(pre))
  Recall <- c(Recall, mean(rec))
  F1 <- c(F1, mean(f1))

  ## ---- homogeneous part
  Precision_hom <- c(Precision_hom, mean(pre))
  Recall_hom <- c(Recall_hom, mean(pre))
  F1_hom <- c(F1_hom, mean(pre))

  ## ---- heterogeneous part
  f1 <- 2 / (1/pre + 1/rec)
  Precision_het <- c(Precision_het, NA)
  Recall_het <- c(Recall_het, NA)
  F1_het <- c(F1_het, NA)

  ## ---- output
  result_mat <- data.frame(
    Precision = Precision, Recall = Recall, F1 = F1,
    Precision_hom = Precision_hom, Recall_hom = Recall_hom, F1_hom = F1_hom,
    Precision_het = Precision_het, Recall_het = Recall_het, F1_het = F1_het,
    G_hat = G_hat, ARI = ARI,
    method = c("SMESH", "True-cluster FE", "SKM + FE", "SHC + FE", "Melody")
  )
  
  save(result_mat, file = data.loc)
  
  
