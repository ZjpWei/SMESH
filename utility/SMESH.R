library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(Matrix)

## G: number of context clusters.
##   NULL (default) - search upward from G = 1 and stop at the first G whose
##     relative GIC improvement over the previous G falls short of GIC_REL_TOL,
##     keeping the previous G.  The search is capped at round(L / 3) so that a
##     cluster has, on average, at least three contexts to be estimated from.
##   a number        - fit that G directly, with no search.  Must not exceed the
##     number of contexts, since each cluster needs at least one.
smesh.meta.summary <- function(summary.stats,
                              G = NULL,
                              tune.type = c("BIC", "HBIC", "KBIC", "EBIC"),
                              tol = 1e-3,
                              NMAX = 20,
                              nperm = 5,
                              doc = "./",
                              verbose = FALSE) {

  tune.type <- match.arg(tune.type)

  ## Relative GIC gain a larger G must deliver to be accepted.
  GIC_REL_TOL <- 2e-2

  study.ID <- names(summary.stats)
  cov.int.ID <- sort(unique(unlist(
    lapply(summary.stats, function(d) colnames(d$est))
  )))
  
  feature.ID <- sort(unique(unlist(
    lapply(summary.stats, function(d) {
      rownames(d$est)[rowSums(!is.na(d$est)) > 0]
    })
  )))
  
  output.result <- vector("list", length(cov.int.ID))
  names(output.result) <- cov.int.ID
  
  for (cov.name in cov.int.ID) {
    if (verbose) {
      message("++ Search for the best model for covariate of interest: ", cov.name, " ++")
    }
    summary.stat.study <- list()
    for (d in study.ID) {
      est.mat <- summary.stats[[d]]$est
      se.mat  <- summary.stats[[d]]$stderr
      
      if (!cov.name %in% colnames(est.mat)) {
        next
      }
      nonna.id <- intersect(feature.ID, rownames(est.mat)[!is.na(est.mat[, cov.name])]) 
      if (length(nonna.id) == 0) {
        next
      }
      
      summary.stat.study[[d]] <- list(
        est = est.mat[nonna.id, cov.name],
        cov = diag(se.mat[nonna.id, cov.name]^2),
        n   = summary.stats[[d]]$n
      )
      colnames(summary.stat.study[[d]]$cov) = rownames(summary.stat.study[[d]]$cov) <- nonna.id
    }
    
    if (length(summary.stat.study) == 0) {
      warning("No valid study found for covariate: ", cov.name)
      next
    }
    
    if(verbose){
      message("++ Stage 1: study clustering with the original SMESH atom dictionary. ++")
    }

    ## ---- how many clusters to try ------------------------------------------
    L <- length(summary.stat.study)
    if (is.null(G)) {
      ## max(1, .) so a study set too small to divide still fits G = 1 rather
      ## than producing the descending sequence 1:0.
      G.grid <- seq_len(max(1, round(L / 3)))
    } else {
      if (length(G) != 1 || is.na(G) || G != round(G) || G < 1) {
        stop("G must be NULL or a single positive integer.")
      }
      if (G > L) {
        stop("G = ", G, " exceeds the number of contexts available for '",
             cov.name, "' (", L, "); each cluster needs at least one context.")
      }
      G.grid <- G
    }

    ## ---- stage 1, over the grid --------------------------------------------
    ## G0 / GIC0 / model0 hold the best accepted fit so far.  A larger G is only
    ## accepted if it improves the GIC by more than GIC_REL_TOL in relative
    ## terms; the first failure ends the search, so the grid is walked upward
    ## and never revisited.  With an explicit G the grid has one element and the
    ## loop simply fits it.
    G0 <- GIC0 <- model0 <- NULL
    for (G.try in G.grid) {
      if(verbose){
        message("++ Search for the best model for G = ", G.try, ". ++")
      }

      ## get lasso mat.  The dictionary must match the one meta_stage1() builds:
      ## a single atom at G = 1, otherwise one atom per cluster plus a shared one.
      lasso.mat <- get_lasso_pre(summary.stat.study = summary.stat.study,
                                 G = G.try,
                                 D_mat = if (G.try == 1) matrix(1, 1, 1)
                                         else rbind(diag(G.try), rep(1, G.try)),
                                 feature.ID = feature.ID,
                                 feature.set = lapply(summary.stat.study, function(d){
                                   taxa.vec.tmp <- rep(FALSE, length(feature.ID))
                                   names(taxa.vec.tmp) <- feature.ID
                                   taxa.vec.tmp[names(d$est)[!is.na(d$est)]] <- TRUE
                                   return(taxa.vec.tmp)
                                 }),
                                 study.ID = study.ID)

      model_stage1 <- meta_stage1(
        summary.stat.study = summary.stat.study,
        lasso.mat = lasso.mat,
        G = G.try,
        feature.ID = feature.ID,
        tune.type = tune.type,
        NMAX = NMAX,
        tol = tol,
        nperm = nperm,
        doc = doc,
        verbose = verbose
      )

      if (is.null(model0)) {
        ## first G on the grid: nothing to compare against yet
        G0 <- G.try; GIC0 <- model_stage1$gic; model0 <- model_stage1
      } else if ((GIC0 - model_stage1$gic) / GIC0 < GIC_REL_TOL) {
        if(verbose){
          message("++ G = ", G.try, " does not improve the GIC by more than ",
                  GIC_REL_TOL * 100, "%; keeping G = ", G0, ". ++")
        }
        break
      } else {
        G0 <- G.try; GIC0 <- model_stage1$gic; model0 <- model_stage1
      }
    }

    ## stage 2
    if(G0 > 2){
      if(verbose){
        message("++ Stage 2: fixed-cluster subset-basis refit for final signature estimation. ++")
      }

      ## preliminary D matrix
      generate_D_mat <- function(G) {
        D_mat <- do.call(rbind, lapply(seq_len(G), function(k) {
          combs <- combn(G, k)
          t(apply(combs, 2, function(idx) {
            x <- rep(0, G)
            x[idx] <- 1
            x
          }))
        }))
        rownames(D_mat) <- NULL
        colnames(D_mat) <- paste0("g", seq_len(G))
        return(D_mat)
      }
      
      lasso.mat <- get_lasso_pre(summary.stat.study = summary.stat.study,
                                 G = G0,
                                 D_mat = generate_D_mat(G0),
                                 feature.ID = feature.ID,
                                 feature.set = lapply(summary.stat.study, function(d){
                                   taxa.vec.tmp <- rep(FALSE, length(feature.ID))
                                   names(taxa.vec.tmp) <- feature.ID
                                   taxa.vec.tmp[names(d$est)[!is.na(d$est)]] <- TRUE
                                   return(taxa.vec.tmp)
                                 }),
                                 study.ID = study.ID)

      model_stage2 <- meta_stage2(
        summary.stat.study = summary.stat.study,
        lasso.mat = lasso.mat,
        tune.type = tune.type,
        D_mat = generate_D_mat(G0),
        w_lg = t(apply(model0$W, 1, function(d){setNames(as.numeric(d == max(d)), seq_along(d))})),
        s_start = model0$s_start
      )

      if(verbose){
        message("++ Stage 2: post-selection refit. ++")
      }
      
      ## Post-selection refit
      theta_refit <- t(sapply(rownames(model_stage2$mu),
                              function(sk,
                                       W = t(apply(model_stage2$W, 1, function(d){setNames(as.numeric(d == max(d)), seq_along(d))})),
                                       B = model_stage2$B,
                                       D_mat = generate_D_mat(G0)){

                                M_sub <- B[sk,]
                                if(all(M_sub == 0)){
                                  return(rep(0, ncol(W)))
                                }else{
                                  D_mat_sub <- D_mat[M_sub!=0,,drop=FALSE]

                                  beta_k <- matrix(sapply(summary.stats, function(d, k = sk){
                                    if(k %in% rownames(d$est)){
                                      return(d$est[k,1])
                                    }else{
                                      return(NA)
                                    }
                                  }), nrow = nrow(W), ncol = ncol(W), dimnames = list(rownames(W), colnames(W)))

                                  v_k <- matrix(sapply(summary.stats, function(d, k = sk){
                                    if(k %in% rownames(d$est)){
                                      return(d$stderr[k,1]^2)
                                    }else{
                                      return(NA)
                                    }
                                  }), nrow = nrow(W), ncol = ncol(W), dimnames = list(rownames(W), colnames(W)))

                                  A_k <- diag(colSums(W / v_k, na.rm = TRUE))
                                  C_k <- colSums(W * beta_k / v_k, na.rm = TRUE)

                                  return(t(D_mat_sub) %*% MASS::ginv(D_mat_sub %*% A_k %*% t(D_mat_sub)) %*% D_mat_sub %*% C_k)
                                }
                              }))

      model_stage2$mu_refit <- theta_refit
      model_stage2$mu_stage1 <- model0$mu
      output.result[[cov.name]] <- model_stage2
    }else{
      output.result[[cov.name]] <- model0
    }
  }
  
  return(output.result)
}

get_lasso_pre <- function(summary.stat.study, 
                          G,
                          D_mat,
                          feature.ID, 
                          feature.set, 
                          study.ID){
  
  L <- length(study.ID)
  if(L == 1){
    tmp.taxa.mat <- feature.set[[study.ID]]
    taxa.mat <- matrix(tmp.taxa.mat, nrow = 1)
    colnames(taxa.mat) <- names(tmp.taxa.mat)
  }else{
    taxa.mat <- NULL
    for(d in study.ID){
      tmp.taxa.mat <- feature.set[[d]]
      taxa.mat <- rbind(taxa.mat, tmp.taxa.mat)
    }
  }
  rownames(taxa.mat) <- study.ID
  nonempty.id <- names(which(colSums(taxa.mat) > 0))
  empty.id <- names(which(colSums(taxa.mat) == 0))
  mu.len <- length(nonempty.id)
  taxa.mat <- taxa.mat[,nonempty.id]
  if(L == 1){
    taxa.mat <- matrix(taxa.mat, nrow = 1)
    colnames(taxa.mat) <- nonempty.id
    rownames(taxa.mat) <- study.ID
  }
  N <- sum(sapply(summary.stat.study, function(d){d$n}))
  k.l <- apply(taxa.mat, 1, sum)
  
  if(G == 1){
    D_sparse <- Matrix::Matrix(1, sparse = TRUE)
  }else{
    D_sparse <- Matrix::Matrix(t(D_mat), sparse = TRUE)
  }
  
  Y_list <- list()
  R_mat <- matrix(0, nrow = L, ncol = length(feature.ID), dimnames = list(study.ID, feature.ID))
  for (d in study.ID) {
    loc.k <- which(d == study.ID)
    Sigma.chol <- Matrix::Diagonal(
      x = sqrt(1 / diag(summary.stat.study[[d]]$cov))
    )
    
    ## Record R matrix
    R_mat[d, names(summary.stat.study[[d]]$est)] <- 1
    
    ## New Y_l = 1_G %x% weighted beta_hat_l
    y_l <- Sigma.chol %*% summary.stat.study[[d]]$est
    Y_list[[loc.k]] <- as.vector(rep(y_l, times = G))
  }
  Y.enlarge <- unlist(Y_list)
  
  ## Number of coefficients in the enlarged design, i.e. ncol(D %x% A_l).
  ## byabess_cluster rebuilds that design itself (with weight-scaled atoms), so only
  ## the count is needed here - it sets the upper end of the sparsity search.
  n.coef <- ncol(D_sparse) * mu.len
  
  lasso.mat <- list(Y.enlarge = Y.enlarge, 
                    n.coef = n.coef,
                    R_mat = R_mat,
                    taxa.mat = taxa.mat,
                    k.l = k.l, 
                    N = N, 
                    mu.len = mu.len,
                    nonempty.id = nonempty.id, 
                    empty.id = empty.id,  
                    feature.ID = feature.ID,
                    study.ID = study.ID)
  
  return(lasso.mat)
}

# Optimization function
byabess_cluster <- function(summary.stat.study, 
                       D_mat,
                       w_lg, 
                       lasso.mat, 
                       support.size,
                       tune.type = c("BIC", "HBIC", "KBIC", "EBIC"),
                       stage = c("1", "2")){
  
  G <- ncol(D_mat)
  M <- nrow(D_mat)
  
  ## customize the X matrix
  W_vec <- list()
  for (d in lasso.mat$study.ID) {
    loc.k <- which(d == lasso.mat$study.ID)
    
    ## W_lst
    W_vec[[loc.k]] <- rep(w_lg[d,], each = length(summary.stat.study[[d]]$est))
  }
  W_hat <- unlist(W_vec)
  
  X_list <- list()
  ## Each atom is scaled by the sqrt of the posterior mass of the clusters it
  ## spans.  DO NOT drop the 1e-5: if a cluster empties during EM its mass is
  ## exactly 0, and an unguarded division puts Inf into the design.  glmtlp
  ## accepted that silently; abess validates its input and stops with
  ## "x has missing value or infinite value!".  This mirrors the guard in
  ## SMESH_v4.R, which divided by (sqrt(pi_h[g]) + 1e-5).
  atom.mass <- sqrt(sapply(1:nrow(D_mat), function(d){
    sum(colMeans(w_lg[, D_mat[d, ] == 1, drop = FALSE]))
  }))
  D_sparse <- Matrix::Matrix(t(D_mat / (atom.mass + 1e-5)), sparse = TRUE)
  for (d in lasso.mat$study.ID) {
    loc.k <- which(d == lasso.mat$study.ID)
    Sigma.chol <- Matrix::Diagonal(
      x = sqrt(1 / diag(summary.stat.study[[d]]$cov))
    )
    mu.mat <- Matrix::Diagonal(lasso.mat$mu.len)
    id.mu <- lasso.mat$taxa.mat[d, ]
    mu.mat <- mu.mat[id.mu, , drop = FALSE]
    A_l <- Sigma.chol %*% mu.mat
    
    ## New X_l = D %x% A_l
    X_list[[loc.k]] <- kronecker(D_sparse, A_l)
    
  }
  X.enlarge <- as.matrix(do.call(rbind, X_list))
  
  ## optim model
  suppressMessages(
    result <- abess::abess(x = X.enlarge,
                           y = as.vector(lasso.mat$Y.enlarge),
                           weight = W_hat,
                           family = "gaussian",
                           tune.type = "bic",
                           tune.path = "sequence",
                           normalize = 0,
                           support.size = support.size,
                           fit.intercept = FALSE,
                           standardize = FALSE)
    
    # result <- glmtlp::glmtlp(X = X.enlarge,
    #                          y = as.vector(lasso.mat$Y.enlarge), 
    #                          family = "gaussian",
    #                          penalty = "l0",
    #                          weights = W_hat,
    #                          kappa = support.size,
    #                          tol = 1e-3,      ## Default: 1e-4
    #                          dc.maxit = 10,   ## Default: 20
    #                          cd.maxit = 1000, ## Default: 10000
    #                          nr.maxit = 10,   ## Default: 20
    #                          standardize = FALSE)
  )
  
  ## Compute residuals and q_loss
  res_beta <- result$beta[,1]
  resid <- as.vector(lasso.mat$Y.enlarge) - X.enlarge %*% res_beta
  q_loss <- sum(W_hat * resid^2) / lasso.mat$N
  
  ## Re-construct estimates
  B_mat <- matrix(res_beta, nrow = lasso.mat$mu.len, ncol = M, dimnames = list(lasso.mat$nonempty.id, as.character(1:M)))
  theta_mat <- B_mat %*% D_mat 
  
  ## Calculate ic
  if(stage == "1"){
    df <- sum(res_beta!=0) + (G - 1)
  }else{
    df <- sum(res_beta!=0)
  }
  
  if(tune.type == "BIC"){
    ic <- df * log(lasso.mat$N)/lasso.mat$N
  }else if(tune.type == "KBIC"){
    ic <- df * log(max(exp(1), log(length(res_beta)))) * log(lasso.mat$N) / lasso.mat$N
  }else if(tune.type == "HBIC"){
    ic <- df * log(length(res_beta)) * log(log(lasso.mat$N))/lasso.mat$N
  }else if (tune.type == "EBIC"){
    ic <- (df * log(lasso.mat$N) + (lfactorial(length(res_beta)) - lfactorial(df) - lfactorial(length(res_beta) - df))) / lasso.mat$N
  }
  return(list(theta_mat = theta_mat, mu = B_mat, res_beta = res_beta, q_loss = q_loss, ic = ic, GIC = q_loss + ic))
}

# Meta-analysis
meta_stage1 <- function(summary.stat.study, 
                        lasso.mat,
                        tune.type, 
                        G, 
                        feature.ID, 
                        tol, 
                        NMAX, 
                        nperm,
                        doc,
                        verbose){
  
  study.ID <- names(summary.stat.study)
  L <- length(study.ID)

  ## Stage-1 atom dictionary: one atom per cluster plus one all-cluster-shared
  ## atom.  This must match the D_mat that smesh.meta.summary() passed to
  ## get_lasso_pre() when it built `lasso.mat`, or byabess_cluster() will index the
  ## design matrix with the wrong number of atoms.
  ## At G > 1 the dictionary is one atom per cluster plus one shared atom.  At
  ## G = 1 those two collapse to the same column, so a single atom is used
  ## instead - which is also what get_lasso_pre() builds in its G == 1 branch.
  D_mat <- if (G == 1) matrix(1, 1, 1) else rbind(diag(G), rep(1, G))

  ## Progress reporting.  Each sparsity level runs nperm restarts of an EM loop,
  ## which is far too many lines to print one by one.  The counter is therefore
  ## rewritten in place with a carriage return, so a whole sparsity level
  ## occupies a single live line; when the level finishes that line is cleared
  ## and one permanent line is committed.
  ##
  ## Padded to a fixed width and cleared with spaces rather than an ANSI escape,
  ## so a redirected log picks up no control codes beyond the "\r" itself.
  PROG_W    <- 56
  prog_line <- function(s, sd, nperm, shift) {
    formatC(sprintf("      [s = %s | restart %d/%d] pi shift %s",
                    s, sd, nperm, signif(shift, 3)),
            width = -PROG_W)
  }
  prog_clear <- function() {
    cat("\r", strrep(" ", PROG_W), "\r", sep = "", file = stderr())
  }

  ## Start to tune
  loop.search <- TRUE
  if(G == 1){
    search.vec <- round( lasso.mat$n.coef * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0) / 2 )
  }else{
    search.vec <- round( lasso.mat$n.coef * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0) / (G+1) )
  }
  search.vec.loc <- rep(TRUE, 4)
  initial_sort <- list()
  gic_sort <- rep(NA, 4)
  
  while(loop.search){
    for(k in 1:sum(search.vec.loc)){
      s.lambda <- (search.vec[search.vec.loc])[k]
      
      if(L == 1){
        tune.deltas <- 1
        names(tune.deltas) <- study.ID
      }else{
        
        cluster_list <- list()
        for(sd in 1:nperm){
          s.lambda.local <- s.lambda
          
          ## Random initial cluster
          set.seed(sd)
          tune.deltas <- sample(1:G, size = L, replace = TRUE)
          names(tune.deltas) <- study.ID
          
          ## initial matrix
          beta_mat <- matrix(
            NA,
            nrow = length(feature.ID),
            ncol = L,
            dimnames = list(feature.ID, study.ID)
          )
          for (d in study.ID) {
            beta_mat[names(summary.stat.study[[d]]$est), d] <- summary.stat.study[[d]]$est
          }
          
          ## initialize mu_g and pi_g
          mu_g <- pi_g <- NULL
          for (l in 1:G) {
            mu_g <- cbind(
              mu_g,
              rowMeans(
                beta_mat[, names(tune.deltas)[tune.deltas == l], drop = FALSE],
                na.rm = TRUE
              )
            )
            pi_g <- c(pi_g, mean(tune.deltas == l))
          }
          colnames(mu_g) <- names(pi_g) <- as.character(1:G)
          
          loops <- TRUE
          loop.count <- 0
          
          while (loops) {
            ## E step
            log_dens <- matrix(
              NA, nrow = L, ncol = G,
              dimnames = list(colnames(beta_mat), as.character(1:G))
            )
            
            for (g in 1:G) {
              for (d in colnames(beta_mat)) {
                mu.tmp <- mu_g[!is.na(beta_mat[, d]), g]
                mu.tmp[is.na(mu.tmp)] <- 0
                log_dens[d, g] <- mvn_density(
                  x = beta_mat[!is.na(beta_mat[, d]), d],
                  mu = mu.tmp,
                  Sigma = summary.stat.study[[d]]$cov
                )
              }
            }
            
            ## apply() returns a bare vector when G == 1, so t() would give a
            ## 1 x L matrix instead of L x G.  Rebuild the shape explicitly.
            w_raw <- apply(log_dens, 1, function(d) {
              log_weighted_normalize(log_vals = d, weights = pi_g)
            })
            w_lg <- matrix(w_raw, nrow = L, ncol = G, byrow = TRUE,
                           dimnames = list(rownames(log_dens), as.character(1:G)))
            
            ## M step
            pi_h <- colMeans(w_lg)
            s.result <- byabess_cluster(
              summary.stat.study = summary.stat.study,
              D_mat = D_mat,
              w_lg = w_lg,
              lasso.mat = lasso.mat,
              support.size = s.lambda.local,
              tune.type = tune.type,
              stage = "1"
            )

            ## theta_mat is B_mat %*% D_mat: the cluster-level effects
            ## (features x G).  `mu` from byabess_cluster() is the basis coefficient
            ## matrix (features x M) and is not what the EM step needs.
            mu_h <- s.result$theta_mat
            
            ## `sd` is the consensus restart (a different random initial cluster
            ## assignment); `loop.count` counts EM iterations within it.
            pi.shift <- sqrt(sum((pi_h - pi_g)^2, na.rm = TRUE))
            if (verbose) {
              cat("\r", prog_line(s.lambda.local, sd, nperm, pi.shift),
                  sep = "", file = stderr())
            }
            if (pi.shift <= tol) {
              #if (sum(s.result$res_beta != 0) == s.lambda.local) {
              loops <- FALSE
              pi_g <- pi_h
              mu_g <- mu_h
            } else if (loop.count > NMAX) {
              warning(sprintf(
                "EM algorithm doesn't converge with given tolerance in s = %s.",
                s.lambda.local
              ))
              loops <- FALSE
              pi_g <- pi_h
              mu_g <- mu_h
            } else {
              pi_g <- pi_h
              mu_g <- mu_h
              loop.count <- loop.count + 1
            }
          }
          
          cluster_list[[sd]] <-   apply(w_lg, 1, which.max)
        }

        ## All restarts for this sparsity level are done: clear the live line
        ## and commit one permanent line showing where it ended.
        if (verbose) {
          prog_clear()
          message(prog_line(s.lambda.local, nperm, nperm, pi.shift))
        }

        # C: consensus matrix (square), row/colnames are study IDs
        ## 1) Hierarchical clustering on consensus (distance = 1 - C)
        C <- build_consensus(cluster_list)
        D  <- as.dist(1 - C)
        stopifnot(is.matrix(C), nrow(C) == ncol(C))
        stopifnot(identical(rownames(C), colnames(C)))
        
        ## Decide cluster way
        hc_ave <- hclust(D, method = "average")
        cluster0 <- cutree(hc_ave, k = G)
        s1 <- within_consensus_score(C, cluster0, weight = "size")
        
        hc_complete <- hclust(D, method = "complete")
        cluster1 <- cutree(hc_complete, k = G)
        s2 <- within_consensus_score(C, cluster1, weight = "size")
        
        if(s1 > s2){
          hc <- hc_ave
        }else{
           hc <- hc_complete
        }
        ord <- hc$labels[hc$order]      # ordered study names
        C_ord <- C[ord, ord, drop = FALSE]
        
        ## 2) Long format for ggplot
        if(verbose){
          df <- as.data.frame(C_ord) %>%
            rownames_to_column("i") %>%
            pivot_longer(-i, names_to = "j", values_to = "p") %>%
            mutate(
              i = factor(i, levels = rev(ord)),  # rev() so top-left is first in order
              j = factor(j, levels = ord)
            )
          
          ## 3) Plot clustered heatmap
          p <- ggplot(df, aes(x = j, y = i, fill = p)) +
            geom_tile(color = "white", linewidth = 0.15) +
            scale_fill_gradientn(
              colours = c("white", "yellow", "#FDB863", "red", "#B2182B"),
              limits = c(0, 1),
              oob = scales::squish,
              name = "Consensus\n(co-cluster)"
            ) +
            coord_fixed() +
            labs(
              title = "Consensus matrix (hierarchical clustering order)",
              x = NULL, y = NULL
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
              axis.text.y = element_text(size = 9),
              panel.grid = element_blank()
            )
        }
        
        
        #=== Search best model by sequence ===#
        tune.deltas <- cutree(hc, k = G)
        
        ## Same shape guard as above.  Width is G, not max(tune.deltas): cutree()
        ## returns labels 1..G so they agree, but G cannot silently under-count.
        w_raw <- sapply(tune.deltas, function(d){
          dd <- rep(0, G)
          dd[d] <- 1
          dd
        })
        w_lg <- matrix(w_raw, nrow = L, ncol = G, byrow = TRUE,
                       dimnames = list(names(tune.deltas), as.character(1:G)))
        
        ## Estimate:
        pi_h <- colMeans(w_lg)
        s.result <- byabess_cluster(summary.stat.study = summary.stat.study,
                               D_mat = D_mat,
                               w_lg = w_lg,
                               lasso.mat = lasso.mat,
                               support.size = s.lambda,
                               tune.type = tune.type,
                               stage = "1")

        mu_h <- s.result$theta_mat
      }
      
      ####
      ct <- which(search.vec == s.lambda)
      initial_sort[[ct]] <- list(pi_g = pi_h, mu_g = mu_h, w_lg = w_lg,
                                 gic = s.result$GIC, ic = s.result$ic,
                                 q_loss = s.result$q_loss,
                                 cluster_list = cluster_list,
                                 s.lambda = s.lambda)
      
      gic_sort[ct] <- initial_sort[[ct]]$gic
    }
    
    ## gsection
    min.gic.id <- which.min(gic_sort)[1]
    if(min.gic.id == 1 | min.gic.id == 2){
      if(abs(search.vec[1] - search.vec[2]) <= 2){
        loop.search <- FALSE
        results.all <- list(mu = initial_sort[[2]]$mu_g, 
                            W = initial_sort[[2]]$w_lg, 
                            Pi = initial_sort[[2]]$pi_g, 
                            gic = initial_sort[[2]]$gic,
                            ic = initial_sort[[2]]$ic,
                            q_loss = initial_sort[[2]]$q_loss,
                            cluster_list = initial_sort[[2]]$cluster_list,
                            ## Support size stage 1 settled on.  meta_stage2()
                            ## takes this as `s_start`, the top of its own
                            ## sparsity search - the role `s_tune` played in the
                            ## original cluster_tune() implementation.
                            s_start = initial_sort[[2]]$s.lambda)
      }else{
        search.vec[3:4] <- search.vec[2:3]
        search.vec[2] <- search.vec[1] + round((1 - 2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[4]] <- initial_sort[[3]]
        initial_sort[[3]] <- initial_sort[[2]]
        initial_sort[[2]] <- NA
        gic_sort[3:4] <- gic_sort[2:3]
        gic_sort[2] <- NA
        search.vec.loc <- c(FALSE, TRUE, FALSE, FALSE)     
      }
    }else{
      if(abs(search.vec[3] - search.vec[4]) <= 2){
        loop.search <- FALSE
        results.all <- list(mu = initial_sort[[3]]$mu_g, 
                            W = initial_sort[[3]]$w_lg, 
                            Pi = initial_sort[[3]]$pi_g, 
                            gic = initial_sort[[3]]$gic,
                            ic = initial_sort[[3]]$ic,
                            q_loss = initial_sort[[3]]$q_loss,
                            cluster_list = initial_sort[[3]]$cluster_list,
                            ## Support size stage 1 settled on.  meta_stage2()
                            ## takes this as `s_start`, the top of its own
                            ## sparsity search - the role `s_tune` played in the
                            ## original cluster_tune() implementation.
                            s_start = initial_sort[[3]]$s.lambda)
      }else{
        search.vec[1:2] <- search.vec[2:3]
        search.vec[3] <- search.vec[1] + round((2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[1]] <- initial_sort[[2]]
        initial_sort[[2]] <- initial_sort[[3]]
        initial_sort[[3]] <- NA
        gic_sort[1:2] <- gic_sort[2:3]
        gic_sort[3] <- NA
        search.vec.loc <- c(FALSE, FALSE, TRUE, FALSE)
      }
    }
  }
  if(verbose){
    ggsave(
      filename = paste0(doc, "consensus_matrix", G,".pdf"),
      plot = p,
      device = cairo_pdf,
      width = 10,
      height = 10,
      units = "in"
    )
  }
  return(results.all)
}

## Stage 2: the cluster assignment from stage 1 is held fixed (w_lg is a hard
## 0/1 membership matrix), so there is no EM loop here - only a golden-section
## search over the sparsity level on the 2^G - 1 subset basis.  That is why this
## function takes no tol/NMAX.
##
## `s_start` sets the upper end of that search: the support size stage 1
## selected, carried over in model_stage1$s_start.  This is the same quantity
## the original cluster_tune() returned as `s_tune`.
meta_stage2 <- function(summary.stat.study,
                        lasso.mat,
                        tune.type,
                        D_mat,
                        w_lg,
                        s_start){


  # Start to tune
  search.vec <- round(s_start * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0))
  search.vec.loc <- rep(TRUE, 4)
  
  ## Step 2: 2^G -1 model fit
  loop.search <- TRUE
  initial_sort <- list()
  gic_sort <- rep(NA, 4)
  while(loop.search){
    for(k in 1:sum(search.vec.loc)){
      s.lambda <- (search.vec[search.vec.loc])[k]
      
      ## M step:
      pi_h <- colMeans(w_lg)
      s.result <- byabess_cluster(summary.stat.study = summary.stat.study,
                             D_mat = D_mat,
                             w_lg = w_lg,
                             lasso.mat = lasso.mat,
                             support.size = s.lambda,
                             tune.type = tune.type,
                             stage = "2")
      
      theta_h <- s.result$theta_mat
      mu_h <- s.result$mu
      
      ## Record model:
      ct <- which(search.vec == s.lambda)
      initial_sort[[ct]] <- list(pi_g = pi_h,
                                 mu_g = mu_h,
                                 theta_g = theta_h, 
                                 w_lg = w_lg,
                                 gic = s.result$GIC,
                                 ic = s.result$ic, 
                                 q_loss = s.result$q_loss,
                                 s.lambda = s.lambda)
      
      gic_sort[ct] <- initial_sort[[ct]]$gic
    }
    
    ## Golden-section search step
    min.gic.id <- which.min(gic_sort)
    extract_result <- function(fit) {
      list(
        mu = fit$theta_g,
        B = fit$mu_g,
        W = fit$w_lg,
        Pi = fit$pi_g,
        gic = fit$gic,
        ic = fit$ic,
        q_loss = fit$q_loss,
        s.lambda = fit$s.lambda
      )
    }
    left.side <- min.gic.id %in% c(1, 2)
    
    if (left.side) {
      ## The minimum is on the left side of the current search interval.
      ## Shrink the interval from the right.
      if (abs(search.vec[1] - search.vec[2]) <= 2) {
        loop.search <- FALSE
        best.id <- which.min(gic_sort)
        results.all <- extract_result(initial_sort[[best.id]])
      } else {
        search.vec[3:4] <- search.vec[2:3]
        search.vec[2] <- search.vec[1] + round((1 - 2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[4]] <- initial_sort[[3]]
        initial_sort[[3]] <- initial_sort[[2]]
        initial_sort[[2]] <- NA
        gic_sort[3:4] <- gic_sort[2:3]
        gic_sort[2] <- NA
        search.vec.loc <- c(FALSE, TRUE, FALSE, FALSE)
      }
    } else {
      ## The minimum is on the right side of the current search interval.
      ## Shrink the interval from the left.
      if (abs(search.vec[3] - search.vec[4]) <= 2) {
        loop.search <- FALSE
        best.id <- which.min(gic_sort)
        results.all <- extract_result(initial_sort[[best.id]])
      } else {
        search.vec[1:2] <- search.vec[2:3]
        search.vec[3] <- search.vec[1] + round((2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[1]] <- initial_sort[[2]]
        initial_sort[[2]] <- initial_sort[[3]]
        initial_sort[[3]] <- NA
        gic_sort[1:2] <- gic_sort[2:3]
        gic_sort[3] <- NA
        search.vec.loc <- c(FALSE, FALSE, TRUE, FALSE)
      }
    }
  }
  
  return(results.all)
}
mvn_density <- function(x, mu, Sigma) {
  n <- length(mu)
  centered <- x - mu
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus
  mahalanobis_sq <- sum(centered^2 / diag(Sigma))
  log_density <- -0.5 * (n * log(2 * pi) + log_det_Sigma + mahalanobis_sq)
  return(as.numeric(log_density))
} 

log_weighted_normalize <- function(log_vals, weights) {
  if (length(log_vals) != length(weights)) {
    stop("log_vals and weights must have the same length")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    warning("Weights do not sum to 1 — normalizing them.")
    weights <- weights / sum(weights)
  }
  log_weighted <- log(weights) + log_vals
  log_max <- max(log_weighted)
  stabilized <- log_weighted - log_max
  unnormalized <- exp(stabilized)
  normalized_weights <- unnormalized / sum(unnormalized)
  return(normalized_weights)
}

build_consensus <- function(cluster_list){
  studies <- names(cluster_list[[1]])
  n <- length(studies)
  S <- length(cluster_list)
  C <- matrix(0, n, n, dimnames = list(studies, studies))
  for (s in seq_along(cluster_list)) {
    z <- cluster_list[[s]][studies]
    for (i in 1:n) {
      for (j in 1:n) {
        if (!is.na(z[i]) && !is.na(z[j])) {
          if (z[i] == z[j]) {
            C[i,j] <- C[i,j] + 1
          }
        }
      }
    }
  }
  C / S
}
within_consensus_score <- function(C, cluster, weight = c("size", "pairs")) {
  weight <- match.arg(weight)
  if (!is.matrix(C) || nrow(C) != ncol(C)) {
    stop("C must be a square matrix.")
  }
  if (length(cluster) != nrow(C)) {
    stop("length(cluster) must equal nrow(C).")
  }
  groups <- split(seq_along(cluster), cluster)
  avg_vals <- numeric(length(groups))
  weights <- numeric(length(groups))
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    nk <- length(idx)
    if (nk < 2) {
      avg_vals[k] <- NA
      weights[k] <- 0
      next
    }
    subC <- C[idx, idx, drop = FALSE]
    avg_vals[k] <- mean(subC[upper.tri(subC)])
    if (weight == "size") {
      weights[k] <- nk
    } else if (weight == "pairs") {
      weights[k] <- choose(nk, 2)
    }
  }
  keep <- weights > 0
  sum(weights[keep] * avg_vals[keep]) / sum(weights[keep])
}
