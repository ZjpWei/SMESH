# =============================================================================
#   SMESH.R  --  structure-learning meta-analysis across contexts
# =============================================================================
#  Implements the revised SMESH algorithms (Algorithms 1-3 of the supplement):
#
#    Algorithm 1  Stage 1 at fixed G and support size s, from ONE start.  The EM
#                 is initialised either from a random partition or, along the
#                 support path, from the same start's fit at a neighbouring s.
#                 Convergence is on the relative change in the observed-data
#                 negative log-likelihood.
#    Algorithm 2  cached discrete golden-section search over s; every support
#                 size in the terminal bracket is evaluated, and the returned
#                 size is the argmin over everything evaluated.
#    Algorithm 3  G by the global minimum of the Stage 1 GIC; Stage 2 over the
#                 sharing patterns allowed for each feature; unpenalized refit.
#
#  Random starts protect against local optima but are NOT averaged: at each
#  (G, s) the best converged start is retained, and the restart agreement is
#  reported only as an optimisation diagnostic.
#
#  OUTPUT CONTRACT.  The returned object keeps the field names the figure and
#  table scripts expect, so Rscript/Figure*.R and Rscript/Make_Table*.R run
#  unchanged:
#      mu           features x G   final cluster effects (unpenalized refit;
#                                   0 where a cluster cannot estimate the feature)
#      W            contexts x G   hard cluster assignment, 0/1
#      Pi           G              mixing proportions
#      gic, ic, q_loss             Stage 1 GIC at the selected (G, s) and its parts
#      B            features x M   selected basis coefficients
#      s.lambda                    selected support size
#      mu_refit     features x G   the same refit with NA where a cluster has no data
#      mu_stage1    features x G   Stage 1 effects before Stage 2
#      W_post       contexts x G   posterior memberships
#      cluster_list                per-start hard assignments at the selected (G, s)
#      cluster_gic                 the GIC each of those starts reached
#      gic_path                    Stage 1 GIC by G, and the support paths
#
#  The legacy consensus implementation that produced the published fits is kept
#  at tmp/legacy/SMESH_consensus_legacy.R.
# =============================================================================

suppressMessages({
  library(Matrix)
  library(abess)
  library(MASS)
})


# =============================================================================
#  Bases, data, likelihood
# =============================================================================

## Stage 1 patterns: the shared atom and one atom per cluster, as a G x M matrix
## D with theta (K x G) = B (K x M) %*% t(D).
stage1_basis <- function(G) if (G == 1) matrix(1, 1, 1) else cbind(1, diag(G))

## All 2^G - 1 nonempty subsets of clusters, one column per pattern.
all_subsets <- function(G) {
  pats <- unlist(lapply(seq_len(G), function(k) combn(G, k, simplify = FALSE)), recursive = FALSE)
  matrix(vapply(pats, function(p) { x <- numeric(G); x[p] <- 1; x }, numeric(G)), nrow = G)
}

## Context x feature matrices of estimates and variances for one covariate.
prep_data <- function(summary.stats, cov) {
  keep <- names(summary.stats)[vapply(summary.stats, function(d)
    cov %in% colnames(d$est) && any(!is.na(d$est[, cov])), logical(1))]
  feature <- sort(unique(unlist(lapply(summary.stats[keep], function(d)
    rownames(d$est)[!is.na(d$est[, cov]) & !is.na(d$stderr[, cov]) & d$stderr[, cov] > 0]))))
  L <- length(keep); K <- length(feature)
  beta <- v <- matrix(NA_real_, L, K, dimnames = list(keep, feature))
  for (l in keep) {
    e <- summary.stats[[l]]$est[, cov]; s <- summary.stats[[l]]$stderr[, cov]
    ok <- intersect(feature, names(e)[!is.na(e) & !is.na(s) & s > 0])
    beta[l, ok] <- e[ok]; v[l, ok] <- s[ok]^2
  }
  R <- !is.na(beta); Rn <- R * 1
  list(context = keep, feature = feature, L = L, K = K, beta = beta, v = v, R = R, Rn = Rn,
       mplus = colSums(Rn),
       n0 = sum(vapply(summary.stats[keep], function(d) d$n, numeric(1))),
       P = ifelse(R, 1 / v, 0), Yb = ifelse(R, beta / v, 0))
}

## log f(beta_l ; V_l, theta_g) for every context and cluster (L x G).
logdens <- function(dat, theta) {
  out <- vapply(seq_len(ncol(theta)), function(g) {
    r <- dat$beta - matrix(theta[, g], dat$L, dat$K, byrow = TRUE)
    -0.5 * rowSums(log(2 * base::pi * dat$v) + r^2 / dat$v, na.rm = TRUE)
  }, numeric(dat$L))
  matrix(out, nrow = dat$L)
}

## Responsibilities and the observed-data negative log-likelihood Q_obs.
estep <- function(dat, theta, pr) {
  lw <- sweep(logdens(dat, theta), 2, log(pr), "+")
  mx <- apply(lw, 1, max)
  ew <- exp(lw - mx); s <- rowSums(ew)
  list(W = matrix(ew / s, nrow = dat$L), Q = -sum(mx + log(s)))
}

## n0^-1 sum_l sum_g w_lg (beta_l - theta_g)' V_l^-1 (beta_l - theta_g).
wrss <- function(dat, W, theta) {
  s <- 0
  for (g in seq_len(ncol(theta))) {
    r <- dat$beta - matrix(theta[, g], dat$L, dat$K, byrow = TRUE)
    s <- s + sum(W[, g] * rowSums(r^2 / dat$v, na.rm = TRUE))
  }
  s / dat$n0
}

## Complexity penalty; BIC by default, with the three alternatives.
ic_pen <- function(df, p, N, type) switch(type,
  BIC  = df * log(N) / N,
  KBIC = df * log(max(exp(1), log(p))) * log(N) / N,
  HBIC = df * log(p) * log(log(N)) / N,
  EBIC = (df * log(N) + lchoose(p, df)) / N)

## Feature-specific calibration c_ka = { sum_{g in A_a} m_kg / m_k+ }^(1/2).
calib <- function(dat, W, D) sqrt((crossprod(dat$Rn, W) %*% D) / dat$mplus)

hard <- function(W) {
  Z <- matrix(0, nrow(W), ncol(W))
  Z[cbind(seq_len(nrow(W)), max.col(W, "first"))] <- 1
  Z
}

## A random partition of L contexts into G nonempty groups (surjective).
rand_partition <- function(L, G) {
  lab <- integer(L); perm <- sample.int(L)
  lab[perm[seq_len(G)]] <- sample.int(G)
  if (L > G) lab[perm[-seq_len(G)]] <- sample.int(G, L - G, replace = TRUE)
  Z <- matrix(0, L, G); Z[cbind(seq_len(L), lab)] <- 1; Z
}

## Algorithm 1, step 3: each cluster-feature effect is the inverse-variance
## weighted mean over the available contexts of that cluster, with a pooled
## feature-level fallback when the cluster has no data for the feature.
init_theta <- function(dat, Z) {
  num <- crossprod(Z, dat$Yb); den <- crossprod(Z, dat$P)
  pooled <- colSums(dat$Yb) / colSums(dat$P)
  th <- num / den
  fill <- matrix(pooled, nrow(th), ncol(th), byrow = TRUE)
  th[den == 0] <- fill[den == 0]
  t(th)
}


# =============================================================================
#  Sparse weighted least squares at a fixed support size
# =============================================================================

## Sparse design for patterns D restricted to the `allowed` (K x M) pairs.
design_template <- function(dat, D, allowed) {
  G <- nrow(D); M <- ncol(D)
  col_id <- matrix(0L, dat$K, M); col_id[allowed] <- seq_len(sum(allowed))
  cols <- which(allowed, arr.ind = TRUE)
  obs <- which(dat$R, arr.ind = TRUE); nobs <- nrow(obs); lk <- obs[, 1]; kk <- obs[, 2]
  sd_obs <- sqrt(dat$v[obs])
  ii <- jj <- aa <- kx <- bb <- list()
  for (g in seq_len(G)) for (a in which(D[g, ] == 1)) {
    ok <- allowed[kk, a]
    if (!any(ok)) next
    ii[[length(ii) + 1]] <- (g - 1) * nobs + which(ok)
    jj[[length(jj) + 1]] <- col_id[cbind(kk[ok], a)]
    aa[[length(aa) + 1]] <- rep(a, sum(ok))
    kx[[length(kx) + 1]] <- kk[ok]
    bb[[length(bb) + 1]] <- 1 / sd_obs[ok]
  }
  list(i = unlist(ii), j = unlist(jj), a = unlist(aa), k = unlist(kx), base = unlist(bb),
       y = rep(dat$beta[obs] / sd_obs, G), row_l = rep(lk, G),
       row_g = rep(seq_len(G), each = nobs),
       nrow = G * nobs, ncol = sum(allowed), colk = cols[, 1], cola = cols[, 2],
       K = dat$K, M = M)
}

abess_fit <- function(X, y, w, s) {
  args <- list(y = y, weight = w, family = "gaussian", tune.type = "bic",
               tune.path = "sequence", support.size = s, normalize = 0,
               fit.intercept = FALSE)
  fit <- tryCatch(suppressMessages(do.call(abess::abess, c(list(x = X), args))),
                  error = function(e)
                    suppressMessages(do.call(abess::abess, c(list(x = as.matrix(X)), args))))
  as.numeric(fit$beta[, 1])
}

## Calibrated sparse fit at support size s.  B is returned on the ORIGINAL
## scale, so theta = B %*% t(D) are the fitted cluster effects.
sparse_fit <- function(tpl, W, Cmat, s) {
  B <- matrix(0, tpl$K, tpl$M)
  s <- min(s, tpl$ncol)
  if (s <= 0) return(B)
  cc <- Cmat[cbind(tpl$k, tpl$a)]
  x <- ifelse(cc > 1e-8, tpl$base / (cc + 1e-5), 0)
  X <- Matrix::sparseMatrix(i = tpl$i, j = tpl$j, x = x, dims = c(tpl$nrow, tpl$ncol))
  bt <- abess_fit(X, tpl$y, W[cbind(tpl$row_l, tpl$row_g)], s)
  ccol <- Cmat[cbind(tpl$colk, tpl$cola)]
  B[cbind(tpl$colk, tpl$cola)] <- ifelse(ccol > 1e-8, bt / (ccol + 1e-5), 0)
  B
}


# =============================================================================
#  Algorithm 1: one start, Stage 1, fixed G and s
# =============================================================================

em_fit <- function(dat, D, tpl, G, s, Z0 = NULL, Psi0 = NULL, tol, NMAX, ic) {
  if (is.null(Psi0)) {
    pr <- colMeans(Z0)
    B  <- init_theta(dat, Z0) %*% MASS::ginv(t(D))
  } else {
    B <- Psi0$B; pr <- Psi0$pr
  }
  es <- estep(dat, B %*% t(D), pr); Q <- es$Q
  h <- 0; conv <- FALSE
  repeat {
    W  <- es$W
    pr <- colMeans(W)
    B  <- sparse_fit(tpl, W, calib(dat, W, D), s)
    es <- estep(dat, B %*% t(D), pr)
    delta <- abs(es$Q - Q) / (1 + abs(Q)); Q <- es$Q; h <- h + 1
    if (delta <= tol) { conv <- TRUE; break }
    if (h >= NMAX) break
  }
  W <- es$W                       # responsibilities at the final parameters
  theta <- B %*% t(D)
  loss <- wrss(dat, W, theta)
  icv <- ic_pen(sum(B != 0) + (G - 1), tpl$ncol, dat$n0, ic)
  list(B = B, theta = theta, pr = pr, W = W, Z = hard(W), q_loss = loss, ic = icv,
       gic = loss + icv, Qobs = Q, converged = conv, iter = h, s = s)
}


# =============================================================================
#  Algorithm 2: cached discrete golden-section search over s
# =============================================================================

support_search <- function(Ffun, smax) {
  phi <- (sqrt(5) - 1) / 2; a <- 0; b <- max(0, round(smax)); cache <- list()
  ev <- function(s) {
    s <- as.integer(min(max(round(s), a), b)); key <- as.character(s)
    if (is.null(cache[[key]])) cache[[key]] <<- Ffun(s)
    s
  }
  gic <- function(s) cache[[as.character(s)]]$gic
  cc <- ev(b - phi * (b - a)); dd <- ev(a + phi * (b - a)); ev(a); ev(b)
  while (b - a > 4) {
    if (cc >= dd) break                       # rounding has collapsed the interior
    if (gic(cc) <= gic(dd)) { b <- dd; dd <- cc; cc <- ev(b - phi * (b - a)) }
    else                    { a <- cc; cc <- dd; dd <- ev(a + phi * (b - a)) }
  }
  for (s in a:b) ev(s)                        # exhaust the terminal bracket
  keys <- as.integer(names(cache)); g <- vapply(cache, function(f) f$gic, numeric(1))
  best <- keys[order(g, keys)[1]]             # ties -> smaller support
  list(s = best, fit = cache[[as.character(best)]], cache = cache,
       path = data.frame(s = keys, gic = unname(g))[order(keys), ])
}


# =============================================================================
#  Algorithm 3, Stage 1 at one G: common starts, warm starts, best start kept
# =============================================================================

stage1_at_G <- function(dat, G, nstart, nrefine, tol, NMAX, ic, seed, verbose) {
  D <- stage1_basis(G)
  tpl <- design_template(dat, D, matrix(TRUE, dat$K, ncol(D)))
  set.seed(seed + G)
  starts <- replicate(nstart, rand_partition(dat$L, G), simplify = FALSE)
  H <- replicate(nstart, list(), simplify = FALSE)   # per-start fits, keyed by s

  F1 <- function(s) {
    res <- lapply(seq_len(nstart), function(r) {
      prev <- H[[r]]
      if (!length(prev)) {
        em_fit(dat, D, tpl, G, s, Z0 = starts[[r]], tol = tol, NMAX = NMAX, ic = ic)
      } else {
        sk <- as.integer(names(prev)); near <- sk[order(abs(sk - s), sk)[1]]
        em_fit(dat, D, tpl, G, s, Psi0 = prev[[as.character(near)]],
               tol = tol, NMAX = NMAX, ic = ic)
      }
    })
    for (r in seq_len(nstart))
      H[[r]][[as.character(s)]] <<- list(B = res[[r]]$B, pr = res[[r]]$pr)
    gics <- vapply(res, function(f) f$gic, numeric(1))
    conv <- vapply(res, function(f) f$converged, logical(1))
    cand <- if (any(conv)) which(conv) else seq_len(nstart)
    rstar <- cand[which.min(gics[cand])]
    best <- res[[rstar]]
    ## restart diagnostics at this (G, s): what each start reached.  Reported
    ## only; the partition is the best converged start, not a consensus.
    best$cluster_list <- lapply(res, function(f) setNames(max.col(f$W, "first"), dat$context))
    best$cluster_gic  <- gics
    best$cluster_converged <- conv
    best$selected_start <- rstar
    if (verbose) message(sprintf("    G = %d, s = %4d: GIC %.5f (start %d), %d/%d converged",
                                 G, s, best$gic, rstar, sum(conv), nstart))
    best
  }

  smax <- if (G == 1) round(dat$K / 2) else dat$K
  if (verbose) message(sprintf("++ Stage 1, G = %d: %d random starts, support 0..%d ++",
                               G, nstart, smax))
  srch <- support_search(F1, smax)
  fit <- srch$fit
  if (nrefine > 0) {
    for (r in seq_len(nrefine)) {
      f <- em_fit(dat, D, tpl, G, srch$s, Z0 = rand_partition(dat$L, G),
                  tol = tol, NMAX = NMAX, ic = ic)
      if (f$converged && f$gic < fit$gic) {
        f$cluster_list <- fit$cluster_list; f$cluster_gic <- fit$cluster_gic
        f$cluster_converged <- fit$cluster_converged; f$selected_start <- NA_integer_
        fit <- f
      }
    }
  }
  list(G = G, s = srch$s, fit = fit, path = srch$path, D = D)
}


# =============================================================================
#  Algorithm 3: Stage 2 and the unpenalized refit at the selected G
# =============================================================================

finish_at_G <- function(dat, st1, G, ic, tau, verbose) {
  f1 <- st1$fit; Z <- f1$Z
  E <- crossprod(dat$Rn, Z) > 0                 # K x G: the clusters E_k
  st2_path <- NULL; s2 <- st1$s; B <- f1$B; D <- st1$D

  if (G > 2) {
    D <- all_subsets(G)
    allowed <- ((1 - E) %*% D) == 0              # A^(2)_k: patterns inside E_k
    tpl <- design_template(dat, D, allowed)
    Cmat <- calib(dat, Z, D)
    F2 <- function(s) {
      B2 <- sparse_fit(tpl, Z, Cmat, s)
      loss <- wrss(dat, Z, B2 %*% t(D))
      list(B = B2, q_loss = loss,
           gic = loss + ic_pen(sum(B2 != 0), tpl$ncol, dat$n0, ic))
    }
    if (verbose) message(sprintf("++ Stage 2: %d allowed feature-pattern pairs ++", sum(allowed)))
    srch <- support_search(F2, st1$s)
    B <- srch$fit$B; st2_path <- srch$path; s2 <- srch$s
  }

  ## Unpenalized post-selection refit on the estimable clusters of each feature.
  A <- crossprod(dat$P, Z); Cc <- crossprod(dat$Yb, Z)
  theta <- matrix(NA_real_, dat$K, G, dimnames = list(dat$feature, as.character(seq_len(G))))
  for (k in seq_len(dat$K)) {
    Ek <- which(E[k, ]); Mk <- which(abs(B[k, ]) > tau)
    if (!length(Ek)) next
    if (!length(Mk)) { theta[k, Ek] <- 0; next }
    Ds <- D[Ek, Mk, drop = FALSE]
    theta[k, Ek] <- Ds %*% (MASS::ginv(t(Ds) %*% diag(A[k, Ek], nrow = length(Ek)) %*% Ds) %*%
                              t(Ds) %*% Cc[k, Ek])
  }
  list(theta = theta, B = B, D = D, s2 = s2, st2_path = st2_path, E = E)
}


# =============================================================================
#  Entry point
# =============================================================================

#' SMESH from context-level summary statistics.
#'
#' @param summary.stats named list, one element per context, each with `est` and
#'   `stderr` (feature x covariate, NA where not estimable) and `n`.
#' @param G fixed number of clusters, or NULL to select it by the global minimum
#'   of the Stage 1 GIC over 1..G.max.
#' @param tune.type complexity penalty: "BIC" (default), "HBIC", "KBIC", "EBIC".
#' @param tol,NMAX EM tolerance on the relative change in the observed-data
#'   negative log-likelihood, and the maximum number of EM iterations.
#' @param nperm number of random starts per candidate G.
#' @param nrefine extra random starts fitted at the selected support size only.
#' @param G.max largest candidate G when `G` is NULL; default min(L-1, [2L/3]).
#' @param seed the starts for candidate G are drawn after set.seed(seed + G).
#' @param doc kept for backward compatibility; unused.
#' @param verbose progress messages.
smesh.meta.summary <- function(summary.stats,
                               G = NULL,
                               tune.type = c("BIC", "HBIC", "KBIC", "EBIC"),
                               tol = 1e-5,
                               NMAX = 30,
                               nperm = 5,
                               nrefine = 0,
                               G.max = NULL,
                               seed = 2026,
                               doc = "./",
                               verbose = FALSE) {

  tune.type <- match.arg(tune.type)
  TAU <- 1e-8

  cov.int.ID <- sort(unique(unlist(lapply(summary.stats, function(d) colnames(d$est)))))
  output.result <- vector("list", length(cov.int.ID))
  names(output.result) <- cov.int.ID

  for (cov.name in cov.int.ID) {
    if (verbose) message("++ Covariate of interest: ", cov.name, " ++")
    dat <- prep_data(summary.stats, cov.name)

    gmax <- if (!is.null(G.max)) G.max else min(dat$L - 1, round(2 * dat$L / 3))
    cand <- if (!is.null(G)) G else seq_len(max(1, gmax))
    if (any(cand > dat$L))
      stop("every candidate G must be at most the number of contexts (", dat$L, ").")

    st1 <- lapply(cand, function(g)
      stage1_at_G(dat, g, nperm, nrefine, tol, NMAX, tune.type, seed, verbose))
    names(st1) <- as.character(cand)
    gics <- vapply(st1, function(x) x$fit$gic, numeric(1))
    Ghat <- cand[order(gics, cand)[1]]           # ties -> smaller G
    if (verbose && length(cand) > 1)
      message(sprintf("++ G = %d selected (global minimum of the Stage 1 GIC) ++", Ghat))

    sel <- st1[[as.character(Ghat)]]
    fin <- finish_at_G(dat, sel, Ghat, tune.type, TAU, verbose)
    f1  <- sel$fit

    ## ---- output, in the legacy field names ---------------------------------
    mu <- fin$theta
    mu[is.na(mu)] <- 0                           # 0 where a cluster has no data
    W      <- matrix(f1$Z, dat$L, Ghat, dimnames = list(dat$context, as.character(seq_len(Ghat))))
    W_post <- matrix(f1$W, dat$L, Ghat, dimnames = list(dat$context, as.character(seq_len(Ghat))))
    mu_stage1 <- f1$theta
    dimnames(mu_stage1) <- list(dat$feature, as.character(seq_len(Ghat)))
    rownames(fin$B) <- dat$feature

    output.result[[cov.name]] <- list(
      mu        = mu,
      W         = W,
      Pi        = setNames(f1$pr, as.character(seq_len(Ghat))),
      gic       = f1$gic,
      ic        = f1$ic,
      q_loss    = f1$q_loss,
      B         = fin$B,
      s.lambda  = fin$s2,
      mu_refit  = fin$theta,                     # NA where a cluster has no data
      mu_stage1 = mu_stage1,
      W_post    = W_post,
      cluster_list = f1$cluster_list,
      cluster_gic  = f1$cluster_gic,
      gic_path  = list(
        by_G   = data.frame(G = cand, s = vapply(st1, function(x) x$s, numeric(1)),
                            gic = unname(gics), selected = cand == Ghat, row.names = NULL),
        stage1 = sel$path,
        stage2 = fin$st2_path))
  }

  output.result
}
