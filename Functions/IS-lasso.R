
# ----------  helpers  -------------------------------------------------

## 1. Third-largest |tᵀx| per factor * 2   (upper bound for your scaling)
lambda_upper_bound <- function(X, R) {
  N  <- nrow(X)
  sv <- svd(X, nu = R, nv = 0)          # first R PC score vectors
  T0 <- sqrt(N) * sv$u                  # ‖t_r‖₂ = √N
  C  <- abs(crossprod(T0, X))           # R × J
  2 * min(apply(C, 1L, function(v) sort(v, TRUE)[3L]))
}

## 2. Count non-zero loadings in every factor (column of P)
enough_items <- function(P, min_items = 3L, eps = 1e-8)
  colSums(abs(P) > eps) >= min_items


# ----------  wrapper around *your* multistart LASSO estimator ---------

fit_once <- function(DATA, R, lambda,
                     MaxIter = 250, eps = 1e-5, nstarts = 1)
{
  res <- MULTISTART_rESEMl1(DATA, R       = R,
                            MaxIter = MaxIter,
                            eps     = eps,
                            nstarts = nstarts,
                            lambda  = lambda)
  res$loadings    # only P is needed for the λ-search
}



#' Title Find the maximum value of lambda so that each factor has at least 3 items
#'
#' @param DATA Standardized data of items X
#' @param R Number of factors
#' @param MaxIter 
#' @param eps 
#' @param nstarts 
#' @param min_items Minimum number of items per factor
#' @param tol 
#' @param max_iter 
#' @param verbose 
#'
#' @returns
#' @export
#'
#' @examples
find_lambda_max <- function(DATA, R,
                            MaxIter   = 250,
                            eps       = 1e-5,
                            nstarts   = 1,
                            min_items = 3L,
                            tol       = 1e-4,
                            max_iter  = 25,
                            verbose   = TRUE)
{
  ## 0. rough upper bound from SVD a-la PCA
  lam_hi <- lambda_upper_bound(DATA, R)
  if (verbose) cat("initial upper bound  =", signif(lam_hi, 4), "\n")
  
  ## 1. make sure that bound is *really* infeasible for your algorithm
  repeat {
    P <- try(fit_once(DATA, R, lam_hi,
                      MaxIter = MaxIter, eps = eps, nstarts = nstarts),
             silent = TRUE)
    ok <- !(inherits(P, "try-error") ||
              any(!enough_items(P, min_items)))
    if (!ok) break               # already infeasible  →  continue below
    lam_hi <- lam_hi * 2         # still feasible     →  go higher
    if (verbose) cat("doubling : λ =", signif(lam_hi, 4), "(still feasible)\n")
  }
  
  lam_lo <- 0                    # 0 is always feasible
  
  ## 2. binary search
  for (k in seq_len(max_iter)) {
    lam_mid <- 0.5 * (lam_hi + lam_lo)
    
    P <- try(fit_once(DATA, R, lam_mid,
                      MaxIter = MaxIter, eps = eps, nstarts = nstarts),
             silent = TRUE)
    ok <- !(inherits(P, "try-error") ||
              any(!enough_items(P, min_items)))
    
    if (verbose)
      cat(sprintf("%2d) λ = %.6g  %s\n",
                  k, lam_mid, if (ok) "✓ feasible" else "✗ too few"))
    
    if (ok)      lam_lo <- lam_mid  else  lam_hi <- lam_mid
    if ((lam_hi - lam_lo) < tol * (1 + lam_lo)) break
  }
  
  lam_lo              # ← the largest λ that still has ≥ 3 items / factor
}

#' Title Model selection strategy using Index of Sparseness 
#'
#' @param X Standardized data matrix of items
#' @param R Number of factors
#' @param card Number of nonzero loadings
#' @param MaxIter 
#' @param eps 
#' @param nstarts Number of starts
#'
#' @returns
#' @export
#'
#' @examples
IS_strategy <- function(X, R, MaxIter, eps, nstarts, lambda_max) {
  J <- dim(X)[2]
  I <- dim(X)[1]
  
  lasso_vec <- seq(0.1*I, lambda_max, length.out = 50)
  K <- length(lasso_vec)
  
  avec <- matrix(nrow = K, ncol = 6)
  colnames(avec) <- c( "lambda", "IS", "PEV", "Prop0", "MinNonZeroL", "MaxSDL")
  for (k in 1:K) {
    # retry loop
    success <- FALSE
    while (!success) {
      res <- tryCatch(
        IS_rESEMl1(
          DATA    = X,
          R       = 3,
          lambda  = lasso_vec[k],
          MaxIter = 20,
          eps     = 1e-4,
          nstarts = nstarts
        ),
        error = function(e) {
          message(sprintf("k = %d: IS_rESEMl1 failed (lambda = %g), retrying...", 
                          k, lasso_vec[k]))
          NULL
        }
      )
      if (!is.null(res)) {
        a       <- res
        success <- TRUE
      }
      # else it loops and tries again
    }
    avec[k, ] <- c(
      lasso_vec[k], a$value, a$vaf, a$propzero, 
      a$smallestP, a$maxsdP
    )
  }
  round(avec, 3)
  
  selcardinality <- lasso_vec[which.max(avec[,2])]
  
  return(selcardinality)
}
