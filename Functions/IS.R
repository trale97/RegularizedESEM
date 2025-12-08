
#' Title Model selection strategy using Index of Sparseness 
#'
#' @param X Standardized data matrix of items
#' @param R Number of factors
#' @param card Number of nonzero loadings
#' @param MaxIter 
#' @param eps 
#' @param nstarts Number of starts
#' @param TOL 
#' @param THR 
#'
#' @returns
#' @export
#'
#' @examples
IS_strategy <- function(X, R, MaxIter, eps, nstarts, TOL, THR) {
  J <- dim(X)[2]
  
  cardvec <- 3:J
  K <- length(cardvec)
  avec <- matrix(nrow = length(cardvec), ncol = 6)
  colnames(avec) <- c( "K", "IS", "PEV", "Prop0", "MinNonZeroL", "MaxSDL")
  for (k in 1:K) {
    a <- IS_rESEM(X = scale(X), R = R, card = rep(cardvec[k], R), MaxIter = 20, eps = 10^-4, nstarts = nstarts)
    avec[k, ] <- c(
      cardvec[k], a$value, a$vaf, a$propzero, 
      a$smallestP, a$maxsdP
    )
  }
  round(avec, 3)
  
  # second check solutions with cardinality at level of the loading matrix
  # second check solutions with cardinality at level of the loading matrix
  cardvec <- (3 * R):(J * R)
  L <- length(cardvec)
  avec <- rbind(avec, matrix(nrow = length(cardvec), ncol = 6))
  
  for (l in 1:L) {
    a <- IS_rESEM(X = scale(X), R = R, card = cardvec[l], MaxIter = 20, eps = 10^-4, nstarts = nstarts)
    avec[K+l, ] <- c(
      cardvec[l], a$value, a$vaf, a$propzero, 
      a$smallestP, a$maxsdP
    )
  }
  round(avec, 4)
  
  placeholder <- seq(1:(K + L))
  indexnonzeroL <- placeholder[avec[, "MinNonZeroL"] < THR][1] - 1 # THR holding selection
  if (indexnonzeroL < K) {
    placeholder2 <- placeholder[(K + 1):(K + L)]
    indexnonzeroL2 <- placeholder2[avec[(K + 1):(K + L), "MinNonZeroL"] < THR][1] - 1 # THR holding selection
    if (avec[indexnonzeroL2, 1] < avec[indexnonzeroL, 1] * R) {
      indexnonzeroL <- indexnonzeroL2
    } else if (avec[indexnonzeroL2, 1] == avec[indexnonzeroL, 1] * R) {
      if (round(avec[indexnonzeroL2, 3], THR) > round(avec[indexnonzeroL, 3], THR)) {
        indexnonzeroL <- indexnonzeroL2
      }
    }
  }
  selcardinality <- avec[indexnonzeroL, 1]
  if (indexnonzeroL < K + 1) {
    selcardinality <- c(rep(avec[indexnonzeroL, 1], R))
  }
  
  #TOL = 2
  samepve <- round(avec[, 3], TOL) == round(avec[indexnonzeroL, 3], TOL)
  
  index <- indexnonzeroL
  if (is.na(indexnonzeroL) || max(round(avec[samepve, 2], 3)) > round(avec[indexnonzeroL, 2], 3)) {
    index <- placeholder[avec[, 2] == max(avec[samepve, 2])][1]
    selcardinality <- avec[index, 1]
    if (index < K) {
      selcardinality <- rep(avec[index, 1], R)
    }
  }
  
  return(selcardinality)
}