findLasso <- function(dat, Ptrue, maxItr, lassou){
  
  percentageZeroesInData <- sum(Ptrue == 0)
  percentageInW <- 0
  i  <- 0
  lassol <- 0 
  converged <- FALSE
  conv0 <- 0
  lasso1 <- 1
  while(conv0 == 0){
    
    lasso <- (lassol + lassou) / 2
    fit <- MULTISTART_rESEMl1(dat, lambda = lasso, R = 3, MaxIter = 200, eps = 10^-6, nstarts = 20)
    percentageInW <- sum(round(fit$loadings,3) == 0)
    if( percentageZeroesInData > percentageInW){
      lassol  <- lasso
      
    } else {
      lassou  <- lasso
    }
    #print(lasso)
    
    if (abs(percentageZeroesInData - percentageInW) < 0){
      conv0 <- 1
    }
    
    
    if(i > maxItr){
      conv0 <- 1
    }
    
    
    if (abs(lasso1 - lasso) < .01){
      conv0 <- 1
    }
    lasso1 <- lasso
    i <- i + 1
  }
  
  if( i < maxItr ){
    converged <- TRUE
  } 
  
  return(list(lasso = lasso, converged = converged))
}