#### Title: Regularized ESEM
#### Author: Tra Le
#### Created: April 10, 2023
#### Last modified: Jan 24, 2025

#########################################################################################
######                         rESEM-CC simulation analysis                         #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

# required package: lslx

# load simulation setup
load("DATA/Info_simulation.RData")
Info_matrix = Infor_simulation$design_matrix_replication
Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication

# register cores
library(doParallel)
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation = foreach(id = 1:length(ind),
                     .combine = rbind)%dopar%{
                       
                       source("functions/CCrESEM.R")
                       source("functions/IS.R")
                       source("functions/num_correct.R")
                       
                       # Load data
                       #i <- ind[id]
                       load(paste0("DATA/data", i, ".RData"))
                       
                       X <- out$X # not standardized
                       y <- out$y # not standardized
                       N <- dim(X)[1]
                       
                       # Population models
                       Ptrue = out$Ptrue
                       Htrue = out$Htrue
                       nfactors <- dim(Ptrue)[2]
                       
                       selcard <- IS_strategy(X = scale(X), R = nfactors, MaxIter, eps, nstarts = 100, TOL = 2, THR = .2)
                       
                       # run model with the selected cardinality
                       model <- MULTISTART_CCrESEM(DATA = scale(X),
                                                   R = nfactors,
                                                   CARD = selcard,
                                                   MaxIter = 50,
                                                   eps = 10^-4,
                                                   nstarts = 100)
                       Pmatrix <- model$loadings
                       Hmatrix <- model$scores
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(nfactors, nfactors)
                       corrate <- vector(length = nrow(perm))
                       
                       for (p in 1:nrow(perm)) {
                         corrate[p] <-num_correct(Ptrue, round(Pmatrix[,perm[p,]],3))
                         
                       }
                       corrate_final <- max(corrate)
                       
                       
                       index <- which.max(corrate)
                       
                       Pmatrix <- Pmatrix[,perm[index,]]
                       corsign <- sign(diag(cor(Ptrue, Pmatrix)))
                       
                       Pmatrix <- Pmatrix%*%diag(corsign)
                       
                       # structural model
                       Hmatrix <- Hmatrix[,perm[index,]]%*%diag(corsign)
                       tucker <- sum(diag(abs(psych::factor.congruence(Htrue,Hmatrix))))/nfactors
                       beta_est <- lm(y ~ 0 + Hmatrix)
                       
                       # convergence
                       converge <- model$converge
                       
                       ######### KNOWN CARDINALITY ##################
                       # run model with the known cardinality
                       model_known <- MULTISTART_CCrESEM(DATA = scale(X),
                                                         R = nfactors,
                                                         CARD = sum(Ptrue != 0),
                                                         MaxIter = 50,
                                                         eps = 10^-4,
                                                         nstarts = 100)
                       Pmatrix_known <- model_known$loadings
                       Hmatrix_known <- model_known$scores
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(nfactors, nfactors)
                       corrate <- vector(length = nrow(perm))
                       
                       for (p in 1:nrow(perm)) {
                         corrate[p] <-num_correct(Ptrue, round(Pmatrix_known[,perm[p,]],3))
                         
                       }
                       corrate_known <- max(corrate)
                       
                       
                       index <- which.max(corrate)
                       
                       Pmatrix_known <- Pmatrix_known[,perm[index,]]
                       corsign <- sign(diag(cor(Ptrue, Pmatrix_known)))
                       
                       Pmatrix_known <- Pmatrix_known%*%diag(corsign)
                       
                       # structural model
                       Hmatrix_known <- Hmatrix_known[,perm[index,]]%*%diag(corsign)
                       tucker_known <- sum(diag(abs(psych::factor.congruence(Htrue,Hmatrix_known))))/nfactors
                       beta_known <- lm(y ~ 0 + Hmatrix_known)
                       
                       # convergence
                       converge_known <- model_known$converge
                       
                       ## measure summary
                       measure_summary <- matrix(c(corrate_final, tucker, converge, corrate_known, tucker_known, converge_known), nrow = 1, byrow = T)
                       colnames(measure_summary) <- c("zero_rec_IS", "tucker_scores_IS", "converge", "zero_rec_known", "tucker_scores_known", "convergence_known")
                       # Saving performance measures
                       output = list(measure_summary = measure_summary, 
                                     Pmatrix = Pmatrix, 
                                     Hmatrix = Hmatrix,
                                     Pmatrix_known = Pmatrix_known,
                                     Hmatrix_known = Hmatrix_known,
                                     beta_est = beta_est,
                                     beta_known = beta_known)
                       
                       save(output, file = paste0("rESEM/rESEM", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)