rm(list = ls())
library(dplyr)
###### Data generation only loadings + no factor scores #####

## Setup
#' Sample size: 50, 100, 500
#' Number of factors: 3, 5
#' Numer of items per factor: 3, 5
#' Correlation between exogenous factors: 0, .3
#' Reliablity: low, high
#' Cross-loadings: moderate


## Simulation setup 
PL = sqrt(.6)
nvar.factor = c(3,5)
nfactors = c(3,5)
N = c(50, 100, 500)
reliability = c(.3, .8)
beta = .1
factor_corr = c(0, .3)
VAFy = c(.5, .9)
CL = c("yes", "no")
n_replications = 50

design_matrix <- expand.grid(nvar.factor = nvar.factor, nfactors = nfactors,
                             sample = N, reliability = reliability,
                             factor_corr = factor_corr, VAFy = VAFy, CL = CL)
design_matrix_replication <- design_matrix %>% slice(rep(1:n(), each = 50))

Infor_simulation = list(n_data_sets = nrow(design_matrix_replication), n_replications = n_replications,
                        design_matrix_replication = design_matrix_replication)

save(Infor_simulation, file = "DATA/Info_simulation.RData")

#### setting the number of cores
rm(list = ls())
library(doParallel) 
no_cores <- detectCores() 
c1 <- makePSOCKcluster(floor( no_cores*.8))
registerDoParallel(c1)

#### start generating data
load("DATA/Info_simulation.RData")
design_matrix_replication <- Infor_simulation$design_matrix_replication

results_sim1_data1 <- foreach(i=1:nrow(design_matrix_replication),
                              .options.RNG = 1997,
                              .packages = c("MASS", "lavaan"),
                              .combine=rbind)%dopar%{
                                source('Functions/generation_function.R', local = T)
                                
                                # set the specific values of the parameters
                                PL = sqrt(.6)
                                nvar.factor = design_matrix_replication$nvar.factor[i]
                                nfactors = design_matrix_replication$nfactors[i]
                                factor_corr = design_matrix_replication$factor_corr[i]
                                beta = rep(.1, nfactors)
                                reliability = design_matrix_replication$reliability[i]
                                N = design_matrix_replication$sample[i]
                                VAFy = design_matrix_replication$VAFy[i]
                                cross_loadings = design_matrix_replication$CL[i]
                                if (cross_loadings == "yes"){
                                  CL = TRUE
                                } else {
                                  CL = FALSE
                                }
                                
                                output = data_generation(PL = PL,
                                                         nvar.factor = nvar.factor,
                                                         nfactors = nfactors,
                                                         reliability = reliability,
                                                         factor_corr = factor_corr,
                                                         N = N,
                                                         beta = beta,
                                                         VAFy = VAFy,
                                                         CL = CL)
                                
                                X <- output$X
                                y <- output$y
                                Ptrue <- output$Ptrue
                                Htrue <- output$Htrue
                                
                                # Saving data
                                out = list(X = X, y = y, Ptrue = Ptrue, Htrue =Htrue, teller = i)
                                save(out, file=paste0("DATA/data",i, ".RData"))
                                
                              }

stopCluster(c1)


