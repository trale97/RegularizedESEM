#### Title: Mplus
#### Author: Tra Le
#### Created: April 10, 2023
#### Last modified: Jan 24, 2025

#########################################################################################
######                         Mplus simulation analysis                         #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

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
Simulation = foreach(i = 1:Ndatasets,
                     .packages = c("MplusAutomation"),
                     .combine = rbind)%dopar%{
                       
                       source('functions/num_correct.R', local = T)
                       
                       # Load data
                       load(paste0("DATA/data", i, ".RData"))
                       
                       X <- scale(out$X)
                       y <- out$y
                       data <- as.data.frame(cbind(X,y))
                       
                       
                       
                       ## Measurement model
                       # population values
                       Ptrue <- out$Ptrue
                       nfactors <- dim(Ptrue)[2]
                       J <- dim(Ptrue)[1]
                       
                       colnames(data) <- c(paste0("v", 1:J), "y")
                       
                       # syntax for Mplus
                       if (nfactors == 3){
                         if (J == 9){
                           varstr <- paste("usevariables = v1 - v9 y;")
                           modelstr <- paste("f1 BY v1-v9*1 (a1-a9);
                                              f2 BY v1-v9*1 (a10-a18);
                                              f3 BY v1-v9*1 (a19-a27);
                                              f1-f3@1;
                                              y ON f1 f2 f3")
                           #prior <- "a1-a27~Geomin(3,1);"
                         } else {
                           varstr <- paste("usevariables = v1 - v15 y;")
                           modelstr <- paste("f1 BY v1-v15*1 (a1-a15);
                                              f2 BY v1-v15*1 (a16-a30);
                                              f3 BY v1-v15*1 (a31-a45);
                                              f1-f3@1;
                                              y ON f1 f2 f3")
                           #prior <- "a1-a45~Geomin(3,1);"
                         }
                       } else {
                         if (J == 15){
                           varstr <- paste("usevariables = v1 - v15 y;")
                           modelstr <- paste("f1 BY v1-v15*1 (a1-a15);
                                              f2 BY v1-v15*1 (a16-a30);
                                              f3 BY v1-v15*1 (a31-a45);
                                              f4 BY v1-v15*1 (a46-a60);
                                              f5 BY v1-v15*1 (a61-a75);
                                              f1-f5@1;
                                              y ON f1 f2 f3 f4 f5")
                           #prior <- "a1-a75~Geomin(5,1);"
                         } else {
                           varstr <- paste("usevariables = v1 - v25 y;")
                           modelstr <- paste("f1 BY v1-v25*1 (a1-a25);
                                              f2 BY v1-v25*1 (a26-a50);
                                              f3 BY v1-v25*1 (a51-a75);
                                              f4 BY v1-v25*1 (a76-a100);
                                              f5 BY v1-v25*1 (a101-a125);
                                              f1-f5@1;
                                              y ON f1 f2 f3 f4 f5")
                           #prior <- "a1-a125~Geomin(5,1);"
                         }
                       }
                       
                       v_values <- c(0.1, 1, 10, 100)
                       success <- FALSE  # flag for success
                       
                       
                       for (v in v_values) {
                         
                         # Define prior string dynamically
                         if (nfactors == 3) {
                           if (J == 9)  prior <- paste0("a1-a27~Geomin(3,", v, ");")
                           else         prior <- paste0("a1-a45~Geomin(3,", v, ");")
                         } else {
                           if (J == 15) prior <- paste0("a1-a75~Geomin(5,", v, ");")
                           else         prior <- paste0("a1-a125~Geomin(5,", v, ");")
                         }
                         
                         # Define Mplus model object
                         m.penefa <- mplusObject(
                           TITLE = paste0("Penalized ESEM with Regression (v=", v, ");"),
                           VARIABLE = varstr,
                           MODEL = modelstr,
                           analysis = "iter=10000;",
                           OUTPUT = "standardized;",
                           MODELPRIORS = prior,
                           rdata = data
                         )
                         
                         # Initialize matrices
                         Pmatrix <- matrix(NA, nrow = J, ncol = nfactors)
                         pvalue <- Pmatrix
                         
                         # Track warnings
                         warned <- FALSE
                         dat_file <- paste0("mplus_files/ESEM_data",  i, ".dat")
                         inp_file <- paste0("mplus_files/ESEM_model", i, "_v", v, ".inp")
                         
                         # Try fitting model
                         m.penefa.fit <- withCallingHandlers(
                           {
                             mplusModeler(
                               m.penefa,
                               dataout   = dat_file,
                               modelout  = inp_file,
                               run       = 1L,
                               writeData = "ifmissing"
                             )
                           },
                           warning = function(w) {
                             warned <<- TRUE
                             invokeRestart("muffleWarning")
                           }
                         )
                         
                         # Extract loadings only if successful
                         if (!warned && !is.null(m.penefa.fit)) {
                           loadingsparams <- m.penefa.fit$results$parameters$unstandardized
                           for (r in seq_len(nfactors)) {
                             strR <- paste0("F", r, ".BY")
                             idx  <- loadingsparams$paramHeader == strR
                             Pmatrix[, r] <- loadingsparams$est[idx]
                             pvalue[, r]  <- loadingsparams$pval[idx]
                           }
                         }
                         
                         # Check if the model worked
                         if (all(is.na(Pmatrix)) == FALSE) {
                           success <- TRUE
                           message("Model succeeded with v = ", v)
                           break
                         } else {
                           message("Model failed with v = ", v, " — trying next value...")
                         }
                       }
                       
                       if (all(is.na(Pmatrix))) {
                         converge <- 0
                         Pmatrix_sim <- Pmatrix
                         corrate_final = corrate_sim_final = NA
                         beta = NA
                       }else{
                         converge <- 1
                         alpha <- 0.05
                         Pmatrix_sim <- Pmatrix
                         Pmatrix_sim[pvalue >= alpha] <- 0
                         # zero/non-zero structure recovery rate + Tucker's congruence
                         perm <- gtools::permutations(nfactors, nfactors)
                         corrate <- vector(length = nrow(perm))
                         corrate_sim <- vector(length = nrow(perm))
                         
                         for (p in 1:nrow(perm)) {
                           corrate[p] <- num_correct(Ptrue, round(Pmatrix[,perm[p,]],nfactors))
                           corrate_sim[p] <- num_correct(Ptrue, round(Pmatrix_sim[,perm[p,]], nfactors))
                         }
                         corrate_final <- max(corrate)
                         corrate_sim_final <- max(corrate_sim)
                         
                         index <- which.max(corrate_sim)
                         
                         Pmatrix <- Pmatrix[,perm[index,]]
                         Pmatrix_sim <- Pmatrix_sim[,perm[index,]]
                         
                         corsign <- sign(diag(cor(Ptrue, Pmatrix)))
                         
                         Pmatrix <- Pmatrix%*%diag(corsign)
                         Pmatrix_sim <- Pmatrix_sim%*%diag(corsign)
                         
                         ## Structural model
                         ind_beta <- loadingsparams$paramHeader == "Y.ON"
                         beta <- loadingsparams$est[ind_beta]
                       }
                       
                       
                       
                       ## measure summary
                       measure_summary <- matrix(c(corrate_final, corrate_sim_final, converge), nrow = 1, byrow = T)
                       colnames(measure_summary) <- c("zero_rec", "zero_rec_sim", "converge")
                       # Saving performance measures
                       output = list(Pmatrix = Pmatrix, 
                                     Pmatrix_sim = Pmatrix_sim,
                                     beta = beta, 
                                     m.penefa.fit,
                                     measure_summary = measure_summary)
                       
                       save(output, file = paste0("mplus_nostd_new/Mplus", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)
