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
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation = foreach(i = 1:Ndatasets,
                     .packages = c("lslx"),
                     .combine = rbind)%dopar%{
                       
                       source("functions/num_correct.R")
                       # Load data
                       load(paste0("DATA/data", i, ".RData"))
                       
                       X <- out$X # not standardized
                       y <- out$y # not standardized
                       X_std <- scale(X)
                       
                       # Population models
                       Ptrue = out$Ptrue
                       Htrue = out$Htrue
                       nfactors <- dim(Ptrue)[2]
                       J <- dim(Ptrue)[1]
                       
                       # syntax for lslx model
                       colnames(X_std) <- paste0("X", 1:J)
                       samplex <- seq(1:J)
                       
                       if (nfactors == 3){
                         line1 <-paste("f1:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line2 <-paste("f2:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line3 <- paste("f3:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         #linex <- paste("f1  <=>  0*f2")
                         #linex1 <- paste0("f2 <=> 0*f3")
                         #linex2 <- paste0("f3 <=> 0*f1")
                         line4 <-paste("f1  <=> fix(1)* f1")
                         line5 <-paste("f2  <=> fix(1)* f2")
                         line6 <-paste("f3  <=> fix(1)* f3")
                         
                         model_efa<- paste(paste(line1,line2,line3,line4,line5,line6, collapse = "\n", sep="\n"),collapse = "\n","\n")
                         lslx_efa <- lslx$new(model = model_efa,
                                              data = X_std)
                         
                         penalty_method <- "MCP"
                         converge <- TRUE
                         warning <- FALSE
                         
                         output <- capture.output({
                           lslx_efa$fit(penalty_method = "mcp")
                         })
                         
                         if (any(grepl("WARNING", output))){
                           warning <- TRUE
                         }
                         
                         # Look for the warning string in the output
                         if (any(grepl("ERROR", output)) | any(grepl("Hessian is not convex", output))) {
                           message("Use elastic net instead")
                           penalty_method <- "elastic_net"
                           elastic_output <- capture.output({lslx_efa$fit(penalty_method = "elastic_net")})
                           
                           if (any(grepl("WARNING", elastic_output))){
                             warning <- TRUE
                           }
                           if (any(grepl("ERROR", elastic_output))) {
                             converge <- FALSE
                           }
                         }
                         
                       } else {
                         line1 <-paste("f1:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line2 <-paste("f2:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line3 <- paste("f3:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line4 <-  paste("f4:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         line5 <-  paste("f5:~> ",paste("X",samplex[1:J],collapse = "+",sep = ""))
                         #linex <- paste("f1  <=>  0*f2")
                         #linex1 <- paste0("f1 <=> 0*f3")
                         #linex2 <- paste0("f1 <=> 0*f4")
                         #linex3<- paste0("f1 <=> 0*f5")
                         #linex4 <- paste0("f2 <=> 0*f3")
                         #linex5 <- paste0("f2 <=> 0*f4")
                         #linex6 <- paste0("f2 <=> 0*f5")
                         #linex7 <- paste0("f3 <=> 0*f4")
                         #linex8 <- paste0("f3 <=> 0*f5")
                         #linex9 <- paste0("f4 <=> 0*f5")
                         
                         line6 <-paste("f1  <=> fix(1)* f1")
                         line7 <-paste("f2  <=> fix(1)* f2")
                         line8 <-paste("f3  <=> fix(1)* f3")
                         line9 <-paste("f4  <=> fix(1)* f4")
                         line10 <-paste("f5  <=> fix(1)* f5")
                         
                         model_efa<- paste(paste(line1,line2,line3,line4,line5,
                                                 line6,line7,line8,line9,line10, collapse = "\n", sep="\n"),collapse = "\n","\n")
                         
                         lslx_efa <- lslx$new(model = model_efa,
                                              data = X_std)
                         
                         penalty_method <- "MCP"
                         converge <- TRUE
                         warning <- FALSE
                         
                         output <- capture.output({
                           lslx_efa$fit(penalty_method = "mcp")
                         })
                         
                         if (any(grepl("WARNING", output))){
                           warning <- TRUE
                         }
                         
                         
                         # Look for the warning string in the output
                         if (any(grepl("ERROR", output)) | any(grepl("Hessian is not convex", output))) {
                           message("Use elastic net instead")
                           penalty_method <- "elastic_net"
                           elastic_output <- capture.output({lslx_efa$fit(penalty_method = "elastic_net")})
                           
                           if (any(grepl("WARNING", elastic_output))){
                             warning <- TRUE
                           }
                           
                           if (any(grepl("ERROR", elastic_output)) | any(grepl("doesn't converge under EVERY", elastic_output)) | any(grepl("Hessian is not convex", elastic_output))) {
                             converge <- FALSE
                           }
                         }
                       }
                       
                       # extract coefficients
                       #lslx_efa$summarize(selector = "bic")
                       if (converge == FALSE){
                         Pmatrix_final <- NA
                         beta_reg <- NA
                         beta_barlett <- NA
                       } else {
                         Pmatrix <- lslx_efa$extract_coefficient_matrix(selector = "bic", block ="y<-f", include_faulty = TRUE)
                         Pmatrix_final <- Pmatrix$g
                       }
                       
                       if(all(is.na(Pmatrix_final))){
                         corrate_final = NA
                       } else {
                         # zero/non-zero structure recovery rate + Tucker's congruence
                         perm <- gtools::permutations(nfactors, nfactors)
                         corrate <- vector(length = nrow(perm))
                         
                         for (p in 1:nrow(perm)) {
                           corrate[p] <-num_correct(Ptrue, round(Pmatrix_final[,perm[p,]],3))
                           
                         }
                         corrate_final <- max(corrate)
                         
                         index <- which.max(corrate)
                         
                         
                         Pmatrix_final <- Pmatrix_final[,perm[index,]]
                         corsign <- sign(diag(cor(Ptrue, Pmatrix_final)))
                         
                         Pmatrix_final <- Pmatrix_final%*%diag(corsign)
                         
                         res <- lslx_efa$extract_coefficient_matrix(selector = "bic", block = "y<->y", include_faulty = TRUE)
                         res <- as.matrix(res$g)
                         
                         # barlett factor scores
                         P <- Pmatrix_final
                         H_barlett <- t(solve(t(P)%*%solve(res)%*%P)%*%t(P)%*%solve(res)%*%t(X_std))
                         beta_barlett <- lm(y ~ 0 + H_barlett)
                         
                         # regression factor scores
                         Phi  <- lslx_efa$extract_coefficient_matrix(selector = "bic", block = "f<->f", include_faulty = TRUE)
                         Phi <- as.matrix(Phi$g)
                         Sigma <- P %*% Phi %*% t(P) + res
                         H_reg <- X_std %*% solve(Sigma) %*% P %*% Phi
                         beta_reg <- lm(y ~ 0 + H_reg)
                         
                       }
                       
                       ## measure summary
                       measure_summary <- matrix(c(corrate_final, penalty_method, converge, warning), nrow = 1, byrow = T)
                       colnames(measure_summary) <- c("zero_rec", "penalty_method", "converge", "warning")
                       # Saving performance measures
                       output = list(measure_summary = measure_summary, 
                                     Pmatrix = Pmatrix_final, 
                                     beta_barlett = beta_barlett,
                                     beta_reg = beta_reg,
                                     lslx_model = lslx_efa)
                       
                       save(output, file = paste0("lslx_Sep/lslx", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)