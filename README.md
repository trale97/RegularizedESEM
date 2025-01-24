Functions, simulation study and empirical data application for Regularized Exploratory SEM method (Regularized ESEM).

## Authors
- Tra T. Le 
- Prof. dr. Jeroen Vermunt
- Dr. Nicola Ballhausen
- Prof. dr. Katrijn Van Deun

## Description
The repository consists of accompanying codes for manuscript *Exploratory Structural Equation Modeling and the Curse of Dimensionality*.
Here, you can find the **functions** developed in R to implement the method, the codes for the **Simulation Study** and **Empirical Data Analysis**.

## Functions
Folder [Functions](Functions) contains all the functions needed for the two methods: 
- [rESEM-CC.R](Functions/rESEM-CC.R) for the cardinality constrained approach.
- [rESEM_l1.R](Functions/rESEM_l1.R) for the one with the LASSO penalty.
- [rESEMl1_undoshrinkage.R](Functions/rESEMl1_undoshrinkage.R) to undo the shrinkage for nonzero loadings.
The first two files also include functions to run **Model Selection** based on the Index of Spareness.

## Simulation Study
- Step 1: Run [corrmatrix_orthogonal.R](corrmatrix_orthogonal.R) to generate the population correlation matrix
- Step 2: Run [datageneration_orthogonal.R](datageneration_orthogonal.R) and [datageneration_correlated](datageneration_correlated.R) to generate data for orthogonal and correlated latent variables scenarios, respectively. 
- Step 3: Run one of the sim_*method*.R file to analyze the simulated data with one of the methods. For example, run [sim_CCrESEM.R](sim_CCrESEM.R) to run the analysis using the proposed method with cardinality constraint.

Additionally, the subfolder [RegSEM-failed](https://github.com/trale97/RegularizedESEM/tree/main/Simulation%20Study/RegSEM-failed) contains codes for RegSEM to estimate the relationships among the factors, which produced problematic results. 

## Empirical Data Application
The folder [Application](Application) contains an R script to analyze the Gene Expression data using the two proposed methods. 
- The Gene Expression data set can be obtained from the NCBI GEO database using GEO accession number GSE7329 (or download directly from [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz)).
- Readers interested in the application using the Stress Recovery Strategies data, please directly contact me at t.t.le_1@tilburguniversity.edu.
