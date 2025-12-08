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
- [IS.R](Functions/IS.R) and [IS-lasso.R](Functions/IS-lasso.R) for model selection with the Index of Spareness for the cardinality-constrained and LASSO approaches, respectively.
- [findLasso.R](Functions/findLasso.R) for binary search procedure to tune LASSO penalty parameter to reach a prespecified number of (non)zero loadings. Note that it does not always work perfectly!

## Simulation Study
Folder [Simulation Study](Simulation Study) contains scripts to generate the simulated data and script for each method to analyze the simulated data:
- Run [data_run.R](Simulation Study/data_run.R) to generate the simulated data.

Additionally, the subfolder [RegSEM-failed](https://github.com/trale97/RegularizedESEM/tree/main/Simulation%20Study/RegSEM-failed) contains codes for RegSEM to estimate the relationships among the factors, which produced problematic results. 

## Empirical Data Application
The folder [Application](Application) contains an R script to analyze the Gene Expression data using the two proposed methods. 
- The Gene Expression data set can be obtained from the NCBI GEO database using GEO accession number GSE7329 (or download directly from [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz)).
- Readers interested in the application using the Stress Recovery Strategies data, please directly contact me at t.t.le_1@tilburguniversity.edu.
