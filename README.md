
Repository with R and C++ routines used in the paper: 

"A Bayesian random partition model for sequential refinement and coagulation"

available at https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13047

All functions are located in the root directory. 

Directory ./simulation contains the results of MCMC algorithm under multiple values of kappa1 and kappa2. 

Directories ./simulation/true_ab contain the results of estimation and simulated data with true values of kappa1 and kappa2 fixed as a and b, respectively (only truth_35 is uploaded here).

- mcmc_rppa.R : Main function. Use source("./your_path_to_mcmc_rppa.R_file/mcmc_rppa.R") in order to load the function in your R session
- choose_cluster.cpp : Auxiliary function for cluster membership estimation a posteriori
- run_scenario_2.R : runs multiple mcmc chains with the objective of choosing kappa1 and kappa2
- simulation_35.R : simulates data in section silmulation1.a of the paper (underlying truth: kappa1 = 3 kappa2 = 5). This file also runs multiple mcmc chains for choice of kappa1 and kappa2.
- GSK 2 PI3K AKT MEK inhib#32.xls: real dataset. Each sheet correspond to a different cell line
