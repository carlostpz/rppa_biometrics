
Repository with R and C++ routines used in the paper: 

"A Bayesian random partition model for sequential refinement and coagulation"

(Include link to Biometrics website)

All functions are located in the root directory. 
Directory ./simulation contains the results of MCMC algorithm under multiple values of kappa1 and kappa2. 
Directories ./simulation/true_ab contain the results of estimation and simulated data with true values of kappa1 and kappa2 fixed as a and b, respectively.

- Main function: mcmc_rppa.R 
- Auxiliary function for cluster membership estimation a posteriori: choose_cluster.cpp
