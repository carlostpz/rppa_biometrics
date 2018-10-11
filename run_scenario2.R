
'
This code
1) calls mcmc main function for first run (i., e., sampling delta) from: mcmc_rppa_ver.R
2) applies it to simulated data under scenario 2 
3) saves each one of the chains, for example:

save( chain_78, file = "./simulation2/chain_78_comparison.Rdata" )

4) runs the mcmc for the best model for more iterations and saves it:

save( chain_46, file = "./siulation2/best_chain_46.Rdata")
'

rm(list=ls(all=TRUE))

set.seed(100)

require(MCMCpack)
require(mvtnorm)
require(foreach)
require(doParallel)
require(doMC)
require(stargazer)

# load main mcmc function (for doing the first run)
source("mcmc_rppa.R")

# Loading the simulated data for scenario 2
load("./simulation2/y_sim.Rdata")


#############################
#     model comparison
#############################

######################
#    DIC function
######################

# calculating theta_hat

calc_dic = function( chain, y, init_iter, end_iter){
  
  CC <- 2
  n <- 55
  T <- 8
  D <- 3
  L <- 4
  J <- 3
  
  # outlier with meaninglessly high expression
  # imputing it with
  y[2, 1, 15, 4, 1, 1] = mean( y[2, 1, 15, 4, 2:3, 1] )
  
  mu_hat = array(0, dim = c(CC, D, L, n, T) )
  Sigma_hat = array(0, dim = c(T, T, CC) )
  
  index = init_iter:end_iter
  
  mu = chain$mu[ , , , , , index ]
  Sigma = chain$Sigma[ , , , index ]
  
  for (cc in 1:CC){
    for (d in 1:D){
      for (i in 1:n){
        for (l in 1:L){
          for (tt in 1:T){
            mu_hat[cc, d, l, i, tt] = mean( mu[cc,d,l,i,tt, ] )
          }
        }
      }
    }
  }
  
  for (cc in 1:CC){
    for (t1 in 1:T){
      for (t2 in 1:T){
        Sigma_hat[t1 , t2, cc] = mean( Sigma[t1, t2, cc, ])
      }
    }
  }
  
  # calculating mean of log-likelihood
  E_log_post = mean( chain$log_lik[ index ] )
  # calculating D_bar
  D_bar = -2*E_log_post
  
  # calculating log likelihood of theta hat
  log_lik_theta_hat = 0
  for(cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        for (i in 1:n){
          for (j in 1:J){
            expres = ( y[cc,d,i,l,j, ] - mu_hat[cc,d,l,i, ] )
            log_lik_theta_hat = log_lik_theta_hat - 0.5 * T * log( 2 * pi ) -
              0.5 * log( det(Sigma_hat[ , , cc]) ) - 
              0.5 * t(expres) %*% solve(Sigma_hat[,,cc]) %*% expres
          }
        }
      }
    }
  }
  
  # calculatin D_theta_bar
  D_theta_bar = -2 * log_lik_theta_hat
  
  # calculating pD
  pD = D_bar - D_theta_bar
  
  # DIC
  DIC = D_bar + pD # the smaller the better
  
  return ( as.numeric(DIC) )
  
}


bic_vec = rep(0, 21)
aic_vec = rep(0, 21)
dic_vec = rep(0, 21)
waic_vec = rep(0, 21)
margin_lik_vec = rep(0, 21)
lpml_vec = rep(0, 21)
ind = 1


cl <- makeCluster(4)
registerDoParallel(cl)

kappa1_vec = c(2,2,2,
               3,3,3,
               4,4,4,
               5,5,5,
               6,6,6,
               7,7,7,
               8,8,8)

kappa2_vec = c(3,4,5,
               4,5,6,
               5,6,7,
               6,7,8,
               7,8,9,
               8,9,10,
               9,10,11)


# This part takes 6.5 hours to run, approximately
system.time(
  mcmc_chains <- foreach( kappa1=kappa1_vec, kappa2=kappa2_vec ) %dopar% {
    mcmc_rppa(n_iter=500, kappa1=kappa1, kappa2=kappa2, y=y, burn_in = 100)
  }
  # returns a list where each element i is equal to vec[i]
)[3]


bic_vec = rep(0, 21)
aic_vec = rep(0, 21)
dic_vec = rep(0, 21)
waic_vec = rep(0, 21)
lpml_vec = rep(0, 21)
ind = 1


# gathering the different chains

###############
#  kappa1 = 2
###############
chain_23 = mcmc_chains[[ind]]
chain = chain_23
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_23, file = "./simulation2/chain_23.Rdata" )
ind = ind + 1

chain_24 = mcmc_chains[[ind]]
chain = chain_24
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_24, file = "./simulation2/chain_24.Rdata" )
ind = ind + 1

chain_25 = mcmc_chains[[ind]]
chain = chain_25
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_25, file = "./simulation2/chain_25.Rdata" )
ind = ind + 1


###############
#  kappa1 = 3
###############
chain_34 = mcmc_chains[[ind]]
chain = chain_34
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_34, file = "./simulation2/chain_34.Rdata" )
ind = ind + 1

chain_35 = mcmc_chains[[ind]]
chain = chain_35
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_35, file = "./simulation2/chain_35.Rdata" )
ind = ind + 1

chain_36 = mcmc_chains[[ind]]
chain = chain_36
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_36, file = "./simulation2/chain_36.Rdata" )
ind = ind + 1


###############
#  kappa1 = 4
###############
chain_45 = mcmc_chains[[ind]]
chain = chain_45
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
save( chain_45, file = "./simulation2/chain_45.Rdata" )
ind = ind + 1

chain_46 = mcmc_chains[[ind]]
chain = chain_46
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_46, file = "./simulation2/chain_46.Rdata" )
ind = ind + 1

chain_47 = mcmc_chains[[ind]]
chain = chain_47
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_47, file = "./simulation2/chain_47.Rdata" )
ind = ind + 1


###############
#  kappa1 = 5
###############
chain_56 = mcmc_chains[[ind]]
chain = chain_56
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_56, file = "./simulation2/chain_56.Rdata" )
ind = ind + 1

chain_57 = mcmc_chains[[ind]]
chain = chain_57
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_57, file = "./simulation2/chain_57.Rdata" )
ind = ind + 1

chain_58 = mcmc_chains[[ind]]
chain = chain_58
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_58, file = "./simulation2/chain_58.Rdata" )
ind = ind + 1


###############
#  kappa1 = 6
###############
chain_67 = mcmc_chains[[ind]]
chain = chain_67
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_67, file = "./simulation2/chain_67.Rdata" )
ind = ind + 1

chain_68 = mcmc_chains[[ind]]
chain = chain_68
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_68, file = "./simulation2/chain_68.Rdata" )
ind = ind + 1

chain_69 = mcmc_chains[[ind]]
chain = chain_69
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_69, file = "./simulation2/chain_69.Rdata" )
ind = ind + 1


###############
#  kappa1 = 7
###############
chain_78 = mcmc_chains[[ind]]
chain = chain_78
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_78, file = "./simulation2/chain_78.Rdata" )
ind = ind + 1

chain_79 = mcmc_chains[[ind]]
chain = chain_79
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_79, file = "./simulation2/chain_79.Rdata" )
ind = ind + 1

chain_710 = mcmc_chains[[ind]]
chain = chain_710
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_710, file = "./simulation2/chain_710.Rdata" )
ind = ind + 1


###############
#  kappa1 = 8
###############
chain_89 = mcmc_chains[[ind]]
chain = chain_89
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_89, file = "./simulation2/chain_89.Rdata" )
ind = ind + 1

chain_810 = mcmc_chains[[ind]]
chain = chain_810
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_810, file = "./simulation2/chain_810.Rdata" )
ind = ind + 1

chain_811 = mcmc_chains[[ind]]
chain = chain_811
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_811, file = "./simulation2/chain_811.Rdata" )
ind = ind + 1


###############################
#    Model comparison: table
###############################

model_comparison = rbind( round( rbind( aic_vec, bic_vec, dic_vec, waic_vec), 3), lpml_vec)
model_comparison = round( model_comparison[,1:21] )

# saving results
write.table(model_comparison, row.names = TRUE, col.names = FALSE, 
            file = "./simulation2/model_comparison_truth_35.txt")

colnames( model_comparison ) = c("23", "24", "25", 
                                 "34", "35", "36",
                                 "45", "46", "47",
                                 "56", "57", "58",
                                 "67", "68", "69",
                                 "78", "79", "710",
                                 "89", "810", "811")

# latex table
stargazer( model_comparison )

save(model_comparison, file = "./simulation2/model_comparison_table.Rdata")

