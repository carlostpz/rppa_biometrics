
# This code simulates the data for scenario 1 under truth (kappa1, kappa2) = (3, 5) 

rm(list=ls(all=TRUE))

options(warn=2)
require(stargazer)

############################################
#          RPPA simulated data
############################################

v0u = 5

set.seed(1000)

require(MCMCpack)
require(mvtnorm)
require(foreach)
require(doParallel)
require(doMC)

CC <- 2
n <- 55
T <- 8
D <- 3
L <- 4
J <- 3

ones_T = rep(1, T)
t_ones_T = t(ones_T)

kappa1 = 3
kappa2 = 5

#initializing parameters

#Sigma
Sigma = array(0, dim = c(T,T,CC) )
for (cc in 1:CC){
  Sigma[,,cc] = 0.1*diag(T)
}
Sigma_inv = Sigma


# delta1
delta1 = array(0, dim=c(CC, D, n) )
for (cc in 1:CC){
  for (d in 1:D){
    delta1[cc,d, ] = sample(x=1:kappa1, size = n, replace = TRUE)
  }
}

# gamma
gamma = 0.9

# delta2
delta2 = array(2, dim=c(CC, D, n) )
for (cc in 1:CC){
  for (d in 1:D){
    for (i in 1:n){
      if ( runif(1) < gamma ){
        delta2[cc,d,i] = delta1[cc,d,i]
      }else{
        delta2[cc,d,i] = sample( x = (kappa1 + 1):kappa2 , size = 1, replace = TRUE)
      }
    }
  }
}


# delta
delta = array( 0, dim = c(2, CC, D, n))
delta[1, , , ] = delta1[ , , ] 
delta[2, , , ] = delta2[ , , ] 

# tau
tau = array(0, dim=c(CC, D, L, 2))
for (l in 1:L){
  tau[ , , l,1] = matrix( c(2,3,4,
                            2,3,4), ncol = 3, nrow = 2, byrow = TRUE )
  tau[ , , l,2] = matrix( c(5,6,7,
                            5,6,7), ncol = 3, nrow = 2, byrow = TRUE )
}

# K_cdu: Number of clusters by time period for each cell line and drug
Kcdu = array(0, dim = c(CC, D, 2) )
for (cc in 1:CC){
  for (d in 1: D){
    for (u in 1:2){
      Kcdu[cc, d, u] = length( unique(delta[u, cc, d, ]) )
    }
  }
}


mu_0u = c( 0.5, 1.5, 0.4 )

# mu*

# mus_cdl
mus = array( 0, dim=c(CC,D,L,3, kappa2) )
for (cc in 1:CC){
  for (d in 1:D){
    for (l in 1:L){
      mus[cc,d,l,1, 1:kappa1] = sort( rnorm( kappa1, mean = mu_0u[1], sd = 1/sqrt(v0u) ) )
      mus[cc,d,l,2, 1:kappa2] = sort( rnorm( kappa2, mean = mu_0u[2], sd = 1/sqrt(v0u) ) )
      mus[cc,d,l,3, 1:kappa1] = sort( rnorm( kappa1, mean = mu_0u[3], sd = 1/sqrt(v0u) ) )
    }
  }
}

# mu_cdli
mu_cdli = array(0, dim = c(CC, D, L, n, T) )

for(cc in 1:CC){
  for ( d in 1:D){
    for (l in 1:L){
      for( i in 1:n){
        mu_cdli [cc,d,l,i, 1:tau[cc, d, l, 1] ]                     = mus[ cc, d, l, 1, delta[1, cc, d, i]]
        mu_cdli [cc,d,l,i, ( tau[cc, d, l, 1] + 1 ):tau[cc,d,l,2] ] = mus[ cc, d, l, 2, delta[2, cc, d, i]]
        mu_cdli [cc,d,l,i, ( tau[cc, d, l, 2] + 1 ):T ]             = mus[ cc, d, l, 3, delta[1, cc, d, i]]
      }
    } 
  }
}

# Simulating the data

y_sim = array( 0, dim=c(CC, D, n, L, J, T) )

for (cc in 1:CC){
  for (d in 1:D){
    for (i in 1:n){
      for (l in 1:L){
        for (j in 1:J){
          y_sim [cc,d,i,l,j, ] = as.numeric( rmvnorm(1, mean = mu_cdli[cc,d,l,i, ], sigma = Sigma[ , ,cc] ) )
        }
      }
    }
  }
}

y = y_sim

# save workspace
save.image( file = "./simulation/workspace_simulation_35.Rdata")

# loading mcmc function
source( file = "mcmc_rppa.R")


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

# This part takes 6.5 hours to run, approximately. Uncomment if you wish to re-run mcmc
#system.time(
#mcmc_chains <- foreach( kappa1=kappa1_vec, kappa2=kappa2_vec ) %dopar% {
#  mcmc_rppa(n_iter=500, kappa1=kappa1, kappa2=kappa2, y=y, burn_in = 100)
#}
# returns a list where each element i is equal to vec[i]
#)[3]


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
save( chain_23, file = "./simulation/truth_35/chain_23.Rdata" )
ind = ind + 1

chain_24 = mcmc_chains[[ind]]
chain = chain_24
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_24, file = "./simulation/truth_35/chain_24.Rdata" )
ind = ind + 1

chain_25 = mcmc_chains[[ind]]
chain = chain_25
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_25, file = "./simulation/truth_35/chain_25.Rdata" )
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
save( chain_34, file = "./simulation/truth_35/chain_34.Rdata" )
ind = ind + 1

chain_35 = mcmc_chains[[ind]]
chain = chain_35
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_35, file = "./simulation/truth_35/chain_35.Rdata" )
ind = ind + 1

chain_36 = mcmc_chains[[ind]]
chain = chain_36
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_36, file = "./simulation/truth_35/chain_36.Rdata" )
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
save( chain_45, file = "./simulation/truth_35/chain_45.Rdata" )
ind = ind + 1

chain_46 = mcmc_chains[[ind]]
chain = chain_46
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_46, file = "./simulation/truth_35/chain_46.Rdata" )
ind = ind + 1

chain_47 = mcmc_chains[[ind]]
chain = chain_47
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_47, file = "./simulation/truth_35/chain_47.Rdata" )
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
save( chain_56, file = "./simulation/truth_35/chain_56.Rdata" )
ind = ind + 1

chain_57 = mcmc_chains[[ind]]
chain = chain_57
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_57, file = "./simulation/truth_35/chain_57.Rdata" )
ind = ind + 1

chain_58 = mcmc_chains[[ind]]
chain = chain_58
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_58, file = "./simulation/truth_35/chain_58.Rdata" )
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
save( chain_67, file = "./simulation/truth_35/chain_67.Rdata" )
ind = ind + 1

chain_68 = mcmc_chains[[ind]]
chain = chain_68
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_68, file = "./simulation/truth_35/chain_68.Rdata" )
ind = ind + 1

chain_69 = mcmc_chains[[ind]]
chain = chain_69
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_69, file = "./simulation/truth_35/chain_69.Rdata" )
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
save( chain_78, file = "./simulation/truth_35/chain_78.Rdata" )
ind = ind + 1

chain_79 = mcmc_chains[[ind]]
chain = chain_79
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_79, file = "./simulation/truth_35/chain_79.Rdata" )
ind = ind + 1

chain_710 = mcmc_chains[[ind]]
chain = chain_710
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_710, file = "./simulation/truth_35/chain_710.Rdata" )
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
save( chain_89, file = "./simulation/truth_35/chain_89.Rdata" )
ind = ind + 1

chain_810 = mcmc_chains[[ind]]
chain = chain_810
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_810, file = "./simulation/truth_35/chain_810.Rdata" )
ind = ind + 1

chain_811 = mcmc_chains[[ind]]
chain = chain_811
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
lpml_vec[ind] = chain$lpml
save( chain_811, file = "./simulation/truth_35/chain_811.Rdata" )
ind = ind + 1


###############################
#    Model comparison: table
###############################

model_comparison = rbind( round( rbind( aic_vec, bic_vec, dic_vec, waic_vec), 3), lpml_vec)
model_comparison = round( model_comparison[,1:21] )

# saving results
write.table(model_comparison, row.names = TRUE, col.names = FALSE, 
            file = "./simulation/truth_35/model_comparison_truth_35.txt")

colnames( model_comparison ) = c("23", "24", "25", 
                                 "34", "35", "36",
                                 "45", "46", "47",
                                 "56", "57", "58",
                                 "67", "68", "69",
                                 "78", "79", "710",
                                 "89", "810", "811")

# latex table
stargazer( model_comparison )
