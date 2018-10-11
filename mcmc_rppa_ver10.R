
rm(list=ls(all=TRUE))

set.seed(100)

mcmc_rppa = function(n_iter, kappa1, kappa2, y, burn_in){
  
  require(Rcpp)
  require(RcppArmadillo)
  require(MCMCpack)
  require(msm)
  require(Matrix)
  require(mvtnorm)
  
  CC <- 2
  n <- 55
  T <- 8
  D <- 3
  L <- 4
  J <- 3
  
  # outlier with meaninglessly high expression
  # imputing it with
  y[2, 1, 15, 4, 1, 1] = mean( y[2, 1, 15, 4, 2:3, 1] )
  
  # number of data points and parameters
  # these are to be used to calculate AIC, BIC, DIC
  n_data = CC*D*L*n*J*T
  # n_param = 11 + CC*D*L*( 2*kappa1 + kappa2 + 2 + 2*n )
  n_param = CC*D*L*( 2*kappa1 + kappa2 )
  
  ones_T = rep(1, T)
  t_ones_T = t(ones_T)
  
  chain_log_lik = rep(0, n_iter)
  
  ##################################
  # Choosing prior hyperparameters
  ##################################
  
  v0 = 10
  mu00 = 0
  v00 = 1/(1.5^2)
  av = 1
  bv = 1
  a_gamma = 1
  b_gamma = 1
  eta1 = rep(1, kappa1)
  eta2 = rep(1, kappa2 - kappa1)
  v_Sigma = 10
  V_Sigma = diag(T)
  
  #############################
  #  initializing parameters  
  #############################
  
  #Sigma
  Sigma = array(0, dim = c(T,T,CC) )
  for (cc in 1:CC){
    Sigma[,,cc] = diag(T)
  }
  Sigma_inv = Sigma
  
  chain_Sigma = array(0, dim = c(T, T, CC, n_iter) )
  chain_Sigma_inv = array(0, dim = c(T, T, CC, n_iter) )
  
  #######################################################
  # starting with every observation in the same cluster
  #######################################################
  
  # delta1
  
  #for (cc in 1:CC){
  #  for (d in 1:D){
  #  # delta1[cc,d, ] = sample(x=1:kappa1, size = n, replace = TRUE)
  #  delta1[cc,d, ] = 1
  #  }
  #}
  
  # delta2
  #delta2 = array(1, dim=c(CC, D, n) )
  #for (cc in 1:CC){
  #  for (d in 1:D){
  #    for (i in 1:n){
  #      # delta2[cc,d, i] = sample( x = c( delta1[ cc, d, i], (kappa1 + 1):kappa2 ) , size = 1, replace = TRUE)
  #      delta2[cc,d, ] = 1
  #    }
  #  }
  #}
  
  ########
  # delta
  ########
  
  delta = array( 0 , dim = c( 2, CC, D, n ) )
  mat_clust = array( 0, dim = c( CC, D, n, T ) )
  
  for (cc in 1:CC){
    for (d in 1:D){
      for (i in 1:n){
        for (tt in 1:T){
          mat_clust[cc, d, i, tt] = mean ( y[cc,d,i, , ,tt] )
          Dist_mat = dist( mat_clust[ cc, d, , ] )
          cluster =  hclust( d = Dist_mat ) 
          delta[1, cc, d, ] = cutree( cluster, k=kappa1 )
        }
      }
    }
  }
  delta[2, , , ] = delta[1, , , ]
  
  chain_delta = array( 0, dim = c(2, CC, D, n, n_iter) )
  
  ############
  #  tau_cdl
  ############
  
  tau = array(1, dim=c(CC, D, L, 2))
  tau[ , , , 2] = 3
  
  chain_tau = array( 0, dim = c( CC, D, L, 2, n_iter ) )
  
  ############
  #  mu*_cdl
  ############
  
  # IMPORTANT: the last kappa2 - kappa1 entries of mus_cdl1 and mus_cdl3 must never be used in the algorithm 
  mus = array(0, dim=c( CC, D, L, 3, kappa2 ) )
  chain_mu_star = array( 0, dim = c( CC, D, L, 3, kappa2, n_iter ) )
  
  ###########
  # mus_cdl1
  ###########
  
  # Initialize as the mean of the protein expressions in the initial time ( t=1 ) 
  for (cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        mus[ cc, d, l, 1, 1:kappa1] = mean( y[cc,d, , l, , 1] )
      }
    }
  }
  
  ###########
  # mus_cdl2
  ###########
  
  # Initialize as the mean of the protein expressions in the "middle" time ( t=4 ) 
  for (cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        mus[ cc, d, l, 2, 1:kappa2] = mean( y[cc,d, , l, , 4] )
      }
    }
  }
  
  ###########
  # mus_cdl3
  ###########
  
  # Initialize as the mean of the protein expressions in the final time ( t=8 ) 
  for (cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        mus[ cc, d, l, 3, 1:kappa1] = mean( y[cc,d, , l, , 8] )
      }
    }
  }
  
  # initializing mu_cdli
  chain_mu_cdli = array( 0, dim = c(CC, D, L, n, T, n_iter) )
  
  
  #########
  # gamma 
  #########
  
  gamma = 0.5
  chain_gamma = rep( 0, n_iter)
  
  #######
  # pi2
  #######
  
  pi2 = rep(1/(kappa2 - kappa1), kappa2 - kappa1)
  chain_pi2 = matrix(0, ncol = length(pi2), nrow = n_iter)
  
  #######
  # pi1
  #######
  
  pi1 = rep(1/kappa1, kappa1)
  chain_pi1 = matrix(0, ncol = length(pi1), nrow = n_iter)
  
  ######
  # v0u 
  ######
  
  v0u = rep(1, 3)
  chain_v0u = matrix(0, ncol = length(v0u), nrow = n_iter)
  
  #######
  # mu0u
  #######
  
  mu0u = c(0, 0, 0)
  chain_mu0u = matrix(0, ncol = length(mu0u), nrow = n_iter)
  
  
  #######################################################################
  #
  #                                MCMC
  #
  #######################################################################
  
  # Profile:
  #Rprof("/home/tadeu/ut_austin/RPPA/depRPPA/inference/rppa_prof.txt")
  
  #################
  #   subrutines  #
  #################
  
  # Pre cashing sum_j y_cdilt
  sum_y_over_j = array(0, dim = c(CC, D, n, L, T) )
  for ( cc in 1:CC){
    for (d in 1:D){
      for (i in 1:n){
        for (l in 1:L){
          for (time in 1:T){
            sum_y_over_j[cc,d, i, l, time] = sum( y[cc,d,i,l, ,time] )
          }
        }
      }
    }
  }
  
  ############################################################################
  #########################   Marginalization part    ########################
  ############################################################################
  
  # function that builds the vector mu_cdli (general version of the function)
  build_mu_cdli_margin = function( mus_cdl1, mus_cdl2, mus_cdl3, 
                                  tau_cdl1, tau_cdl2, 
                                  delta_1cdi, delta_2cdi, T=8){
    
    new_mu_cdli = rep(0, T)
    
    new_mu_cdli [1:tau_cdl1 ]              = mus_cdl1[ delta_1cdi ]
    new_mu_cdli [( tau_cdl1 + 1 ):tau_cdl2] = mus_cdl2[ delta_2cdi ]
    new_mu_cdli [( tau_cdl2 + 1 ):T ]      = mus_cdl3[ delta_1cdi ]
    
    return(new_mu_cdli)
    
  }
  
  
  # building the design matrix: X_cdil
  build_X = function( delta_1cdi, delta_2cdi, tau_cdl1, tau_cdl2, kappa1, kappa2, T=8 ){
    
    tau1 = tau_cdl1
    tau2 = tau_cdl2
    delta1 = delta_1cdi
    delta2 = delta_2cdi
    
    X = matrix(0, ncol = 2*kappa1 + kappa2, nrow = T)
    X [ 1:tau1, delta1] = 1
    X [ (tau1+1):tau2, delta2 + kappa1 ] = 1
    X [ (tau2+1):T, delta1 + kappa1 + kappa2 ] = 1

    return(X)
    
  }
  
  ######################
  #  sampling tau_cdl1
  ######################
  sampling_tau_cdl1 = function( l, mu0u, v0u,  
                                delta_1cd, delta_2cd, Sigma_inv_c, tau_cdl2,
                                kappa1, kappa2, sum_y_over_j_cd, y_cd, T=8){
    
    if (tau_cdl2 == 2){
      return(1)
    }else{
      
      pmf_tau_cdl1 = rep( 0, tau_cdl2-1 )
      
      for ( m in 1:( tau_cdl2 - 1 ) ){
        
        #  prior mean and precision matrix
        m_cdl = c( rep(mu0u[1], kappa1), rep(mu0u[2], kappa2), rep(mu0u[3], kappa1) )
        S_cdl_inv = diag( c( rep(v0u[1], kappa1), rep(v0u[2], kappa2), rep(v0u[3], kappa1) ) ) 
        
        S_cdl_star = S_cdl_inv
        m_cdl_star = S_cdl_inv %*% m_cdl
        
        for (i in 1:n){ # loop through proteins "i"

            X_cdil = build_X( delta_1cdi = delta_1cd[i], 
                              delta_2cdi = delta_2cd[i], 
                              tau_cdl1 = m, 
                              tau_cdl2 = tau_cdl2,
                              kappa1 = kappa1, 
                              kappa2 = kappa2 )
            
            # posterior covariance matrix
            S_cdl_star = S_cdl_star + J * t(X_cdil) %*% Sigma_inv_c %*% X_cdil
            
            # posterior mean vector
            m_cdl_star = m_cdl_star + t(X_cdil) %*% Sigma_inv_c %*% sum_y_over_j_cd[i,l, ]

        } # end of loop through proteins "i"
        
        # finalizing 
        S_cdl_star = solve(S_cdl_star) 
        m_cdl_star = S_cdl_star %*% m_cdl_star
        
        # calculating the "likelihood"
        sum_aux = 0
        for (i in 1:n){
          mu_cdli_aux = build_mu_cdli_margin( mus_cdl1 = m_cdl_star[1:kappa1],
                                              mus_cdl2 = m_cdl_star[ (kappa1 + 1) :( kappa1 + kappa2) ],
                                              mus_cdl3 = m_cdl_star[ ( kappa1 + kappa2 + 1): (2*kappa1 + kappa2)],
                                              tau_cdl1 = m,
                                              tau_cdl2 = tau_cdl2,
                                              delta_1cdi = delta_1cd[i],
                                              delta_2cdi = delta_2cd[i] )
          
          for (j in 1:J){
            expres = y_cd[i, l, j, ] - mu_cdli_aux
            sum_aux = sum_aux + t( expres ) %*% Sigma_inv_c %*% expres
          }
          
        }
        log_lik_part = as.numeric(sum_aux)
        
        pmf_tau_cdl1[m]=
          0.5 * log( det( S_cdl_star ) ) - 0.5* (
            v0u[1] * sum( (m_cdl_star[ 1:kappa1] - mu0u[1] )^2 ) +
              v0u[2] * sum( (m_cdl_star[ (kappa1 + 1):(kappa1 + kappa2)] - mu0u[2])^2 ) +
              v0u[3] * sum( (m_cdl_star[ (kappa1 + kappa2 + 1):(2*kappa1 + kappa2)] - mu0u[3])^2 ) +
              log_lik_part
          )
        
      }
      
      # normalizing the weights
      pmf_tau_cdl1 = exp(pmf_tau_cdl1 - max(pmf_tau_cdl1) )
      pmf_tau_cdl1 = pmf_tau_cdl1/sum(pmf_tau_cdl1)
      
      # sampling tau_cd1
      tau_sample = sample( size = 1, x = 1:(tau_cdl2-1), prob = pmf_tau_cdl1 )
      
      return( tau_sample )
      
    }
    
  }
  
  ######################
  #  sampling tau_cdl2
  ######################
  
  sampling_tau_cdl2 = function( l, mu0u, v0u, delta_1cd, delta_2cd, 
                                Sigma_inv_c, tau_cdl1,
                                kappa1, kappa2, sum_y_over_j_cd, y_cd, T=8){
    
    if(tau_cdl1 == 6){
      return(7)
    }else{
      
      pmf_tau_cdl2 = rep( 0, T - tau_cdl1 - 1 )
      
      for ( m in ( tau_cdl1 + 1 ):(T-1) ){
        
        #  prior mean and precision matrix
        m_cdl = c( rep(mu0u[1], kappa1), rep(mu0u[2], kappa2), rep(mu0u[3], kappa1) )
        S_cdl_inv = diag( c( rep(v0u[1], kappa1), rep(v0u[2], kappa2), rep(v0u[3], kappa1) ) )
        
        S_cdl_star = S_cdl_inv
        m_cdl_star = S_cdl_inv %*% m_cdl
        for (i in 1:n){
          X_cdil = build_X( delta_1cdi = delta_1cd[i], 
                            delta_2cdi = delta_2cd[i], 
                            tau_cdl1 = tau_cdl1, 
                            tau_cdl2 = m,
                            kappa1 = kappa1, kappa2 = kappa2 )
            
          # posterior covariance matrix
          S_cdl_star = S_cdl_star + J * t(X_cdil) %*% Sigma_inv_c %*% X_cdil
          
          # posterior mean vector
          m_cdl_star = m_cdl_star + t(X_cdil) %*% Sigma_inv_c %*% sum_y_over_j_cd[i,l, ]
        }
        # finalizing 
        S_cdl_star = solve(S_cdl_star) 
        m_cdl_star = S_cdl_star %*% m_cdl_star
        
        # calculating the "likelihood"
        sum_aux = 0
        for (i in 1:n){
          mu_cdli_aux = build_mu_cdli_margin( mus_cdl1 = m_cdl_star[ 1 : kappa1 ],
                                              mus_cdl2 = m_cdl_star[ (kappa1 + 1) : ( kappa1 + kappa2) ],
                                              mus_cdl3 = m_cdl_star[ ( kappa1 + kappa2 + 1) : (2*kappa1 + kappa2) ],
                                              tau_cdl1 = tau_cdl1,
                                              tau_cdl2 = m,
                                              delta_1cdi = delta_1cd[i],
                                              delta_2cdi = delta_2cd[i] )
          
          for (j in 1:J){
            expres = y_cd[i, l, j, ] - mu_cdli_aux
            sum_aux = sum_aux + t( expres ) %*% Sigma_inv_c %*% expres
          }
          
        }
        log_lik_part = as.numeric(sum_aux)
        
        pmf_tau_cdl2[ m - tau_cdl1 ] = 
          0.5 * log( det(S_cdl_star) ) - 0.5* (
            v0u[1] * sum( (m_cdl_star[ 1:kappa1] - mu0u[1] )^2 ) +
              v0u[2] * sum( (m_cdl_star[ (kappa1 + 1):(kappa1 + kappa2)] - mu0u[2])^2 ) +
              v0u[3] * sum( (m_cdl_star[ (kappa1 + kappa2 + 1):(2*kappa1 + kappa2)] - mu0u[3])^2 ) +
              log_lik_part
          )
        
      }
      
      # normalizing the weights
      pmf_tau_cdl2 = exp(pmf_tau_cdl2 - max(pmf_tau_cdl2) )
      pmf_tau_cdl2 = pmf_tau_cdl2/sum(pmf_tau_cdl2)
      
      # sampling tau_cd2
      tau_sample = sample( size = 1, x = (tau_cdl1+1):(T-1), prob = pmf_tau_cdl2 )
      
      return( tau_sample )
    }    
  }
  
  ##################
  # sample_mus_cd1
  ##################
  
  sample_mus_cdl1 = function (l, m, Sigma_c_inv, 
                              mus_cdl1, mus_cdl2, mus_cdl3, 
                              tau_cdl1, tau_cdl2,  
                              delta_1cd, delta_2cd, 
                              y_cd, kappa1, T=8 ){
    
    ones_T = rep(1, T)
    
    u1 = rep(0, T)
    u2 = rep(0, T)
    u3 = rep(0, T)
    
    u1[ 1 : tau_cdl1 ] = 1
    u2[ (tau_cdl1 + 1) : tau_cdl2] = 1
    u3[ (tau_cdl2 + 1) : T] = 1
    
    t_u1 = t(u1)
    
    u = 1
    mus_cdl3_m = mus_cdl3[m]
    
    set_P = which( delta_1cd == m )
    hash_P = length( set_P )
    
    bfy = rep(0, T)
    for (time in 1:T){
      bfy[time] =  sum( y_cd[set_P, l, , time] )
    }
    
    b1 = 1/( v0u[u] + J * hash_P * t_u1 %*% Sigma_c_inv %*% u1 )
    a1 = b1 * ( t_u1 %*% Sigma_c_inv %*%  
                  (bfy - J * u2 * sum( mus_cdl2[ delta_2cd[ set_P ] ] ) - 
                   J * hash_P * u3 * mus_cdl3_m ) + mu0u[u]*v0u[u] )
    
    if( m == 1 ){
      return( rtnorm(1, mean=a1, sd = sqrt(b1), upper=mus_cdl1[ m+1 ] ) )
    }else{
      if (m == kappa1){
        return( rtnorm(1, mean=a1, sd = sqrt(b1), lower=mus_cdl1[ m-1 ] ) )
      }else{
        return( rtnorm(1, mean=a1, sd = sqrt(b1), lower=mus_cdl1[ m-1 ], upper=mus_cdl1[ m+1 ] ) )
      } 
    }
    
  }
  
  ####################
  #  sample_mus_cd2
  ####################
  
  sample_mus_cdl2 = function ( l, m, Sigma_c_inv, 
                              mus_cdl1, mus_cdl3, 
                              tau_cdl1, tau_cdl2,  
                              delta_1cd, delta_2cd,
                              y_cd, T=8 ){
    
    ones_T = rep(1, T)
    
    u1 = rep(0, T)
    u2 = rep(0, T)
    u3 = rep(0, T)
    
    u1[1: tau_cdl1] = 1
    u2[ (tau_cdl1 + 1):tau_cdl2] = 1
    u3[ (tau_cdl2 + 1):T] = 1
    
    u = 2
    t_u2 = t(u2)
    
    set_P = which( delta_2cd ==m )
    hash_P = length( set_P )
    
    bfy = rep(0, T)
    for (time in 1:T){
      bfy[time] =  sum( y_cd[ set_P, l, , time ] )
    }
    
    b2 = 1/( v0u[u] + J*hash_P*t_u2%*%Sigma_c_inv%*%u2 )
    a2 = b2 * ( t_u2 %*% Sigma_c_inv %*% 
                  ( bfy - 
                      J * u1 * sum( mus_cdl1[ delta_1cd[set_P] ] ) - 
                      J * u3 * sum( mus_cdl3[ delta_1cd[set_P] ] ) ) + mu0u[u]*v0u[u] )
    
    return( rnorm(1, mean=a2, sd = sqrt(b2) ) )
    
  }
  
  #################
  # sample_mus_cd3
  #################
  
  sample_mus_cdl3 = function (l, m, Sigma_c_inv, 
                              mus_cdl1, mus_cdl2, 
                              tau_cdl1, tau_cdl2,  
                              delta_1cd, delta_2cd,
                              y_cd, T=8 ){
    
    ones_T = rep(1, T)
    
    u1 = rep(0, T)
    u2 = rep(0, T)
    u3 = rep(0, T)
    
    u1[ 1 : tau_cdl1 ] = 1
    u2[ ( tau_cdl1 + 1):tau_cdl2 ] = 1
    u3[ ( tau_cdl2 + 1):T ] = 1
    
    t_u3 = t(u3)
    
    u = 3
    mus_cdl1_m = mus_cdl1[m]
    
    set_P = which( delta_1cd == m )
    hash_P = length( set_P )
    
    bfy = rep(0, T)
    for (time in 1:T){
      bfy[time] =  sum( y_cd[set_P, l, , time] )
    }
    
    b3 = 1/( v0u[u] + J*hash_P*t_u3%*%Sigma_c_inv%*%u3 )
    a3 = b3 * (  t_u3%*% Sigma_c_inv %*% 
                   ( bfy - 
                       J * u2 * sum( mus_cdl2[ delta_2cd[ set_P ] ] ) - 
                       J * hash_P * u1 * mus_cdl1_m ) + mu0u[u]*v0u[u] )
    
    return( rnorm(1, mean=a3, sd = sqrt(b3) ) )
    
  }
  
  ##########################
  #  sampling delta_1cdi 
  ##########################
  
  
  sample_delta_1cdi = function( kappa1, kappa2, 
                                mus_cd1, mus_cd2, mus_cd3, 
                                tau_cd1, tau_cd2, delta_2cdi, 
                                Sigma_inv_c, y_cdi, pi1){
    
    # mus_cd1 has to include dose l and time t ( in this order )
    # idem for mus_cd2 and mus_cd3
    # tau_cd1 and tau_cd2 need to include doses l
    
    if ( delta_2cdi < kappa1 + 1){
      return(delta_2cdi)
    }else{
      
      ones_T = rep(1, T)
      
      log_N_1cdi = rep(0, kappa1)
      
      for ( m in 1:kappa1){
        
        for (l in 1:L){
          
          mu_cdli_m = build_mu_cdli_margin( mus_cdl1 = mus_cd1[l, ], 
                                            mus_cdl2 = mus_cd2[l, ],
                                            mus_cdl3 = mus_cd3[l, ],
                                            tau_cdl1  = tau_cd1[l],
                                            tau_cdl2  = tau_cd2[l],
                                            delta_1cdi = m,
                                            delta_2cdi = delta_2cdi )
        
          for ( j in 1:J){
            expres = y_cdi[l,j, ] - mu_cdli_m  
            log_N_1cdi[m] = log_N_1cdi[m] - 0.5 * t( expres ) %*% Sigma_inv_c %*% expres
          }
          
        }
        
      }
      
      prob_delta_1cdi = exp( log_N_1cdi + log(pi1) - max( log_N_1cdi + log(pi1) ) )
      prob_delta_1cdi = prob_delta_1cdi/ sum(prob_delta_1cdi)
      
      new_delta1 = sample( size = 1, x = 1:kappa1, prob = prob_delta_1cdi )
      
      return( new_delta1 )
      
    }
    
  }
  
  ##########################
  #  sampling delta_2cdi
  ##########################
  
  sample_delta_2cdi = function( kappa1, kappa2,
                                mus_cd1, mus_cd2, mus_cd3, 
                                tau_cd1, tau_cd2, delta_1cdi, 
                                Sigma_inv_c, y_cdi, gamma, pi2){
    
    # mus_cd1 has to include dose l and time t ( in this order )
    # idem for mus_cd2 and mus_cd3
    # tau_cd1 and tau_cd2 need to include doses l
    
    ones_T = rep(1, T)
    
    log_N_2cdi = rep(0, kappa2 )
    prob_delta_2cdi = rep(0, kappa2)
    
    for ( m in c( (kappa1 + 1) : kappa2, delta_1cdi ) ){
      
      for (l in 1:L){
        
        mu_cdli_m  = build_mu_cdli_margin( mus_cdl1 = mus_cd1[l, ], 
                                          mus_cdl2 = mus_cd2[l, ],
                                          mus_cdl3 = mus_cd3[l, ],
                                          tau_cdl1 = tau_cd1[ l],
                                          tau_cdl2 = tau_cd2[ l],
                                          delta_1cdi = delta_1cdi,
                                          delta_2cdi = m )
      
      
        for ( j in 1:J){
          expres = y_cdi[l,j, ] - mu_cdli_m
          log_N_2cdi[m] = log_N_2cdi[m] - 0.5 * t( expres ) %*% Sigma_inv_c %*% expres
        }
      }
      
    }
    
    index = (kappa1 + 1): kappa2
    prob_delta_2cdi[ index ] = log_N_2cdi[ index ] + log(pi2) + log( 1 - gamma )
    
    index = delta_1cdi
    prob_delta_2cdi[ index ] = log_N_2cdi[ index ] + log( gamma )
    
    index = c( delta_1cdi, (kappa1 + 1): kappa2 )
    prob_delta_2cdi[ index ] = prob_delta_2cdi[ index ] - max( prob_delta_2cdi[ index ] )
    prob_delta_2cdi[ index ] = exp( prob_delta_2cdi[ index ] )/ sum ( exp( prob_delta_2cdi[ index ] ) )
    
    new_delta2 = sample( size = 1, x = 1:kappa2, prob = prob_delta_2cdi )
    
    return( new_delta2 )
    
  }
  
  
  ############################################################################
  ##########################      Updates      ###############################
  ############################################################################
  
  start_time = proc.time()
  
  for (iter in 1:n_iter){
    
    #######################
    #    Updating v_0u
    #######################
    
    #v_01
    u = 1 
    v0u[u] = rgamma( 1, 
                     shape = av + 0.5 * L*CC*D*kappa1, 
                     rate = 0.5 * sum( (mus[ , , , u, ] - mu0u[u])^2 ) + bv ) 
    
    #v_02                               
    u = 2 
    v0u[u] = rgamma( 1, 
                     shape = av + 0.5 * L*CC*D*kappa2, 
                     rate = 0.5 * sum( (mus[ , , , u, ] - mu0u[u])^2 ) + bv ) 
    
    
    #v_03
    u = 3 
    v0u[u] = rgamma( 1, 
                     shape = av + 0.5 * L*CC*D*kappa1, 
                     rate = 0.5 * sum( (mus[ , , , u, ] - mu0u[u])^2 ) + bv ) 
    

    #####################
    #   Updating mu_0u
    #####################
    
    # mu_01
    u = 1
    mu0u[u] = rnorm ( 1, 
                      sd = sqrt( 1/ (v00 + L*CC*D*kappa1 * v0u[u]) ), 
                      mean = ( v0u[u] * sum( mus[ , , , 1, 1:kappa1] ) + mu00*v00 ) / (v00 + L*CC*D*kappa1*v0u[u]) )
    
    # mu_02
    u = 2
    mu0u[u] = rnorm ( 1, 
                      sd = sqrt( 1/ (v00 + CC*D*kappa2 * v0u[u]) ), 
                      mean = ( v0u[u] *  sum( mus[ , , , 2, 1:kappa2] ) +   mu00*v00 ) / ( v00 + L*CC*D*kappa2* v0u[u] ) )
    
    # mu_03
    u = 3
    mu0u[u] = rnorm ( 1, 
                      sd = sqrt( 1/ (v00 + CC*D*kappa1 * v0u[u]) ), 
                      mean = ( v0u[u] * sum( mus[ , , , 3, 1:kappa1] ) +  mu00*v00 ) / (v00 + L*CC*D*kappa1 * v0u[u]) )
    
    ###################
    #  Updating tau1
    ###################
    
    for (cc in 1:CC){
      for (d in 1:D){
        for (l in 1:L){
        
          tau[cc, d, l, 1] = sampling_tau_cdl1( l = l, mu0u = mu0u, 
                                               v0u = v0u, 
                                               delta_1cd = delta[1,cc,d, ], 
                                               delta_2cd = delta[2,cc,d, ], 
                                               Sigma_inv_c = Sigma_inv[,,cc], 
                                               tau_cdl2 = tau[cc,d,l,2],
                                               kappa1 = kappa1, 
                                               kappa2 = kappa2, 
                                               sum_y_over_j_cd = sum_y_over_j[cc,d, , , ], 
                                               y_cd = y[cc,d, , , , ])
        }
      }
    }
    
    ##################
    #  Updating tau2
    ##################
    
    for (cc in 1:CC){
      for (d in 1:D){
        for (l in 1:L){
          
          tau[cc, d, l, 2] = sampling_tau_cdl2(l = l,
                                               mu0u = mu0u, 
                                               v0u = v0u, 
                                               delta_1cd = delta[1,cc,d, ], 
                                               delta_2cd = delta[2,cc,d, ], 
                                               Sigma_inv_c = Sigma_inv[,,cc], 
                                               tau_cdl1 = tau[cc,d,l,1],
                                               kappa1 = kappa1, 
                                               kappa2 = kappa2, 
                                               sum_y_over_j_cd = sum_y_over_j[cc,d, , , ], 
                                               y_cd = y[cc,d, , , , ])
        }
      }
    }
    
    ########################
    #   Updating mu*_cdlu
    ########################
    
    for (cc in 1:CC){
      
      Sigma_c_inv = Sigma_inv[,,cc] 
      
      for (d in 1:D){
        
        for (l in 1:L){
        
          tau_cdl1 = tau[cc,d,l,1] 
          tau_cdl2 = tau[cc,d,l,2]
          
          delta_1cd = delta[1,cc,d, ] 
          delta_2cd = delta[2,cc,d, ]
          
          y_cd = y[cc,d, , , , ]
          
          for (m in 1:kappa1 ){
            
            mus[cc,d,l,1,m] = sample_mus_cdl1( l=l, m=m, Sigma_c_inv = Sigma_c_inv, 
                                              mus_cdl1 = mus[cc,d,l,1, ], 
                                              mus_cdl2 = mus[cc,d,l,2, ], 
                                              mus_cdl3 = mus[cc,d,l,3, ], 
                                              tau_cdl1 = tau_cdl1, 
                                              tau_cdl2 = tau_cdl2,
                                              delta_1cd = delta_1cd, 
                                              delta_2cd = delta_2cd,
                                              y_cd = y_cd,
                                              kappa1 = kappa1 )
          }
          for (m in 1:kappa2){
            mus[cc,d,l,2,m] = sample_mus_cdl2( l=l, m=m, Sigma_c_inv = Sigma_c_inv, 
                                              mus_cdl1 = mus[cc,d,l,1, ], 
                                              mus_cdl3 = mus[cc,d,l,3, ], 
                                              tau_cdl1 = tau_cdl1, 
                                              tau_cdl2 = tau_cdl2,
                                              delta_1cd = delta_1cd, 
                                              delta_2cd = delta_2cd,
                                              y_cd = y_cd)
          }
          for (m in 1:kappa1){
            mus[cc,d,l,3,m] = sample_mus_cdl3( l=l, m=m, Sigma_c_inv = Sigma_c_inv, 
                                              mus_cdl1 = mus[cc,d,l,1, ], 
                                              mus_cdl2 = mus[cc,d,l,2, ], 
                                              tau_cdl1 = tau_cdl1, 
                                              tau_cdl2 = tau_cdl2,
                                              delta_1cd = delta_1cd, 
                                              delta_2cd = delta_2cd,
                                              y_cd = y_cd)
          }
          
        } # end for loop in l
      } # end for loop in d
    } # end for loop in cc
    
    
    #################################
    #  building the vectors mu_cdli
    #################################
    
    mu_cdli = array(0, dim = c(CC, D, L, n, T) )
    
    for (cc in 1:CC){
      for (d in 1:D){
        for (l in 1:L){
          
          mus_cdl1 = mus[cc,d,l,1, ]
          mus_cdl2 = mus[cc,d,l,2, ]
          mus_cdl3 = mus[cc,d,l,3, ]
          
          tau_cdl1 = tau[cc,d,l,1]
          tau_cdl2 = tau[cc,d,l,2]
          
          for (i in 1:n){
            mu_cdli[cc,d,l,i, ] = build_mu_cdli_margin(mus_cdl1 = mus_cdl1,
                                                      mus_cdl2 = mus_cdl2,
                                                      mus_cdl3 = mus_cdl3,
                                                      tau_cdl1 = tau_cdl1,
                                                      tau_cdl2 = tau_cdl2,
                                                      delta_1cdi = delta[1,cc,d,i],
                                                      delta_2cdi = delta[2,cc,d,i] )
          }
        }
      }
    }
    
    ######################
    #   Udating Sigma_c
    ######################
    
    for( cc in 1:CC ){
      
      scale_IW = V_Sigma
      
      for( d in 1:D){
        for (i in 1:n){
          for (j in 1:J){
            for (l in 1:L){
              scale_IW = scale_IW + ( y[cc, d, i, l, j, ] - mu_cdli[cc, d, l, i, ] ) %*% 
                t( y[cc, d, i, l, j, ] - mu_cdli[cc, d, l, i, ]  ) 
            }
          }
        }
      }
      
      Sigma[,,cc] = riwish(v = n*D*L*J + v_Sigma, S=scale_IW)
      Sigma_inv[,,cc] = solve( Sigma[,,cc] )
      
    }
    
    ######################
    #  Updating gamma
    ######################
    
    how_many_equal = sum( delta[ 2, , , ] == delta[ 1, , , ] )
    gamma = rbeta( 1, shape1 = a_gamma + how_many_equal, shape2 = b_gamma + n*CC*D - how_many_equal )
    

    ########################
    #  Updating delta_1cdi
    ########################
    
    for (cc in 1:CC){
      
      Sigma_inv_c = Sigma_inv[,,cc]
      
      for (d in 1:D){
        
        mus_cd1 = mus[cc,d, ,1, ] 
        mus_cd2 = mus[cc,d, ,2, ] 
        mus_cd3 = mus[cc,d, ,3, ]
        
        tau_cd1 = tau[cc,d, ,1]
        tau_cd2 = tau[cc,d, ,2]
        
        for (i in 1:n){
          delta[1, cc, d, i] = sample_delta_1cdi( kappa1 = kappa1, kappa2 = kappa2, 
                                                  mus_cd1 = mus_cd1, 
                                                  mus_cd2 = mus_cd2, 
                                                  mus_cd3 = mus_cd3, 
                                                  tau_cd1 = tau_cd1, 
                                                  tau_cd2 = tau_cd2,
                                                  delta_2cdi = delta[2, cc, d, i], 
                                                  Sigma_inv_c = Sigma_inv_c, 
                                                  y_cdi = y[cc,d,i, , , ], 
                                                  pi1=pi1 )
        }
      }
    }
    
    ########################
    #  Updating delta_2cdi
    ########################
    
    for (cc in 1:CC){
      
      Sigma_inv_c = Sigma_inv[,,cc]
      
      for (d in 1:D){
        
        mus_cd1 = mus[cc,d, ,1, ] 
        mus_cd2 = mus[cc,d, ,2, ] 
        mus_cd3 = mus[cc,d, ,3, ]
        
        tau_cd1 = tau[cc,d, ,1]
        tau_cd2 = tau[cc,d, ,2]
        
        for (i in 1:n){
          delta[2, cc, d, i] = sample_delta_2cdi( kappa1=kappa1, kappa2=kappa2, 
                                                  mus_cd1 = mus_cd1, 
                                                  mus_cd2 = mus_cd2,
                                                  mus_cd3 = mus_cd3,
                                                  tau_cd1 = tau_cd1, 
                                                  tau_cd2 = tau_cd2,
                                                  delta_1cdi = delta[1, cc, d, i], 
                                                  Sigma_inv_c = Sigma_inv_c, 
                                                  y_cdi = y[cc,d,i, , , ], 
                                                  pi2=pi2, 
                                                  gamma = gamma )
        }
      }
    }
    
    
    ###################
    #  Updating pi1
    ###################
    
    dir_param = rep(0, kappa1)
    
    for ( a in 1: kappa1){
      dir_param[a] = sum( delta[ 1, , , ] == a) + eta1[a]
    }
    
    pi1 = as.numeric( rdirichlet( 1, alpha = dir_param ) )
    
    ###################
    #  Updating pi2
    ###################
    
    dir_param = rep(0, kappa2 - kappa1)
    
    for ( a in 1: (kappa2 - kappa1) ){
      dir_param[a] = sum( delta[ 2, , , ] == a + kappa1 ) + eta2[a]
    }
    
    pi2 = as.numeric( rdirichlet( 1, alpha = dir_param ) )
    
    #################################
    #  building the vectors mu_cdli
    #################################
    
    mu_cdli = array(0, dim = c(CC, D, L, n, T) )
    
    for (cc in 1:CC){
      for (d in 1:D){
        for (l in 1:L){
        
          mus_cdl1 = mus[cc,d,l,1, ]
          mus_cdl2 = mus[cc,d,l,2, ]
          mus_cdl3 = mus[cc,d,l,3, ]
          
          tau_cdl1 = tau[cc,d,l,1]
          tau_cdl2 = tau[cc,d,l,2]
          
          for (i in 1:n){
            mu_cdli[cc,d,l,i, ] = build_mu_cdli_margin( mus_cdl1 =  mus_cdl1,
                                                        mus_cdl2 =  mus_cdl2,
                                                        mus_cdl3 =  mus_cdl3,
                                                        tau_cdl1 = tau_cdl1,
                                                        tau_cdl2 = tau_cdl2,
                                                        delta_1cdi = delta[1,cc,d,i],
                                                        delta_2cdi = delta[2,cc,d,i] )
          }
        }
      }
    }
    
    ############################
    #       Log likelihood
    ############################ 
    
    log_lik = 0
    
    for(cc in 1:CC){
      for (d in 1:D){
        for (l in 1:L){
          for(i in 1:n){
            for (j in 1:J){
              expres = ( y[cc,d,i,l,j, ] - mu_cdli[cc,d,l,i, ] )
              log_lik = log_lik - 0.5 * T * log( 2 * pi ) +
                        0.5 * log( det(Sigma_inv[ , , cc]) ) - 
                        0.5 * t(expres) %*% Sigma_inv[,,cc] %*% expres
            }
          }
        }
      }
    }
    
    #################################
    #   Storing the sampled values
    #################################
    
    chain_log_lik[iter] = log_lik
    
    chain_mu_cdli[ , , , , , iter ] = mu_cdli
    chain_mu_star[ , , , , , iter ] = mus
    
    chain_Sigma[ , , , iter] = Sigma
    chain_Sigma_inv[ , , , iter] = Sigma_inv
    chain_delta[ , , , , iter ] = delta
    chain_tau [ , , , , iter] = tau
    
    chain_gamma[ iter ] = gamma
    chain_pi2[ iter, ] = pi2
    chain_pi1[ iter, ] = pi1
    chain_v0u[ iter, ] = v0u
    chain_mu0u[ iter, ] = mu0u
    
    print(iter)
    
  }
  
  # index after convergence
  index = (burn_in+1):n_iter
  
  ####################
  # calculating BIC
  ####################
  BIC = max( chain_log_lik ) - n_param * log( n_data )/2 # the bigger the better
  
  ####################
  # calculating AIC
  ####################
  AIC = max( chain_log_lik ) - n_param # the bigger the better
  
  ####################
  # calculating WAIC
  ####################
  
  # log "predictive density"
  outside_sum = 0
  for(cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        for (i in 1:n){
          for (j in 1:J){
            inside_sum = 0
            for (iter in index){
              inside_sum = inside_sum + dmvnorm( x = y[cc, d, i, l ,j, ], 
                                                 mean = chain_mu_cdli[cc,d,l,i, , iter],
                                                 sigma = chain_Sigma[ , , cc,iter]) 
            }
            outside_sum = outside_sum + log( max( inside_sum/length(index), 10^(-20) ) )
          }
        }
      }
    }
  }
  log_predictive = outside_sum
  outside_sum = 0

  pw1 = -2*mean( chain_log_lik[ index ] ) + 2 * log_predictive
  WAIC = log_predictive - pw1 # the bigger the better
  
  ###########################
  #    marginal likelihood
  ###########################
  
  max_log_lik = max( chain_log_lik[index] )
  margin_lik = 1 / ( mean( exp( max_log_lik - chain_log_lik[index] ) ) / exp(max_log_lik) )
  
  ###########################
  #          LPML
  ###########################
  
  outside_sum = 0
  for(cc in 1:CC){
    for (d in 1:D){
      for (l in 1:L){
        for (i in 1:n){
          for (j in 1:J){
            inside_sum = 0
            for (iter in index){
              inside_sum = inside_sum + 1/max( dmvnorm( x = y[cc, d, i, l ,j, ], 
                                                        mean = chain_mu_cdli[cc,d,l,i, , iter],
                                                        sigma = chain_Sigma[ , , cc,iter]), 10^(-20) )
            }
            outside_sum = outside_sum + log( inside_sum/length(index) )
          }
        }
      }
    }
  }
  
  lpml = -1 * outside_sum
  
  
  #####################################
  # cluster membership point estimate
  #####################################
  
  # loading cpp function for point cluster estimation
  # server
  sourceCpp("./functions/choose_cluster.cpp")
  #local
  #sourceCpp("./choose_cluster.cpp")
  
  #load("./mcmc_rppa_ver10.Rdata")
  
  # delta 1
  ind1 = matrix(0, nrow = CC, ncol = D)
  est_delta1 = array(0, dim=c( CC, D, n ) )
  for(cc in 1:CC){
    for (d in 1:D){
      res = as.numeric( choose_cluster( chain_delta[1, cc, d, , ] ) )
      min_res = min(res)
      ind1[cc, d] = max( which( res == min_res ) )
      est_delta1[cc, d, ] = chain_delta[1,cc,d, ,ind1[cc,d] ]
    }
  }
  
  # delta 2
  ind2 = matrix(0, nrow = CC, ncol = D)
  est_delta2 = array(0, dim=c( CC, D, n ) )
  for(cc in 1:CC){
    for (d in 1:D){
      res = as.numeric( choose_cluster( chain_delta[2, cc, d, , ] ) )
      min_res = min(res)
      ind2[cc, d] = max( which( res == min_res ) )
      est_delta2[cc, d, ] = chain_delta[2,cc,d, ,ind1[cc,d] ]
    }
  }
  
  # time elapsed
  end_time = proc.time()
  print ( end_time - start_time )
  
  return( list( "mu" = chain_mu_cdli,
                "mu_star" = chain_mu_star,
                "Sigma" = chain_Sigma,
                "delta" = delta,
                "tau" = chain_tau,
                "gamma" = chain_gamma,
                "pi1" = chain_pi1,
                "pi2" = chain_pi2,
                "v0" = chain_v0u,
                "mu0" = chain_mu0u,
                "log_lik" = chain_log_lik,
                "AIC" = AIC,
                "BIC" = BIC,
                "WAIC" = WAIC,
                "est_delta1" = est_delta1,
                "est_delta2" = est_delta2,
                "margin_lik" = margin_lik,
                "lpml" = lpml,
                "kappa1" = kappa1,
                "kappa2" = kappa2) )
  
}

# Saving the mcmc function
save( mcmc_rppa, file ="/home/tadeu/MEGA/RPPA/depRPPA/inference/functions/mcmc_rppa_ver10.Rdata")

# Running the mcmc
load("/home/tadeu/ut_austin/RPPA/depRPPA/inference/real_data/y_sim3.Rdata")

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


chain_46 = mcmc_rppa(n_iter=5, kappa1=4, kappa2=6, y=y, burn_in = 2)
save( chain_46, file = "/home/tadeu/MEGA/RPPA/depRPPA/inference/real_data/mcmc_chains/best_chain_64.Rdata" )


###############
#  kappa1 = 2
###############

chain = mcmc_rppa(n_iter=500, kappa1=2, kappa2=3, y=y, burn_in = 100)
chain_23 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1


chain = mcmc_rppa(n_iter=500, kappa1=2, kappa2=4, y=y, burn_in = 100)
chain_24 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=2, kappa2=5, y=y, burn_in = 100)
chain_25 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

###############
#  kappa1 = 3
###############

chain = mcmc_rppa(n_iter=500, kappa1=3, kappa2=4, y=y, burn_in = 100)
chain_34 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=3, kappa2=5, y=y, burn_in = 100)
chain_35 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=3, kappa2=6, y=y, burn_in = 100)
chain_36 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

###############
#  kappa1 = 4
###############

chain = mcmc_rppa(n_iter=500, kappa1=4, kappa2=5, y=y, burn_in = 100)
chain_45 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=4, kappa2=6, y=y, burn_in = 100)
chain_46 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=4, kappa2=7, y=y, burn_in = 100)
chain_47 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1


###############
#  kappa1 = 5
###############

chain = mcmc_rppa(n_iter=500, kappa1=5, kappa2=6, y=y, burn_in = 100)
chain_56 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=5, kappa2=7, y=y, burn_in = 100)
chain_57 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=5, kappa2=8, y=y, burn_in = 100)
chain_58 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1


###############
#  kappa1 = 6
###############

chain = mcmc_rppa(n_iter=500, kappa1=6, kappa2=7, y=y, burn_in = 100)
chain_67 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=6, kappa2=8, y=y, burn_in = 100)
chain_68 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=6, kappa2=9, y=y, burn_in = 100)
chain_69 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1


###############
#  kappa1 = 7
###############

chain = mcmc_rppa(n_iter=500, kappa1=7, kappa2=8, y=y, burn_in = 100)
chain_78 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=7, kappa2=9, y=y, burn_in = 100)
chain_79 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=7, kappa2=10, y=y, burn_in = 100)
chain_710 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1


###############
#  kappa1 = 8
###############

chain = mcmc_rppa(n_iter=500, kappa1=8, kappa2=9, y=y, burn_in = 100)
chain_89 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=8, kappa2=10, y=y, burn_in = 100)
chain_810 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1

chain = mcmc_rppa(n_iter=500, kappa1=8, kappa2=11, y=y, burn_in = 100)
chain_811 = chain
aic_vec[ind] = chain$AIC
bic_vec[ind] = chain$BIC
dic_vec[ind] = calc_dic( chain=chain, y=y, init_iter=101, end_iter=500)
waic_vec[ind] = chain$WAIC
margin_lik_vec[ind] = chain$margin_lik
lpml_vec[ind] = chain$lpml
print(chain$AIC); print(chain$BIC); print(chain$WAIC); print( dic_vec[ind] ); print(chain$margin_lik); print(chain$lpml) 
ind = ind + 1



model_comparison = rbind( round( rbind( aic_vec, bic_vec, dic_vec, waic_vec), 3), margin_lik_vec, lpml_vec)

# saving results
write.table(model_comparison, row.names = TRUE, col.names = FALSE, file = "/home/tadeu/ut_austin/RPPA/depRPPA/writing/model_comparison.txt")

# kappa1 = 2
save( chain_23, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_23_comparison.Rdata" )
save( chain_24, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_24_comparison.Rdata" )
save( chain_25, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_25_comparison.Rdata" )

# kappa1 = 3
save( chain_34, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_34_comparison.Rdata" )
save( chain_35, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_35_comparison.Rdata" )
save( chain_36, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_36_comparison.Rdata" )

# kappa1 = 4
save( chain_45, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_45_comparison.Rdata" )
save( chain_46, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_46_comparison.Rdata" )
save( chain_47, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_47_comparison.Rdata" )

# kappa1 = 5
save( chain_56, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_56_comparison.Rdata" )
save( chain_57, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_57_comparison.Rdata" )
save( chain_58, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_58_comparison.Rdata" )


# kappa1 = 6
save( chain_67, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_67_comparison.Rdata" )
save( chain_68, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_68_comparison.Rdata" )
save( chain_69, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_69_comparison.Rdata" )

# kappa1 = 7
save( chain_78, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_78_comparison.Rdata" )
save( chain_79, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_79_comparison.Rdata" )
save( chain_710, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_710_comparison.Rdata" )

# kappa1 = 8
save( chain_89, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_89_comparison.Rdata" )
save( chain_810, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_810_comparison.Rdata" )
save( chain_811, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_811_comparison.Rdata" )


# save( chain, file = "/home/tadeu/ut_austin/RPPA/depRPPA/inference/chain_real_ver10_identify.Rdata" )

# Profile
# Rprof(NULL)
# summarize the profile
# summaryRprof("/home/tadeu/ut_austin/RPPA/depRPPA/inference/rppa_prof.txt")

