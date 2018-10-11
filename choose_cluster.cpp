
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// This script contains the C++ function used to estimate the cluster membership indicators based on
// a given input. The input is a matrix of dimension nxm with m representing the number of MCMC iterations
// and n the number of units being clustered

// [[Rcpp::export]]
arma::colvec choose_cluster( arma::mat delta ) {
  
  int n = delta.n_rows;
  int m = delta.n_cols;
  arma::mat phat(n, n); phat.fill(0.0);
  double aux=0;
  
  for( int i=1; i<n; i++ ){
    for( int j=0; j<i; j++ ){
      for ( int k=0; k<m; k++ ){
       if( delta(i, k) == delta(j, k) ){ 
          aux = aux + 1; 
         }
       }
    phat(i, j) = aux/m;
    aux = 0;
    }
  }
  
  arma::mat matrix_10(n,n); 
  arma::colvec distances(m);
  
  for (int iter=0; iter<m; iter++){
    matrix_10.fill(0.0);
    for( int i=1; i<n; i++ ){
      for( int j=0; j<i; j++ ){
        if( delta(i, iter) == delta(j, iter) ){
          matrix_10(i,j) = 1;
        }
      }
    }
  
    double aux2 = 0;
  
    for( int i=1; i<n; i++ ){
      for( int j=0; j<i; j++ ){
        aux2 = aux2 + std::abs( matrix_10(i,j) - phat(i,j) ) ;
      }
    }
    distances(iter) = aux2;
  }  
    
  return distances;
}

