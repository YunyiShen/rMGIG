// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#include "Riccati_Solver.h"

/*
Rejection sampling to sample MGIG distribution following
    Farideh Fazayeli and Arindam Banerjee 2016
    "The Matrix Generalized Inverse Gaussian Distribution: 
        Properties and Applications"
    
*/
// [[Rcpp::export]]
arma::mat rMGIG_cpp(int n, const double & nu, 
                    const arma::mat & phi, 
                    const arma::mat & psi,
                    const double & df,
                    int maxit){

    int N = phi.n_rows;
    // will save COLUMN vectors that are upper triangular entries of sample
    arma::mat res(0.5 * N * (N+1),n,arma::fill::zeros);
    res += NA_REAL; 

    arma::mat A = eye(size(phi));
    A *= (nu-(N+1)/2);
    //Rcout << A <<endl;

    // find mode of the MGIG by solving the CARE:
    arma::mat Lambda_star;
    CARE_ArimotoPotter_cpp(Lambda_star, A, phi, psi);
    
    // proposal distribution
    arma::mat Sigma_star = Lambda_star/(df-N-1);
    arma::mat chol_Sigma_star = chol(Sigma_star);
    arma::mat invSigma_star = inv(Sigma_star);
    arma::mat prop, samp;
    double detprop,signprop, logweight;
    int i_sample = 0; 
    for(int i = 0 ; i<maxit ; ++i){
        // draw proposal sample
        prop = arma::wishrnd( Sigma_star, df, chol_Sigma_star );
        arma::log_det(detprop,signprop,prop); 
        logweight = (nu-df/2) * detprop - 
            0.5 * arma::trace(psi * inv(prop)+(phi-invSigma_star)*prop); 
        if(log(R::runif(0,1))<=logweight){
            //Rcout << "flag" <<endl;
            res.col(i_sample) = prop(trimatu_ind(size(prop)));
            ++i_sample;
        }
        if(i_sample==n) break;
    }
    return(res);
}

// [[Rcpp::export]]
arma::mat mMGIG_cpp(const double & nu, 
                    const arma::mat & phi, 
                    const arma::mat & psi){

    int N = phi.n_rows;
    arma::mat A = eye(size(phi));
    A *= (nu-(N+1)/2);

    // find mode of the MGIG by solving the CARE:
    arma::mat Lambda_star;
    CARE_ArimotoPotter_cpp(Lambda_star, A, phi, psi);
    return(Lambda_star);
}

// [[Rcpp::export]]
double fMGIG_cpp(const arma::mat & X, 
                 const double & nu, 
                 const arma::mat & phi, 
                 const arma::mat & psi, 
                 bool logrithm = true){
    int d = X.n_cols;
    arma::mat X_inv = inv_sympd(X);
    double logf = (nu-(d+1.0)/2) * log_det_sympd(X);
    logf -= .5*trace(psi * X_inv + phi * X);
    return logrithm ? logf : exp(logf);  
}
