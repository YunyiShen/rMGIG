#ifndef RICCATI_SOLVER_H
#define RICCATI_SOLVER_H
// [[Rcpp::depends(RcppArmadillo)]]
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
/* ARE_Solver_cpp solves Algebric Riccati Equations with the form
 - Continuous time
    A'*P + P*A - PBR^-1 B'P+Q = 0, where P is symmetric.
*/

// [[Rcpp::export]]
void CARE_ArimotoPotter_cpp(arma::mat &X,
                            const arma::mat &A,
                            const arma::mat &R,
                            const arma::mat &Q);

#endif