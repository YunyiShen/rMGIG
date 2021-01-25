// [[Rcpp::depends(RcppArmadillo)]]
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
#include <Riccati_Solver.h>
/* ARE_Solver_cpp solves Algebric Riccati Equations with the form
 - Continuous time
    A'*X + X*A + XRX+Q = 0, where X is symmetric.

  follow notations in Farideh Fazayeli and Arindam Banerjee
*/

// [[Rcpp::export]]
void CARE_ArimotoPotter_cpp(arma::mat &X,
                            const arma::mat &A,
                            const arma::mat &Q,
                            const arma::mat &R)
{
    // Hamiltonian matrix:
    int n = A.n_rows;
    arma::mat Z(2 * n, 2 * n, arma::fill::zeros);
    Z.submat(0, 0, n - 1, n - 1) = A;
    Z.submat(0, n, n - 1, 2 * n - 1) = R;
    Z.submat(n, 0, 2 * n - 1, n - 1) = -Q;
    Z.submat(n, n, 2 * n - 1, 2 * n - 1) = -A.t();
    //Rcout << Z <<endl;
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    eig_gen(eigval, eigvec, Z);
    arma::vec real_eigval = arma::real(eigval);
    //Rcout << eigval <<endl;
    // stable eigen vectors, i.e. with negative eigen values
    arma::cx_mat stable_eigenvec(2 * n, n, arma::fill::zeros);
    int j = 0;
    for (int i = 0; i < 2 * n; ++i)
    {
        if (real_eigval(i) < 0)
        {
            stable_eigenvec.col(j) = eigvec.col(i);
            ++j;
        }
    }
    //Rcout << stable_eigenvec << endl;
    arma::cx_mat U1 = stable_eigenvec.rows(0, n - 1);
    arma::cx_mat U2 = stable_eigenvec.rows(n, 2 * n - 1);
    X = arma::real(U2 * arma::inv(U1));
    return;
}
