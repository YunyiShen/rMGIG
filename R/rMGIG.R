get_mat <- function(i,mat,n){
    ut <- mat[,i]
    temp <- matrix(0,n,n)
    temp[upper.tri(temp,T)] <- ut
    temp <- temp + t(temp)
    diag(temp) <- diag(temp)/2
    return(temp)
}

is.positive.definite <- function(mat,tol = 1e-8){
    mat[abs(mat)<=tol] <- 0
    all(eigen(mat,T,T,F)$values>=0)
}



#' Sampling from Matrix Generalized Inverse Gaussian distribution
#'
#' Method follow Farideh Fazayeli and Arindam Banerjee 2016
#' "The Matrix Generalized Inverse Gaussian Distribution:
#' Properties and Applications"
#' 
#'
#' @param n number of sample needed
#' @param nu double, nu parameter for the MGIG
#' @param phi phi parameter positive (semi)defined matrix, reduce to inverse-Wishart if set to 0
#' @param psi psi parameter positive (semi)defined matrix, reduce to Wishart if set to 0
#' @param df positive double, should be greater than dimision-1 degree of freedom for the proposal Wishart distribution during importance sampling
#' @param list bool, whether want list output, if ture will return list of full matrix otherwise return matrix of upper triagular part
#' @param maxit integer maximum iterations for importance sampling algorithm
#' @return a matrix with upper triangular part of the positive defined matrix as columns or a list whose entries are the sampled matrix
#' @export
#' @examples
#' my.sample <- rMGIG(n = 1, nu = 6, phi = diag(3), psi = diag(3), df = 6, list = FALSE, maxit = 1000)
#' my.sample.list <- rMGIG(n = 1, nu = 6, phi = diag(3), psi = diag(3), df = 6, list = TRUE, maxit = 1000)



rMGIG <- function(n = 1, nu = 6, phi = diag(3), 
    psi = diag(3), df = 6, list = F, maxit = 10000) {
    if (nrow(phi) != ncol(phi))
        stop("phi must be square")
    if (nrow(psi) != ncol(psi))
        stop("psi must be square")
    if (nrow(psi) != nrow(phi))
        stop("dimension of phi and psi must match")
    if(!isSymmetric(phi))
        stop("phi must be symmetric")
    if(!isSymmetric(psi))
        stop("psi must be symmetric")
    if (!is.positive.definite(phi, tol = 1e-8))
        stop("phi must be positive defined")
    if (!is.positive.definite(psi, tol = 1e-8))
        stop("psi must be positive defined")
    if (df <= nrow(phi))
        stop("df must be greater than dimension-1")

    res <- rMGIG_cpp(n,nu,phi,psi,df,maxit)

    if(sum(is.na(res))>0) 
        warning("failed to reach desired sample size, consider larger maxiter, or change df")
    
    if(list){
        return(lapply(1:n,get_mat,res,nrow(psi)))
    }
    else {
       return(res)
    }
}



