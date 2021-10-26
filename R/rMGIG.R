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
#' Sampling from Matrix Generalized Inverse Gaussian distribution, 
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
#' my.sample <- rMGIG(n = 1, nu = 6, phi = diag(3), 
#' psi = diag(3), df = 6, 
#' list = FALSE, maxit = 1000)
#' my.sample.list <- rMGIG(n = 1, nu = 6, phi = diag(3), 
#' psi = diag(3), df = 6, 
#' list = TRUE, maxit = 1000)



rMGIG <- function(n = 1, nu = 6, phi = diag(3), 
    psi = diag(3), df = 6, list = TRUE, maxit = 10000) {
    if (n<1 | (n%%1 != 0)){
        stop("n must be a positive integer")
    }
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
    if (df < nrow(phi)-1)
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





#' Return the mode of Matrix Generalized Inverse Gaussian distribution
#'
#' Mode of Matrix Generalized Inverse Gaussian distribution, 
#' Method follow Farideh Fazayeli and Arindam Banerjee 2016
#' "The Matrix Generalized Inverse Gaussian Distribution:
#' Properties and Applications" 
#' by solving the ARE 
#' 
#'
#' @param nu double, nu parameter for the MGIG
#' @param phi phi parameter positive (semi)defined matrix, reduce to inverse-Wishart if set to 0
#' @param psi psi parameter positive (semi)defined matrix, reduce to Wishart if set to 0
#' @return a matrix, the mode of the MGIG
#' @export
#' @examples
#' my.mode <- mMGIG(nu = 6, phi = diag(3), psi = diag(3))


mMGIG <- function(nu = 6, phi = diag(3), psi = diag(3)) {
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


    res <- mMGIG_cpp(nu,phi,psi)


    return(res)
}


#' Return the unnormalized density of MGIG
#'
#' The unnormalized density of MGIG given X
#' 
#' @param X a positive definite matrix where to calculate density
#' @param nu double, nu parameter for the MGIG
#' @param phi phi parameter positive (semi)defined matrix, reduce to inverse-Wishart if set to 0
#' @param psi psi parameter positive (semi)defined matrix, reduce to Wishart if set to 0
#' @param logirithm whether to take log
#' @return a matrix, the mode of the MGIG
#' @export
#' @examples
#' my.mode <- fMGIG(X = diag(3), nu = 6, phi = diag(3), psi = diag(3))


fMGIG <- function(X, nu = 6, phi = diag(3), psi = diag(3), logirithm = TRUE) {
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
    
    if(!isSymmetric(X))
        stop("X")
    if (!is.positive.definite(X, tol = 1e-8))
        stop("X must be positive defined")


    res <- fMGIG_cpp(X, nu,phi,psi, logirithm)
    return(res)
}


#' Estimate the KL divergence between two MGIG distribution
#'
#' Calculate the KL(MGIG(nu1,phi1,psi1)||MGIG(nu2,phi2,psi2))
#' 
#'
#' @param nu1 double, nu parameter for the first MGIG
#' @param phi1 phi parameter for the first positive (semi)defined matrix, reduce to inverse-Wishart if set to 0
#' @param psi1 psi parameter for the first positive (semi)defined matrix, reduce to Wishart if set to 0
#' @param nu2 double, nu parameter for the second MGIG
#' @param phi2 phi parameter for the second positive (semi)defined matrix, reduce to inverse-Wishart if set to 0
#' @param psi2 psi parameter for the second positive (semi)defined matrix, reduce to Wishart if set to 0
#' @param n_samples samples to take to estimate the KL divergence
#' @param df degree of freedom for the proposed Wishart when sampling the MGIG
#' @param maxit maximum iteration to sample MGIG
#' @return a matrix, the mode of the MGIG
#' @export

KLdiv <- function(nu1, phi1, psi1, nu2, phi2, psi2, n_samples = 5000 , df = 10*nrow(psi1), maxit = 1e6){
    samples <- rMGIG(n_samples, nu1, phi1, psi1, df, list = TRUE, maxit = maxit)
    logfoverg <- lapply(samples, function(X,nu1, phi1, psi1, nu2, phi2, psi2 ){
        fMGIG(X,nu1, phi1, psi1) - fMGIG(X,nu2, phi2, psi2)
    },nu1, phi1, psi1, nu2, phi2, psi2)

    term1 <- log(mean(sapply(logfoverg, function(x){exp(-x)})))
    term2 <- Reduce(mean, logfoverg)
    
    term1+term2
}
