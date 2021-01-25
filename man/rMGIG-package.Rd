\name{rMGIG-package}
\alias{rMGIG-package}
\docType{package}
\title{
\packageTitle{rMGIG}
}
\description{
\packageDescription{rMGIG}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rMGIG}
\packageIndices{rMGIG}
`rMGIG` draw sample from the matrix Generalized Inverse Gaussian (MGIG) Distribution
}
\author{
\packageAuthor{rMGIG}

Maintainer: \packageMaintainer{rMGIG}
}
\references{
    Fazayeli, Farideh, and Arindam Banerjee. "The matrix generalized inverse Gaussian distribution: Properties and applications." Joint European Conference on Machine Learning and Knowledge Discovery in Databases. Springer, Cham, 2016.
}

\keyword{ package }
\examples{
my.sample <- rMGIG(n = 1, nu = 6, phi = diag(3), psi = diag(3), df = 6, list = FALSE, maxit = 1000)
}
