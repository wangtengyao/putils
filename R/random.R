########### Random element generation ##########

#' Generate n rademacher random variables
#' @param n length of random vector
#' @export
random.rademacher <- function(n=1){
    sample(c(-1,1),n,replace=T)
}

#' Generate n random Bernoulli variables
#' @param n length of random vector
#' @param p probability (can be a vector)
#' @export
random.bernoulli <- function(n=1, p=0.5) {
    rbinom(n,1,p)
}

#' Generate a random unit vectors in R^n
#' @param n length of random vector
#' @export
random.UnitVector <- function(n){
    v = rnorm(n)
    v/vector.norm(v)
}

#' Generate a Gaussian random matrix
#' @param nrow number of rows
#' @param ncol number of columns
#' @export
random.GaussianMatrix <- function(nrow, ncol=nrow){
    matrix(rnorm(nrow*ncol),nrow,ncol)
}

#' Generate a random nxn orthogonal matrix
#' @param n dimension
#' @export
random.OrthogonalMatrix <- function(n){
    A = matrix(rnorm(n*n),n,n)
    qr.Q(qr(A))
}

#' Generate a random nxn Wishart matrix
#' @param n dimension
#' @export
random.WishartMatrix <- function(n){
    A = matrix(rnorm(n*n),n,n)
    t(A)%*%A
}

#' Generate a random nxn Wigner matrix, where off diagonals are N(0,1) and diagonals are N(0,2)
#' @param n dimension
#' @export
random.WignerMatrix <- function(n){
    M = matrix(0,n,n)
    M[upper.tri(M, diag = TRUE)] = rnorm(n*(n+1)/2)
    M + t(M)
}

#' Generate a random positive semidefinite matrix with operator norm <=1
#' @param n dimension
#' @export
random.psdMatrix <- function(n){
    Lambda = diag(sort(runif(n),decreasing = TRUE))
    V = random.OrthogonalMatrix(n)
    V%*%Lambda%*%t(V)
}

#' Generate a random symmetric matrix
#' @param n dimension
#' @export
random.SymmetricMatrix <- function(n){
    Lambda = diag(sort(runif(n,min=-1,max=1),decreasing = TRUE))
    V = random.OrthogonalMatrix(n)
    V%*%Lambda%*%t(V)
}
