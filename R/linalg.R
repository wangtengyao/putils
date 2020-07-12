###########   Matrices and vectors   ##############

#' Norm of a vector
#' @description Calculate the entrywise L_q norm of a vector or a matrix
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return the entrywise L_q norm of a vector or a matrix
#' @export
vector.norm <- function(v, q = 2, na.rm = FALSE){
    if (na.rm) v <- na.omit(v)
    M <- max(abs(v))
    if (M == 0) return(0) else v <- v/M
    if (q == Inf) {
      nm <- max(abs(v))
    } else if (q > 0) {
      nm <- (sum(abs(v)^q))^(1/q)
    } else if (q == 0) {
      nm <- sum(v!=0)
    } else {
      return(NaN)
    }
    return(nm * M)
}

#' Norm of a matrix
#' @param M real matrix
#' @param norm name of matrix norm
#' @param type type of matrix norm
#' @param pnorms indexed by p in the same type
#' @param q norms indexed by q in the same type
#' @export
matrix.norm <- function(M, type=c('operator', 'entrywise', 'Schatten', 'rowcol', 'colrow'), p=2, q=p){
  if (length(type) > 1) type <- type[1]
  if (type == 'operator'){
    if (p==2 && q==2) {
      type<-'Schatten'; p<-Inf
    } else if (p==1) {
      type<-'colrow'
    } else if (q==Inf){
      type<-'rowcol'; p<-1/(1-1/p)
    } else {
      print('I do not know how to compute this norm.')
      return(NA)
    }
  }
  if (type == 'Frobenius') {type <- 'entrywise'; p <- 2}
  if (type == 'nuclear') {type <- 'Schatten'; p <- 1}

  if (type == 'entrywise') return(vector.norm(M, p))
  if (type == 'Schatten') return(vector.norm(svd(M)$d, p))
  if (type == 'rowcol') return(vector.norm(apply(M, 1, function(v) vector.norm(v,p)), q))
  if (type == 'colrow') return(vector.norm(apply(M, 2, function(v) vector.norm(v,p)), q))

  stop('norm type not recognised')
}


#' Normalise a vector
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return normalised version of this vector
#' @export
vector.normalise <- function(v, q = 2, na.rm = FALSE){
  return(v / vector.norm(v, q, na.rm))
}

#' Clipping a vector from above and below
#' @description Clipping vector or matrix x from above and below
#' @param x a vector of real numbers
#' @param upper clip above this value
#' @param lower clip below this value
#' @return the entrywise L_q norm of a vector or a matrix
#' @export
vector.clip <- function(x, upper = Inf, lower = -upper){
  if (upper < lower)  stop("upper limit cannot be below lower limit")
  x[x<lower]<-lower;
  x[x>upper]<-upper;
  x
}

#' Soft thresholding a vector
#' @param x a vector of real numbers
#' @param lambda soft thresholding value
#' @return a vector of the same length
#' @description entries of v are moved towards 0 by the amount lambda until they hit 0.
#' @export
vector.soft.thresh <- function(x, lambda){
  sign(x)*pmax(0,(abs(x)-lambda))
}

#' Hard thresholding a vector
#' @param x a vector of real numbers
#' @param lambda hard thresholding value
#' @return a vector of the same length
#' @description entries of v that are below lambda are set to 0.
#' @export
vector.hard.thresh <- function(x, lambda){
  x[abs(x)<lambda] <- 0; x
}

#' Apply Gram Schmidt operation to a matrix
#' @param A a matrix of real values
#' @param tolerance tolerance level for 0
#' @return a matrix of the same size with orthonormal columns
#' @details the resulting matrix is uniquely identified by asserting the diagonal to be non-negative
#' @export
matrix.GramSchmidt <- function(A, tolerance = 1e-14){
  A <- as.matrix(A); n = dim(A)[1]; p = dim(A)[2];
  qrdecomp <- qr(A)
  v <- diag(qr.R(qrdecomp))
  B <- qr.Q(qrdecomp)
  B[,abs(v)<tolerance] <- 0
  w <- sign(diag(B)); w[w==0] <- 1;
  W <- matrix(rep(c(w, rep(1,p-length(w))), each = n), n, p)
  B <- B*W
  B
}

#' Randomly shuffle a vector
#' @param v a vector
#' @return a shuffled vector of the same length
#' @export
vector.scramble <- function(v){
  v[sample(length(v))]
}

#' Clean up floating point epsilon terms
#' @param v a vector
#' @param tolerance level below which we set entries to 0
#' @return a vector of the same length
#' @export
vector.cleanup <- function(v, tolerance = 1e-14){
  v[abs(v)<tolerance] <- 0
  v
}

#' Trace of a square matrix
#' @param M a square matrix of real numbers
#' @return a real number
#' @export
matrix.trace <- function(M) {
    M <- as.matrix(M)
    if (dim(M)[1] != dim(M)[2]) stop('M is not a square matrix.')
    sum(diag(M))
}

#' Rank of a matrix
#' @param A a matrix of real numbers
#' @param tolerance level below which we treat a singular value to be 0
#' @return an nonnegative integer
#' @export
matrix.rank <- function(A, tolerance = 1e-10){
    v <- diag(qr.R(qr(A)))
    sum(abs(v)>tolerance)
}

#' Sine angle loss
#' @param U an orthonormal matrix
#' @param V an orthonormal matrix
#' @return sine angle loss between Vhat and V (Frobenius norm)
#' @export
sinThetaLoss <- function(U,V){
  sqrt(sum((U%*%t(U) - V%*%t(V))^2)/2)
}


#' Power method for leading eigenvectors and eigenvalues
#' @param A a square matrix
#' @param k number of leading eigenvectors/eigenvalues
#' @param eps tolerance for convergence (in Frobenius norm)
#' @param maxiter maximum iteration
#' @export
powerMethod <- function(A, k = 1, eps = 1e-10, maxiter = 1000){
  if (nrow(A) != ncol(A)) stop('powerMethod requires a square matrix')
  d <- nrow(A)
  V <- random.OrthogonalMatrix(d)[,seq_len(k)]
  for (i in seq_len(maxiter)){
    V_new <- matrix.GramSchmidt(A%*%V)
    if ((vector.norm(V_new - V)) < eps) break
    V <- V_new
    printPercentage(i, maxiter)
  }
  if (i == maxiter) warning('max iter reached without convergence.')
  evals <- diag(t(V_new)%*%A%*%V_new)
  return(list(values = evals, vectors = V_new))
}

#' Standardise the columns of a matrix to have mean zero and length 1
#' @param X a matrix
#' @return a standardised matrix
#' @export
matrix.standardise <- function(X){
  X <- sweep(X, 2, colMeans(X))
  X <- apply(X, 2, vector.normalise)
  return(X)
}


#' QR decomposition of a matrix
#' @param A an nxp matrix for QR decomposition
#' @return a list of two matrices Q and R so that A = QR
#' \itemize{
#'   \item Q - nxp matrix with orthonormal columns with the same span as A
#'   \item R - a upper triangular pxp matrix with nonnegative diagonal entries
#' }
#' @export
QR <- function(A){
  tmp <- qr(A)
  Q <- qr.Q(tmp)
  R <- qr.R(tmp)
  sign_flips <- sign(diag(R))
  sign_flips[sign_flips == 0] <- 1
  Q <- Q %*% diag(sign_flips)
  R <- diag(sign_flips) %*% R
  return(list(Q=Q, R=R))
}

#' QL decomposition of a matrix
#' @param A an nxp matrix for QL decomposition
#' @return a list of two matrices Q and L so that A = QL
#' \itemize{
#'   \item Q - nxp matrix with orthonormal columns with the same span as A
#'   \item L - a lower triangular pxp matrix with nonnegative diagonal entries
#' }
#' @export
QL <- function(A){
  B <- A[,ncol(A):1]
  tmp <- QR(B)
  Q <- tmp$Q
  R <- tmp$R
  Q <- Q[,ncol(Q):1]
  L <- R[nrow(R):1,ncol(R):1]
  return(list(Q=Q, L=L))
}

#' Extracting the off-diagonal part of A
#' @param A a rectangular matrix
#' @param as.vector if FALSE, a matrix with main diagonal set to 0 is returned
#' otherwise, a vector of offdiagonal entries is returned
#' @export
offdiag <- function(A, as.vector=FALSE){
  mask <- matrix(TRUE, nrow=nrow(A), ncol=ncol(A))
  diag(mask) <- FALSE
  if (as.vector){
    return(A[mask])
  } else {
    A[!mask] <- 0
    return(A)
  }
}

#' Compute matrix power of a symmetric matrix
#' @param A a square matrix
#' @param power power exponent
#' @param pseudoinverse whether to use pseudoinverse if power is negative
#' @export
matrix.power <- function(A, power, pseudoinverse=TRUE){
  if (nrow(A)!=ncol(A)) stop('A need to be a square matrix.')
  symm <- isSymmetric(A)
  tmp <- eigen(A, symmetric=symm)
  evecs <- tmp$vectors
  evals <- tmp$values
  evals[abs(evals) < 1e-12] <- 0

  if (sum(evals==0) > 0 && power < 0){
    if (!pseudoinverse){
      stop()
    } else {
      power <- -power
      evals[evals!=0] <- 1/evals[evals!=0]
    }
  }

  if (symm && min(Re(evals)) >=0) {
    return(evecs %*% diag(evals^power, nrow=nrow(A)) %*% t(evecs))
  } else {
    return(evecs %*% diag(as.complex(evals)^power, nrow=nrow(A)) %*% t(evecs))
  }
}

#' Compute Moore-Penrose pseudoinverse of a real rectangular matrix
#' @param M real matrix to be inverted
#' @param tolerance singular values below this level are treated as 0
#' @return matrix M.inv satisfying M %*% M.inv %*% M = M and M.inv %*% M %*% M.inv = M.inv
#' @export
pseudoinverse <- function(M, tolerance = 1e-10){
  tmp <- svd(M)
  u <- tmp$u; d <- tmp$d; v <- tmp$v
  filter <- abs(d) < tolerance
  d[filter] <- 0; d[!filter] <- 1/d[!filter]
  return(v%*%diag(d,nrow=length(d))%*%t(u))
}

#' Column maximum
#' @param M matrix
#' @return a vector of column maxima
#' @export
colMax <- function(M){
  apply(M, 2, max)
}

#' Column minimum
#' @param M matrix
#' @return a vector of column maxima
#' @export
colMin <- function(M){
  apply(M, 2, min)
}

#' Row maximum
#' @param M matrix
#' @return a vector of column maxima
#' @export
rowMax <- function(M){
  apply(M, 1, max)
}

#' Row minimum
#' @param M matrix
#' @return a vector of column maxima
#' @export
rowMin <- function(M){
  apply(M, 1, min)
}
