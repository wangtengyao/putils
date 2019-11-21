#' putils: A package for personal utility functions.
#'
#' The putil package provides personal utility functions.
#' 
#' @section mathematical functsions:
#' * expit
#' * logit
#' * lambertW
#' 
#' @section matrices and vectors:
#' * vector.norm
#' * vector.normalise
#' * vector.clip
#' * vector.soft.thresh
#' * vector.hard.thresh
#' * vector.scramble
#' * vector.cleanup
#' * matrix.GramSchmidt
#' * matrix.trace
#' * matrix.rank
#' * matrix.standardise
#' * matrix.power
#' * sinThetaLoss
#' * powerMethod
#' * ql
#' * offdiag
#'
#' @section random element generation:
#' * random.rademacher
#' * random.bernoulli
#' * random.UnitVector
#' * random.OrthogonalMatrix
#' * random.WishartMatrix
#' * random.WignerMatrix
#' * random.psdMatrix
#' * random.SymmetricMatrix

#' @section auxiliary functions:
#' * printPercentage
#' * visualise
#' * snippet
#' * find.first
#' * find.last
#' * sf
#' * sf_exp
#' * sim.params
#' * show.params
#' * '%=%'

#' @section statistical functions:
#' * CvM.test
#' 
#' @section string operations:
#' * strlen
#' * strstr
#' * prefix
#' * suffix
#' 
#' @section NA handling:
#' * setNA
#' 
#' @docType package
#' @name putils
#' 
#' @importFrom graphics image
#' @importFrom stats na.omit rbinom rnorm runif
NULL

###########   Mathematical functions   ##############

#' Inverse of the logistic function
#' @param a a real number 
#' @return a real number of value 1/(exp(-a)+1)
#' @export
#' @seealso \code{\link{logit}}
expit <- function(a){ 1/(exp(-a)+1) }


#' The logistic function
#' @param a a real number 
#' @return a real number of value log(a/(1-a))
#' @export
#' @seealso \code{\link{expit}}
logit <- function(a){ -log(1/a-1) }

#' The Lambert W function
#' @description the Lambert W function is the inverse of f(x) = x*exp(x). The two main branches are b = 0 and b = -1. 
#' @param z a complex number 
#' @param b an integer, determining the branch of the W function
#' @param maxiter maximum number of Halley iteration
#' @param eps floating point tolerance
#' @param min.imag floating point tolerance for imaginary part
#' @return a real number of value log(a/(1-a))
#' @export
lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,
                    min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}

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
    if (q > 0) (sum(abs(v)^q))^(1/q)
    else if (q == 0) sum(v!=0)
    else if (q == Inf) max(abs(v))
    else NaN
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
matrix.standardise <- function(X){
  X <- sweep(X, 2, colMeans(X))
  X <- apply(X, 2, vector.normalise)
  return(X)
}

#' QL decomposition of a matrix
#' @param A an nxp matrix for QL decomposition
#' @return a list of two matrices Q and L so that A = QL
#' \itemize{
#'   \item Q - nxp matrix with orthonormal columns with the same span as A
#'   \item L - a lower triangular pxp matrix 
#' }
ql <- function(A){
  B <- A[,ncol(A):1]
  tmp <- qr(B)
  Q <- qr.Q(tmp)
  R <- qr.R(tmp)
  Q <- Q[,ncol(Q):1]
  L <- R[nrow(R):1,ncol(R):1]
  return(list(Q=Q, L=L))
}

#' extracting the off-diagonal part of A
#' @param A a rectangular matrix 
#' @param as.vector if FALSE, a matrix with main diagonal set to 0 is returned
#' otherwise, a vector of offdiagonal entries is returned
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

#' Compute matrix square-root of a symmetric matrix
#' @param A a square matrix
#' @param power power exponent
#' @param pseudoinverse whether to use pseudoinverse if power is negative
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

############### Auxiliary functions ###############

#' Print percentage
#' @param ind a vector of for loop interator
#' @param tot a vector of for loop lengths
#' @return on screen output of percentage
#' @export
printPercentage <- function (ind, tot){
    ind <- as.vector(ind); tot <- as.vector(tot)
    if ((length(tot) > 1) & (length(ind) == 1)) {ind <- match(ind, tot); tot <- length(tot)}
    len <- length(ind)
    contrib <- rep(1,len)
    if (len > 1) {
        for (i in (len-1):1) contrib[i] <- contrib[i+1] * tot[i+1]
    }
    grand_tot <- contrib[1] * tot[1]
    count <- (sum(contrib * (ind - 1)) + 1) 
    out <- ""
    if (sum(ind-1)>0) out <- paste0(rep("\b", nchar(round((count-1)/grand_tot * 100))+1), collapse = "")
    out <- paste0(out, round(count/grand_tot*100), "%")
    if (identical(ind, tot)) out <- paste0(out, '\n')
    cat(out)
    return(NULL)
}

#' Visualise a matrix X
#' @param X a matrix
#' @param aspect.ratio if automatic, it will be calculated automatically to fit screen, otherwise, the actual dimension of the matrix will be used.
#' @param axes whether to display axes
#' @param frame.plot whether to draw a frame
#' @return a color plot of matrix value magnitude 
#' @export
visualise <- function(X, aspect.ratio = c('automatic', 'actual'), axes = FALSE, frame.plot = FALSE){
    aspect.ratio = match.arg(aspect.ratio)
    n = dim(X)[1]; p = dim(X)[2]
    if (aspect.ratio == 'actual') {
        image(t(X[n:1,]),asp=n/p, axes = axes, frame.plot = frame.plot)
    }
    else {
        image(t(X[n:1,]), axes = axes, frame.plot = frame.plot)
    }
}

#' Show snippet of a large vector/matrix
#' @param A a vector, matrix or array
snippet <- function(A, nrow=5, ncol=nrow){
  if (is.vector(A)){
    cat('Vector of length ', length(A), ', with leading entries:\n', sep='')
    print(A[seq_len(min(length(A), nrow))])
  } else if (is.matrix(A)) {
    cat('Matrix with shape (', paste(as.character(dim(A)), collapse=', '), 
        '), with leading entries:\n')
    print(A[seq_len(min(nrow, nrow(A))), seq_len(min(ncol, ncol(A)))])
  } else if (is.array(A)) {
    dims <- dim(A); d <- length(dims); 
    shape <- paste(as.character(dim(A)), collapse=', ')
    if (d == 1){
      cat('1-d array of length ', dims, ', with leading entries:\n', sep='')
      print(A[seq_len(min(length(A), nrow))])
    } else if (d == 2){
      cat('2-d array with shape (', shape, '), with leading entries:\n')
      print(A[seq_len(min(nrow, nrow(A))), seq_len(min(ncol, ncol(A)))])
    } else {
      frames <- rep(0, d-2); starting_index <- 0
      for (i in seq_len(d-2)){
        frames[d-1-i] <- sample(dims[d+1-i], 1)
        starting_index <- starting_index + prod(head(dims, d-i)) * (frames[d-1-i] - 1)
      }
      cat(d, '-d array with shape (', shape, '), with leading entries in frame [:, :, ',
          paste(as.character(frames), collapse=', '), ']:\n', sep='')
      M <- matrix(A[starting_index + seq_len(dims[1]*dims[2])], dims[1], dims[2])
      print(M[seq_len(min(nrow, nrow(M))), seq_len(min(ncol, ncol(M)))])
    }
  } else {
    stop('A need to be a vector or a matrix or an array.')
  }
}

#' Find the location of first TRUE value in a boolean vector
#' @param v a logical vector
#' @return an integer denotating the location, return NA if not found.
#' @export
find.first <- function(v){
  match(TRUE, v, nomatch = NA)
}

#' Find the location of final TRUE value in a boolean vector
#' @param v a logical vector
#' @return an integer denotating the location, return NA if not found.
#' @export
find.last <- function(v){
  n <- length(v)
  n + 1L - match(TRUE, rev(v), nomatch = NA)
}

#' display signif of exponentiated number nicely
#' @details significant figure computed after subtracting 1. keep trailing zeros, not use scientific notation
#' @param x a real number
#' @param digits positive integer, number of significant figures
#' @return a string
#' @export
sf_exp <- function(x, digits){
  as.character(as.numeric(sf(x-1,digits))+1)
}

#' display signif nicely
#' @details keep trailing zeros, not use scientific notation
#' @param x a real number
#' @param digits number of significant figures to keep
#' @return a string
#' @export
sf <- function(x, digits){
  formatC(signif(x, digits=digits), digits=digits, format="fg", flag="#")
}

#' Simulation parameter data frame generation
#' @description  create a dataframe of all possible parameter combinations in lexicographic order (if tags are supplied, use tag for column names)
#' @param ... each argument should be of the form of tag = vector, meaning the variable named 'tag' takes values in 'vector'.
#' @details A sample usage is sim.params(tag1 = vec1, tag2 = vec2, tag3 = vec3).
#' @export
sim.params <- function(...){
  x <- list(...)
  n <- length(x)
  vnames <- names(x); no.vn <- !nzchar(vnames)
  vnames[no.vn] <- paste0('Var', seq_len(n))[no.vn]
  df <- expand.grid(rev(x))[,rev(seq_len(n))]
  colnames(df) <- vnames
  return(df)
}

#' Show parameter values
#' @description Print out parameters in a list in a nice format
#' @param x a list of parameters
#' @export
show.params <- function(x){
  n <- length(x)
  vnames <- names(x); no.vn <- !nzchar(vnames)
  vnames[no.vn] <- paste0('Var', seq_len(n))[no.vn]
  names(x) <- vnames
  str <- ""
  for (i in seq_len(n)){
    str <- paste(str, vnames[i], '=', paste(x[[i]]), ',')
  }
  substr(str, 2, nchar(str)-2)
}

#' Multiple assignment 
#' @description assign multiple items in a list on RHS to multiple items in a list on LHS
#' @details A sample usage is  bunch(a,b,c) %=% list('hello', 123, list('apple','orange')), or
#' bunch(a,b,c) %=% 1:3
#' @param l left side list, enclosed by the \code{bunch} function
#' @param r right side list
#' @export

'%=%' <- function(l, r) UseMethod('%=%')  # Generic form

# Binary Operator
'%=%.lbunch' = function(l, r) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (i in seq_along(l)) {
    do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
bunch = function(...) {
  List <- as.list(substitute(list(...)))[-1L]
  class(List) <- 'lbunch'
  return(List)
}


########### Random element generation ##########

#' generate n rademacher random variables
#' @param n length of random vector
#' @export
random.rademacher <- function(n=1){
    sample(c(-1,1),n,replace=T)
}

#' generate n random Bernoulli variables
#' @param n length of random vector
#' @param p probability (can be a vector)
#' @export
random.bernoulli <- function(n=1, p=0.5) {
    rbinom(n,1,p)
}

#' generate a random unit vectors in R^n
#' @param n length of random vector
#' @export
random.UnitVector <- function(n){
    v = rnorm(n)
    v/vector.norm(v)
}

#' generate a random nxn orthogonal matrix 
#' @param n dimension
#' @export
random.OrthogonalMatrix <- function(n){
    A = matrix(rnorm(n*n),n,n)
    qr.Q(qr(A))
}

#' generate a random nxn Wishart matrix 
#' @param n dimension
#' @export
random.WishartMatrix <- function(n){
    A = matrix(rnorm(n*n),n,n)
    t(A)%*%A
}

#' generate a random nxn Wigner matrix, where off diagonals are N(0,1) and diagonals are N(0,2)
#' @param n dimension
#' @export
random.WignerMatrix <- function(n){
    M = matrix(0,n,n)
    M[upper.tri(M, diag = TRUE)] = rnorm(n*(n+1)/2)
    M + t(M)
}

#' generate a random positive semidefinite matrix with operator norm <=1
#' @param n dimension
#' @export
random.psdMatrix <- function(n){
    Lambda = diag(sort(runif(n),decreasing = TRUE))
    V = random.OrthogonalMatrix(n)
    V%*%Lambda%*%t(V)
}

#' generate a random symmetric matrix
#' @param n dimension
#' @export
random.SymmetricMatrix <- function(n){
    Lambda = diag(sort(runif(n,min=-1,max=1),decreasing = TRUE))
    V = random.OrthogonalMatrix(n)
    V%*%Lambda%*%t(V)
}

########## statistical functions ##########

#' Cramer--von-Mise two sample test
#' @param x a vector
#' @param y a vector
#' @return a list containing test statistic and significance code
#' @export
CvM.test <- function(x,y){
    a <- sort(x); b <- sort(y)
    N <- length(a); M <- length(b)
    z <- c(a,b); temp <- rank(z);
    r <- temp[1:N]; s <- temp[(N+1):(N+M)];
    test.stat <- N*M/(N+M)^2*(sum((r/M - (1:N)*(1/M+1/N))^2) + sum((s/N - (1:M)*(1/M+1/N))^2))
    test.mean <- 1/6+1/6/(M+N)
    test.sd <- sqrt(1/45*(M+N+1)/(M+N)^2*(4*M*N*(M+N)-3*(M^2+N^2)-2*M*N)/(4*M*N))
    test.stat.normalised <- (test.stat - test.mean)/test.sd*(1/sqrt(45)) + 1/6
    signif <- ' '
    if (test.stat.normalised > 0.3473) signif <- '.'
    if (test.stat.normalised > 0.46136) signif <- '*'
    if (test.stat.normalised > 0.74346) signif <- '**'
    if (test.stat.normalised > 1.16786) signif <- '***'
    ret = list(t = test.stat.normalised, signif = signif)
    ret
}

########## string operations ##########

#' String length
#' @param str a string
#' @return its length
#' @export
strlen = function(str){
  return(nchar(str))
}

#' Substring search
#' @description find the position of the nth occurrence of needle in haystack, returns 0 if not found
#' @param haystack a string
#' @param needle substring to search for
#' @param startpos start position for search
#' @param n the nth occurrence
#' @return an integer
#' @export
strstr <- function(haystack, needle, startpos=1, n=1){
  aa <- unlist(strsplit(substring(haystack, startpos), needle))
  if(length(aa) < n + 1 ) return(0);
  return(sum(nchar(aa[1:n])) + startpos + (n-1) * nchar(needle) )
}

#' Extract left substring
#' @param str a string
#' @param len length of substring
#' @export
#' 
prefix <- function(str, len){
  substr(str, 1, len)
}

#' Extract right substring
#' @param str a string
#' @param len length of substring
#' @export
#' 
suffix <- function(str, len){
  substr(str, strlen(str) - len + 1, strlen(str))
}

########### NA handling ##########
#' Change all NA values in v to a
#' @param v a vector
#' @param a target value
#' @return updated vector
#' @export
setNA <- function(v, a){
  v[is.na(v)] <- a;
  return(v)
}