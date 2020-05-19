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
