########## Statistical functions ##########

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

#' EM algorithm
#' @description EM algorithm with identity covariance matrix for all groups
#' @param X data matrix of dimension n x p
#' @param num_clusters number of clusters
#' @param tol tolerance for convergence
#' @return a list: mu_hat - the estimated cluster centroids and soft_label -
#' the estimated cluster membership
EM <- function(X, num_clusters, tol=1e-10){
  # randomly initialise cluster centroids
  mu_hat <- matrix(rnorm(ncol(X)*num_clusters), num_clusters)
  sumX2 <- rowSums(X^2) # precompute squared length of X to speed up dist calc

  # EM iterations
  repeat{
    # E-step: compute soft labels
    # construct squared distance matrix between X and mu_hat
    dist2 <- outer(sumX2, rowSums(mu_hat^2), '+') - 2 * X %*% t(mu_hat)
    lik <- exp(-dist2/2)
    soft_labels <- lik / rowSums(lik)

    # M-step:update cluster centroids
    old_mu_hat <- mu_hat
    mu_hat <- t(soft_labels) %*% X / colSums(soft_labels)

    if (sum((mu_hat - old_mu_hat)^2) < tol) break
  }

  return(list(mu_hat=mu_hat,
              labels = apply(soft_labels, 1, which.max),
              soft_labels=soft_labels))
}
