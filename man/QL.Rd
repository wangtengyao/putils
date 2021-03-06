% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linalg.R
\name{QL}
\alias{QL}
\title{QL decomposition of a matrix}
\usage{
QL(A)
}
\arguments{
\item{A}{an nxp matrix for QL decomposition}
}
\value{
a list of two matrices Q and L so that A = QL
\itemize{
  \item Q - nxp matrix with orthonormal columns with the same span as A
  \item L - a lower triangular pxp matrix with nonnegative diagonal entries
}
}
\description{
QL decomposition of a matrix
}
