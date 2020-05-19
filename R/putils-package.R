#' putils: A package for personal utility functions.
#'
#' The putils package provides personal utility functions.
#'
#' @section mathematical functsions:
#' \itemize{
#'   \item expit
#'   \item logit
#'   \item lambertW
#' }
#'
#' @section matrices and vectors:
#' \itemize{
#'   \item vector.norm
#'   \item vector.norm
#'   \item vector.normalise
#'   \item vector.clip
#'   \item vector.soft.thresh
#'   \item vector.hard.thresh
#'   \item vector.scramble
#'   \item vector.cleanup
#'   \item matrix.GramSchmidt
#'   \item matrix.trace
#'   \item matrix.rank
#'   \item matrix.standardise
#'   \item matrix.power
#'   \item sinThetaLoss
#'   \item powerMethod
#'   \item QR
#'   \item QL
#'   \item offdiag
#' }
#'
#' @section random element generation:
#' \itemize{
#'   \item random.rademacher
#'   \item random.bernoulli
#'   \item random.UnitVector
#'   \item random.OrthogonalMatrix
#'   \item random.WishartMatrix
#'   \item random.WignerMatrix
#'   \item random.psdMatrix
#'   \item random.SymmetricMatrix
#' }

#' @section auxiliary functions:
#' \itemize{
#'   \item printPercentage
#'   \item println
#'   \item write.latextable
#'   \item visualise
#'   \item snippet
#'   \item find.first
#'   \item find.last
#'   \item dp
#'   \item sf
#'   \item sf_exp
#'   \item sim.params
#'   \item show.params
#'   \item "%=%"
#' }

#' @section statistical functions:
#' \itemize{
#'   \item CvM.test
#' }
#'
#' @section string operations:
#' \itemize{
#'   \item strlen
#'   \item strstr
#'   \item prefix
#'   \item suffix
#' }
#'
#' @section NA handling:
#' \itemize{
#'   \item setNA
#' }
#'
#' @docType package
#' @name putils
#'
#' @importFrom graphics image
#' @importFrom stats na.omit rbinom rnorm runif
NULL
