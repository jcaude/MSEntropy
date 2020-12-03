## Requirements
require(parallelMap)

#' Multi Scale Entropy (MSE)
#'
#' This function calculates the multiscale entropy (MSE) of a univariate signal x
#'
#' @param x a univariate signal (vector)
#' @param m the embedding dimension (default: 2)
#' @param r the tolerance factor for "similarity" (default: 0.15)
#' @param scales the scale factors as a vector (default 1:10)
#'
#' @return a vector of size 1 * length(scales) - the RCMDE of x
#'
#' @export
#' @importFrom stats sd
#' @importFrom parallelMap parallelLapply
#'
#' @references
#' \enumerate{
#'    \item  Costa M, Goldberger AL, Peng CK (2005) "Multiscale entropy analysis of biological signals". Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 71.
#' }
#'
#' @examples
#' data(EG_181117)
#' MSE(EG_181117)
MSE <- function(x,m=2,r=0.15,scales=1:10) {

  # init.
  mu <- mean(x)
  sigma <- sd(x)

  # calc.
  mse <- sapply(scales, function(s) {
    if(s == 1)
      xs <- x
    else
      xs <- CoarseGraining(x = x, scale = s)
    sampen <- SampEn(x = xs, m = m, r= r, sd = sigma)
    return(sampen)
  })
  names(mse) <- scales

  # eop
  return(mse)
}
