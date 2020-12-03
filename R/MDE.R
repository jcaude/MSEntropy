## Requirements
require(parallelMap)

#' @useDynLib MSEntropy, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppArmadillo armadillo_version
NULL

#' Refined Composite Multiscale Dispersion Entropy (RCMDE)
#'
#' This function calculates the refined composite multiscale dispersion entropy (RCMDE) of a univariate signal x
#'
#' @param x a univariate signal (vector)
#' @param m the embedding dimension (default: 1)
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
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
#'    \item  Azami H, Rostaghi M, Abasolo D, Escudero J (2017) "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical Signals". IEEE transactions on bio-medical engineering 64:2872–2879.
#'    \item Rostaghi M, Letters HAISP, (2016) "Dispersion entropy: A measure for time-series analysis". IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
#' }
#'
#' @examples
#' data(EG_181117)
#' RCMDE(EG_181117)
RCMDE <- function(x,m=1,nc=6,tau=1,scales=1:10) {

  # init.
  nx <- length(x)
  mu <- mean(x)
  sigma <- sd(x)

  # other scales
  rcmde <- sapply(scales, function(s) {
    if(s == 1) {
      dispen <- DispEn(x = x, ma = "NCDF", m = m, nc = nc, tau = tau)
      mde <- dispen$disp.en
    }
    else {
      pdf <- parallelLapply(1:s, function(ss) {
        xs <- x[ss:nx]
        xs <- CoarseGraining(x = xs, scale = s)
        dispen <- DispEn(x = xs, ma = "NCDF", m = m, nc = nc,
                                   tau = tau, mu = mu, sigma = sigma)
        return(dispen$pdf)
      })
      pdf <- unlist(pdf)
      pdf <- matrix(pdf,nrow = s,ncol = length(pdf)/s, byrow = TRUE)
      pdf <- apply(pdf,2,mean)
      pdf <- pdf[pdf != 0]
      mde <- -sum(pdf * log(pdf))
    }
    return(mde)
  })
  names(rcmde) <- scales

  # eop
  return(rcmde)
}


#' MultiScale Dispersion Entropy (MDE)
#'
#' This function calculates the multiscale dispersion entropy (MDE) of a univariate signal x
#'
#' @param x a univariate signal (vector)
#' @param m the embedding dimension (default: 1)
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
#' @param scales the scale factors as a vector (default 1:10)
#'
#' @return a vector of size 1 * length(scales) - the MDE of x
#'
#' @export
#' @importFrom stats sd
#'
#' @references
#' \enumerate{
#'    \item  Azami H, Rostaghi M, Abasolo D, Escudero J (2017) "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical Signals". IEEE transactions on bio-medical engineering 64:2872–2879.
#'    \item Rostaghi M, Letters HAISP, (2016) "Dispersion entropy: A measure for time-series analysis". IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
#' }
#'
#' @examples
#' data(EG_181171)
#' MDE(EG_181117)
MDE <- function(x,m=1,nc=6,tau=1,scales=1:10) {

  # init.
  mu <- mean(x)
  sigma <- sd(x)

  # calc.
  mde <- sapply(scales, function(s) {
    if(s == 1)
      xs <- x
    else
      xs <- CoarseGraining(x = x, scale = s)
    dispen <- DispEn(x = xs, ma = "NCDF", m = m, nc = nc, tau = tau,
                               mu = mu, sigma = sigma)
    return(dispen$disp.en)
  })
  names(mde) <- scales

  # eop
  return(mde)
}
