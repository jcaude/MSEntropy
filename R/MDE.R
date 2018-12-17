## #'    \item Azami H, Rostaghi M, Abasolo D, Escudero J (2017) "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical Signals". IEEE transactions on bio-medical engineering 64:2872–2879.


#' Refined Composite Multiscale Dispersion Entropy
#'
#' This function calculates the refined composite multiscale dispersion entropy (RCMDE) of a univariate signal x
#'
#' @param x a univariate signal (vector)
#' @param m the embedding dimension
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
#' @param nscale the number of scale factors
#'
#' @return a vector of size 1 * Scale - the RCMDE of x
#' @export
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
RCMDE <- function(x,m,nc=6,tau=1,nscale=5) {

  # init.
  mde <- vector(mode = "numeric",length = nscale)
  mde[1:nscale] <- NaN
  nx <- length(x)
  mu <- mean(x)
  sigma <- sd(x)

  # scale 1
  disen <- DispersionEntropy(x = x, ma = "NCDF", m = m, nc = nc, tau = tau)
  mde[1] <- disen$disp.en

  # other scales
  if(nscale > 2) {
    for(s in 2:nscale) {
      pdf <- lapply(1:s, function(ss) {
        xs <- x[ss:nx]
        xs <- ScaleX(x = xs, scale = s)
        disen <- DispersionEntropy(x = xs, ma = "NCDF", m = m, nc = nc,
                                   tau = tau, mu = mu, sigma = sigma)
        return(disen$pdf)
      })
      pdf <- unlist(pdf)
      pdf <- matrix(pdf,nrow = s,ncol = length(pdf)/s, byrow = TRUE)
      pdf <- apply(pdf,2,mean)
      pdf <- pdf[pdf != 0]
      mde[s] <- -sum(pdf * log(pdf))
    }
  }

  # eop
  return(mde)
}


#' MultiScale Dispersion Entropy
#'
#' This function calculates the multiscale dispersion entropy (MDE) of a univariate signal x
#'
#' @param x a univariate signal (vector)
#' @param m the embedding dimension
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
#' @param nscale the number of scale factors
#'
#' @return a vector of size 1 * Scale - the MDE of x
#' @export
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
MDE <- function(x,m,nc=6,tau=1,nscale=5) {

  # init.
  mde <- vector(mode = "numeric",length = nscale)
  mu <- mean(x)
  sigma <- sd(x)

  # scale 1
  disen <- DispersionEntropy(x = x, ma = "NCDF", m = m, nc = nc, tau = tau)
  mde[1] <- disen$disp.en

  # scale 2...
  if(nscale > 2) {
    for(s in 2:nscale) {
      xs <- ScaleX(x = x, scale = s)
      disen <- DispersionEntropy(x = xs, ma = "NCDF", m = m, nc = nc, tau = tau,
                                 mu = mu, sigma = sigma)
      mde[s] <- disen$disp.en
    }
  }

  # eop
  return(mde)
}


ScaleX <- function(x,scale) {
  xc <- seq_along(x)
  Rtrim <- max(which((xc %% scale) == 0))
  xc <- xc[1:Rtrim]
  x <- x[xc]
  xs <- split(x,ceiling(seq_along(x)/scale))
  xs <- sapply(xs,mean)
  names(xs) <- NULL
  return(xs)
}
