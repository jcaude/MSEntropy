RCMDE <- function(x,m,nc=6,tau=1,nscale=5) {

  # init.
  mde <- vector(mode = "numeric",length = nscale)
  mde[1:nscale] <- NaN
  nx <- length(x)
  mu <- mean(x)
  sigma <- sd(x)

  # scale 1
  disen <- DispersionEntropy(x = x, ma = "NCDF", m = m, nc = nc, tau = tau)
  mde[1] <- disen$Out_DisEn

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


#' Title
#'
#' @param x
#' @param m
#' @param nc
#' @param tau
#' @param nscale
#'
#' @return
#' @export
#'
#' @examples
MDE <- function(x,m,nc=6,tau=1,nscale=5) {

  # init.
  mde <- vector(mode = "numeric",length = nscale)
  mu <- mean(x)
  sigma <- sd(x)

  # scale 1
  disen <- DispersionEntropy(x = x, ma = "NCDF", m = m, nc = nc, tau = tau)
  mde[1] <- disen$Out_DisEn

  # scale 2...
  if(nscale > 2) {
    for(s in 2:nscale) {
      xs <- ScaleX(x = x, scale = s)
      disen <- DispersionEntropy(x = xs, ma = "NCDF", m = m, nc = nc, tau = tau,
                                 mu = mu, sigma = sigma)
      mde[s] <- disen$Out_DisEn
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
