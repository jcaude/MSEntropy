## This file is an R port of the matlab source code from Azami H. & al. 2017,2018


#' This function calculates dispersion entropy of a univariate signal
#' x, using different mapping approaches (MA)
#'
#' @param x
#' @param ma
#' @param m
#' @param nc
#' @param tau
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
DispersionEntropy <- function(x,ma="NCDF",m=3,nc=6,tau=1,mu,sigma) {

  # local
  .mapminmax <- function(x,ymin,ymax) {
    xmin <- min(x)
    xmax <- max(x)
    if(xmin == xmax)
      y <- min(ymax,max(x,ymin))
    else
      y <- (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin
    return(y)
  }

  # init.
  N <- length(x)
  sigma_x <- ifelse(missing(sigma),sd(x),sigma)
  mu_x <- ifelse(missing(mu),mean(x),mu)

  #
  switch (ma,

          LM = {
            y <- .mapminmax(x,0,1)
          },

          NCDF = {
            y <- pnorm(x,mean = mu_x, sd = sigma_x)
          },

          LOGSIG = {
            y <- (x-mu_x)/sigma_x
            y <- 1 / (1 + exp(-y))
            y <- .mapminmax(y,0,1)
          },

          TANSIG = {
            y <- (x-mu_x)/sigma_x
            y <- 2/(1+exp(-2*y))-1
            y <- .mapminmax(y,0,1)
          },

          SORT = {
            x <- x[1:(nc*floor(N/nc))]
            N <- length(x)
            osx <- order(x)
            cx <- rep(1:nc,each=N/nc)
            z <- sapply(1:N,function(i) {cx[osx == i]})
          },

          {
           stop("Invalid Mapping Approach (MA)")
          }
  )

  #
  if(ma != "SORT") {
    y[y == 1] <- 1-(1e-10)
    y[y == 0] <- 1e-10
    z <- round(y*nc+0.5);
  }

  #
  all_patterns <- matrix(0,nrow = nc^m, ncol = m)
  for(idx in 1:m)
    all_patterns[,idx] <- rep.int(rep(1:nc,each=nc^(idx-1)),times = nc^(m-idx))

  #
  key_coef <- 10^(m:1-1)
  key <- apply(all_patterns,1,function(v){sum(v*key_coef)})

  #
  embd2 <- vector(mode = "numeric", length = N-(m-1)*tau)
  for(idx in 1:m) {
    embd2 <- z[(1+(idx-1)*tau):(N-(m-idx)*tau)]*10^(m-idx) + embd2
  }

  #
  pdf <- sapply(1:nc^m, function(idx) {sum(embd2 == key[idx])})

  #
  npdf <- pdf/(N-(m-1)*tau)
  p <- npdf[npdf > 0]
  Out_DisEn <- -sum(p * log(p))

  # eop
  return(list(Out_DisEn=Out_DisEn, pdf=npdf))
}
