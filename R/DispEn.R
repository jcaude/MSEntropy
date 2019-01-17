#' Dispersion Entropy
#'
#' This function calculates dispersion entropy of a univariate signal
#' x, using different mapping approaches (MA)
#'
#' @param x a univariate signal (vector)
#' @param ma the mapping approach, possible values are:
#' \itemize{
#'    \item \code{LM}: linear mapping
#'    \item \code{NCDF}: normal cumulative distribution function (default)
#'    \item \code{LOGSIG}: logarithm sigmoid
#'    \item \code{TANSIG}: tangent sigmoid
#'    \item \code{SORT}: sorting method
#' }
#' @param m the embedding dimension
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
#' @param mu an optional defined mean (the observed mean by default)
#' @param sigma an option defined std-deviation (the observed std-deviation by default)
#'
#' @return a named list:
#' \itemize{
#'    \item \code{disp.en}: the Dispersion Entropy value of the signal x
#'    \item \code{pdf}: a vector of length nc^m, showing the normalized number of disersion patterns of x
#' }
#' @export
#'
#' @references
#' \enumerate{
#'    \item Azami H, Entropy JE, (2018) "Amplitude-and Fluctuation-Based Dispersion Entropy". Entropy 20:210.
#'    \item Rostaghi M, Letters HAISP, (2016) "Dispersion entropy: A measure for time-series analysis". IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
#' }
#'
#' @examples
#' data(EG_181117)
#' DispEn(EG_181117,m=1)
DispEn <- function(x,ma="NCDF",m=3,nc=6,tau=1,mu,sigma) {

  # init.
  N <- length(x)
  sigma_x <- ifelse(missing(sigma),NA,sigma)
  mu_x <- ifelse(missing(mu),NA,mu)

  #
  switch (ma,

          LM = {
            z <- dispen_map(x,1,nc,mu_x,sigma_x)
          },

          NCDF = {
            z <- dispen_map(x,2,nc,mu_x,sigma_x)
          },

          LOGSIG = {
            z <- dispen_map(x,3,nc,mu_x,sigma_x)
          },

          TANSIG = {
            z <- dispen_map(x,4,nc,mu_x,sigma_x)
          },

          SORT = {
            z <- dispen_map(x,5,nc,mu_x,sigma_x)
            N <- length(z)
          },

          {
            stop("Invalid Mapping Approach (MA)")
          }
  )

  # compute the normalize PDF using step2
  npdf <- dispen_npdf(z,nc,m,tau)

  # calc. final results
  p <- npdf[npdf > 0]
  Out_DispEn <- -sum(p * log(p))

  # eop
  return(list(disp.en=Out_DispEn, pdf=npdf))
}


#' Fluctuation Based Dispersion Entropy
#'
#' This function calculates fluctuation-based dispersion entropy (FDispEn) of a
#' univariate signal x, using different mapping approaches (MA)
#'
#' @param x a univariate signal (vector)
#' @param ma the mapping approach, possible values are:
#' \itemize{
#'    \item \code{LM}: linear mapping
#'    \item \code{NCDF}: normal cumulative distribution function (default)
#'    \item \code{LOGSIG}: logarithm sigmoid
#'    \item \code{TANSIG}: tangent sigmoid
#'    \item \code{SORT}: sorting method
#' }
#' @param m the embedding dimension
#' @param nc the number of classes (it is usually equal to a number between 3 and 9, we use 6 by default)
#' @param tau the time lag (it is usually equal to 1)
#' @param mu an optional defined mean (the observed mean by default)
#' @param sigma an option defined std-deviation (the observed std-deviation by default)
#'
#' @return a named list:
#' \itemize{
#'    \item \code{fdisp.en}: the Fluctuation-Based Dispersion Entropy value of the signal x
#'    \item \code{pdf}: a vector of length nc^m, showing the normalized number of disersion patterns of x
#' }
#' @export
#'
#' @references
#' \enumerate{
#'    \item Azami H, Entropy JE, (2018) "Amplitude-and Fluctuation-Based Dispersion Entropy". Entropy 20:210.
#'    \item Rostaghi M, Letters HAISP, (2016) "Dispersion entropy: A measure for time-series analysis". IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
#' }
#'
#' @examples
#' data(EG_181117)
#' FDispEn(EG_181117,m=1)
FDispEn <- function(x,ma="NCDF",m=3,nc=6,tau=1,mu,sigma) {

  # init.
  N <- length(x)
  sigma_x <- ifelse(missing(sigma),NA,sigma)
  mu_x <- ifelse(missing(mu),NA,mu)

  #
  switch (ma,

          LM = {
            z <- dispen_map(x,1,nc,mu_x,sigma_x)
          },

          NCDF = {
            z <- dispen_map(x,2,nc,mu_x,sigma_x)
          },

          LOGSIG = {
            z <- dispen_map(x,3,nc,mu_x,sigma_x)
          },

          TANSIG = {
            z <- dispen_map(x,4,nc,mu_x,sigma_x)
          },

          SORT = {
            z <- dispen_map(x,5,nc,mu_x,sigma_x)
            N <- length(z)
          },

          {
            stop("Invalid Mapping Approach (MA)")
          }
  )

  # compute the normalize PDF using step2
  npdf <- fdispen_npdf(z,nc,m,tau)

  # calc. final results
  p <- npdf[npdf > 0]
  Out_FDispEn <- -sum(p * log(p))

  # eop
  return(list(fdisp.en=Out_FDispEn, pdf=npdf))
}
