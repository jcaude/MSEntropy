library(tuneR)

#' Generate a random signal
#'
#' This function generate a random signal using a very simple random walk
#'
#' @param n the length of the signal to generate
#'
#' @return a vector
#' @export
#'
#' @examples
#' x <- randomSignal(1000)
#' \donttest{plot(x,type='l')}
randomSignal <- function(n, model="white", alpha=1) {

  # forge the signal
  if(model == "power")
    x <- noise(kind = model, duration = n, alpha=alpha, stereo = FALSE)
  else
    x <- noise(kind = model, duration = n, stereo = FALSE)
  x <- x@left

  # eop
  return(x)
}
