
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
randomSignal <- function(n) {

  # init.
  moves <- round(runif(n,min = -1,max = 1))
  x <- vector(mode = "integer",length = n)
  x[1] <- 0;

  # generate a random path
  for(i in 2:n) {
    x[i] <- x[i-1] + moves[i]
  }

  # eop
  return(x)
}
