#' Mode estimation of a continuous variable
#'
#' Kernel density estimation of modes of a continuous variable.
#'
#' @param bin  the number of bins for the histogram.
#' @param by increment of the bins i.e. the width of bins.
#' @param plot logical for ploting the histogram.
#' @inheritParams graphics::hist
#' @return Data frame with modes values and densities.
#' @export

find_mode <- function(x, xlab=NULL, main="", by=NULL, bin=20, plot=TRUE) {

  lim.inf <- min(x)-1
  lim.sup <- max(x)+1

  if (is.null(by)) by <- (lim.sup - lim.inf)/bin
  if (is.null(xlab)) xlab <- deparse(substitute(x))

  s <- stats::density(x, from=lim.inf, to=lim.sup, bw=by)
  n <- length(s$y)
  v1 <- s$y[1:(n-2)];
  v2 <- s$y[2:(n-1)];
  v3 <- s$y[3:n]
  ix <- 1+which((v1 < v2) & (v2 > v3))

  if (plot) {
    graphics::hist(x, freq=FALSE, breaks=seq(lim.inf,lim.sup, by=by), main=main, xlab=xlab, col="grey80")
    graphics::lines(s$x, s$y, col="red")
    graphics::points(s$x[ix], s$y[ix],col="blue")
  }

  return(data.frame(mds.vals=s$x[ix], mds.dens=s$y[ix]))
}
