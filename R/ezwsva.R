#' Weighted Surrogate Variable Analysis
#'
#' Weighted surrogate variable analysis, with selection of number of surrogate variables.
#'
#' @param num An integer scalar specifying the number of surrogate variables at the elbow point.
#' @inheritParams limma::wsva
#' @return Data frame with statistics from rank products test.
#' @export

ezwsva <- function (y, design, n.sv = 1L, weight.by.sd = FALSE, plot = FALSE, num = NULL,
                    ...)
{
  y <- as.matrix(y)
  ngenes <- nrow(y)
  narrays <- ncol(y)
  p <- ncol(design)
  d <- narrays - p
  n.sv <- max(n.sv, 1L)
  n.sv <- min(n.sv, d)
  if (n.sv <= 0L)
    stop("No residual df")
  if (weight.by.sd) {
    if (plot)
      message("Plot not available with weight.by.sd=TRUE")
    for (i in 1L:n.sv) {
      Effects <- limma:::.lmEffects(y, design, ...)[, -1L]
      s <- sqrt(rowMeans(Effects^2))
      Effects <- s * Effects
      u <- drop(svd(Effects, nu = 1L, nv = 0L)$u)
      u <- u * s
      sv <- colSums(u * y)
      design <- cbind(design, sv)
    }
    SV <- t(design[, -(1:p), drop = FALSE])
  }
  else {
    Effects <- limma:::.lmEffects(y, design, ...)[, -1L]
    SVD <- svd(Effects, nu = n.sv, nv = 0L)
    SV <- crossprod(SVD$u, y)
    if (plot) {
      lambda <- SVD$d^2
      lambda <- lambda/sum(lambda)
      if (is.null(num)) {
        num <-  PCAtools::findElbowPoint(lambda)
      }
      plot(lambda, xlab = "Surrogate variable number",
           ylab = "Proportion variance explained")
      graphics::abline(v = num, col = "red")
    }
  }
  A <- rowMeans(SV^2)
  SV <- t(SV/sqrt(A))
  colnames(SV) <- paste0("SV", 1L:n.sv)
  return(SV=SV, num=num)
}
