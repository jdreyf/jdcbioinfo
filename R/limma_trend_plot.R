#' Mean-variance trend plot
#'
#' Estimate the mean-variance relationship and plot its trend using \code{limma} and \code{ggplot2}.
#' @param name file name for plot
#' @param alpha Transparency, passed to \code{\link[ggplot2]{geom_point}}.
#' @param ... Passed to \code{\link[limma]{lmFit}}.
#' @inheritParams limma::lmFit
#' @inheritParams limma::voom
#' @inheritParams grDevices::pdf
#' @export

limma_trend_plot <- function(object, design=NULL, name=NA, span=0.5, alpha=0.5, width=7, height=7, ...){

  if (is.null(design)) {
    design <- matrix(1, ncol(object), 1)
    rownames(design) <- colnames(object)
    colnames(design) <- "GrandMean"
  }else{
    stopifnot(rownames(design)==colnames(object))
  }

  fit <- limma::lmFit(object, design=design, ...)
  if (is.null(fit$Amean)) {
    fit$Amean <- rowMeans(object, na.rm=TRUE)
  }

  sx <- fit$Amean
  sy <- sqrt(fit$sigma)
  dat <- data.frame(sx=sx, sy=sy)

  l <- stats::lowess(sx, sy, f=span)
  l <- as.data.frame(l)

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"), width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  gp <- ggplot2::ggplot(dat=dat, mapping=ggplot2::aes(x=sx, y=sy)) + ggplot2::geom_point(alpha=alpha)
  gp <- gp + ggplot2::geom_line(mapping=ggplot2::aes(x=x, y=y), data=l, color="blue", linewidth=1)
  gp <- gp + ggplot2::labs(title="Mean-variance Trend", x="Mean Log2 Values", y="Sqrt( Srandard Deviation )")

  graphics::plot(gp)
  return(invisible(gp))
}
