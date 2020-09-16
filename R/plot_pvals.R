#' Plot  p-values
#'
#' Plot p-values
#'
#' @param pvals A matrix of p-values.
#' @param name Name of the plot file.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @return NULL

plot_pvals <- function(pvals, name=NA, width=8, height=7) {

  stopifnot(pvals>=0, pvals<=1)

  signifSymbols <- c("***", "**", "*", "")
  signifCutpoints <- c(0, 0.001, 0.01, 0.05, 1)

  signif <- pvals
  for (i in 1:ncol(pvals)) {
    signif[, i] <- c(stats::symnum(pvals[, i], corr=FALSE, na=FALSE, cutpoints=signifCutpoints, symbols=signifSymbols))
  }

  plotLabels <- pvals
  for (i in 1:nrow(pvals)) {
    for (j in 1:ncol(pvals)) {
      plotLabels[i, j] <- paste(round(pvals[i, j], 3), signif[i, j], sep="")
    }
  }

  posLab=1
  axisTicks=c(1, 0)
  cols <- grDevices::colorRampPalette(rev(c("white", "cornsilk1", "gold", "forestgreen", "darkgreen")))

  maxp <- max(pvals)
  minp <- min(pvals)
  iUpperRange <- min(1, maxp + 0.01)
  iLowerRange <- max(0, minp - 0.01)


  labels <- function(x, y, z, ...) {
    lattice::panel.levelplot(x, y, z, ...)
    lattice::ltext(x, y, labels=plotLabels, cex=1, col="black", font=1)
  }

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"), width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  lattice::levelplot(pvals, xlab=list(label="",cex=1, rot=0, col="black", font=2),
                     ylab=list(label="", cex=1, rot=0, col="black", font=2), panel=labels,
                     pretty=TRUE, par.settings=list(panel.background=list(col="white")),
                     scales=list(x=list(cex=1, rot=0, col="black", font=2), y=list(cex=1, rot=0, col="black", font=2),
                                 tck=axisTicks, alternating=posLab), aspect="fill",
                     col.regions=cols, cuts=100, at=seq(iLowerRange, iUpperRange, 0.01),
                     main=list(label="p-Values", cex=2, rot=0, col="black", font=2),
                     colorkey=list(space="right", labels=list(cex=1)))
}
