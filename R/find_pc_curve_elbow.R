#' Identify and plot the elbow point of the PC curve
#'
#' Identify and plot the elbow point of the PC curve
#'
#' @param method 1 or 2. 1 = use \code{\link[pathviewr]{find_curve_elbow}}, 2 = use \code{\link[PCAtools]{findElbowPoint}}, default is 1.
#' @param maxPC maximum PCs to be considered. Default is 50.
#' @param width figure width.
#' @param height figure height.
#' @inheritParams ezheat
#' @return the PC number at elbow point.
#' @export

find_pc_curve_elbow <- function(object, method=1, maxPC=50, plot=TRUE, name=NA, width=6, height=4) {

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale. = FALSE)
  numPC <- min(ncol(object), maxPC)
  pve <- summary(pca)$importance["Proportion of Variance", 1:numPC]*100
  pcaDf <- data.frame(PC = 1:length(pve), PVE = pve)

  if (method==1) {
    elbow <- pathviewr::find_curve_elbow(pcaDf)
  } else if (method==2) {
    elbow <- PCAtools::findElbowPoint(pve)
  }

  ggp <- ggplot2::ggplot(pcaDf, mapping=ggplot2::aes(x=PC, y=PVE)) +ggplot2::geom_point(shape=1, size=2)
  ggp <- ggp + ggplot2::geom_vline(xintercept=elbow, linetype="dashed", linewidth=0.75, color="blue")
  ggp <- ggp + ggplot2::labs(title="Importance of PCs", x="Principal Components", y="Percent of Variance Explained")

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"), width=width, height=height)
    on.exit(grDevices::dev.off())
  }

  if (plot) graphics::plot(ggp)

  return(elbow)
}
