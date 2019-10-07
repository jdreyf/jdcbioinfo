#' Hierarchical Clustering in each group
#'
#' Hierarchical cluster analysis in each group.
#'
#' @param dist.method The distance measure to be used.
#' @param hc.method The agglomeration method to be used.
#' @inheritParams ezlimmaplot::ezheat
#' @inheritParams ezlimma::limma_contrasts
#' @return Vector of the labels.
#' @export

hclust_in_grp <- function(object, grp, sc="ctr", clip=NA, dist.method="euclidean", hc.method="ward.D2"){

  stopifnot(nrow(object)==length(grp))

  # sc
  stopifnot(sc %in% c("ctr", "z", "none"))
  if (sc=="ctr") {
    object.sc <- t(scale(x=t(object), scale=FALSE))
  } else if (sc=="z") {
    object.sc <- t(scale(x=t(object), center=TRUE, scale=TRUE))
  } else {
    object.sc <- object
  }

  # clip
  if (!is.na(clip)) {
    stopifnot(length(clip)==1, clip>0)
    object.sc[object.sc < -clip] <- -clip
    object.sc[object.sc > clip] <- clip
  }

  # cluster
  hc.labels <- list()
  groups <- unique(grp)
  for (group in groups) {
    if (sum(grp==group) > 1) {
      object.sc.tmp <- object.sc[grp==group, ]
      object.dist.tmp <- dist(object.sc.tmp, method=dist.method)
      hc.tmp <- hclust(object.dist.tmp, method=hc.method)
      hc.labels[[group]] <- hc.tmp$labels[hc.tmp$order]
    } else {
      hc.labels[[group]] <- rownames(object.sc)[grp==group]
    }
  }
  hc.labels <- unlist(hc.labels)
  return(hc.labels)
}
