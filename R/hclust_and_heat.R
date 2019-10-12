#' Hierarchical Clustering & visulization by heatmap
#'
#' Hierarchical cluster analysis and visualization by \pkg{ezlimmaplot} \code{ezheat}.
#'
#' @param annot Feature annotation.
#' @inheritParams hclust_in_grp
#' @inheritParams ezlimmaplot::ezheat
#' @inheritParams dynamicTreeCut::cutreeDynamic
#' @return A data.frame of the clusters assingment and feature annotation.
#' @export
#'
hclust_and_heat <- function(object, annot, sc="z", clip=NA, dist.method="euclidean", hc.method="ward.D2",
                            deepSplit=1, minClusterSize=10, verbose=TRUE,
                            pheno.df=NULL, labrows="", labcols="", color.v=NULL,
                            main="Log2 Expression", name="hclust_heat", reorder_cols=FALSE,
                            gaps_col=NULL, width=NA, height=NA, plot=TRUE){

  stopifnot(rownames(object) %in% rownames(annot))

  # sc
  stopifnot(sc %in% c("ctr", "z", "none"))
  if (sc=="ctr"){
    object.sc <- t(scale(x=t(object), center=TRUE, scale=FALSE))
    main <- paste("Centered", main)
  } else if (sc=="z"){
    object.sc <- t(scale(x=t(object), center=TRUE, scale=TRUE))
    main <- paste("Z-scored", main)
  } else{
    object.sc <- object
  }

  # clip
  if (!is.na(clip)){
    stopifnot(length(clip)==1, clip > 0)
    object.sc[object.sc < -clip] <- -clip
    object.sc[object.sc > clip] <- clip
  }

  # cluster
  dist_mat <- dist(object.sc, method=dist.method)
  hc <- hclust(dist_mat, method=hc.method)
  lbs <- hc$labels[hc$order]

  # cut tree
  clus <- unname(dynamicTreeCut::cutreeDynamic(hc, distM=as.matrix(dist_mat), deepSplit=deepSplit,
                               minClusterSize=minClusterSize, verbose=0))
  if (verbose){
    print(knitr::kable(table(clus), caption="Clusters"))
  }

  # sort
  names(clus) <- rownames(object.sc)
  clus <- clus[lbs]
  clus <- sort(clus)
  lbs <- names(clus)

  clus_df <- as.data.frame(clus)
  colnames(clus_df) <- "Cluster"
  clus_df$Cluster <- factor(paste0("clus_", clus_df$Cluster), levels=paste0("clus_", sort(unique(clus))))

  # heatmap
  gaps_row <- which(diff(clus, lag=1) != 0)
  annotation_colors <-
    list(Cluster=setNames(colorRampPalette(RColorBrewer::brewer.pal(n=min(max(clus), 9), name="Set1"))(max(clus)),
                          nm=levels(clus_df$Cluster)))
  ph <- ezheat(object.sc[lbs, ], pheno.df=pheno.df, sc="none", reorder_rows=FALSE, reorder_cols=reorder_cols,
               labrows=labrows, labcols=labcols, color.v=color.v, gaps_col=gaps_col,
               gaps_row=gaps_row, annotation_row=clus_df, annotation_colors=annotation_colors,
               main=main, name=name, height=height, width=width, plot=plot)[["mat"]]

  clus_df <- data.frame(clus_df[rownames(ph), ,drop=FALSE], annot[rownames(ph), , drop=FALSE], check.names=FALSE)

  return(clus_df)
}
