#' Correlation test in each group
#'
#' Two-sided correlation test in each group and combine the directional p-values.
#'
#' @inheritParams ezlimma::ezcor
#' @inheritParams ezlimma::limma_contrasts
#' @return A data.frame of test statistics for each group and the combined statistics.
#' @export
#'

cor_by_grp <- function(object, phenotype, grp, method="pearson", reorder.rows=TRUE, prefix=NULL, check.names=TRUE){

  stopifnot(length(phenotype)==ncol(object), length(grp)==ncol(object))
  if (check.names){
    stopifnot(names(phenotype)==colnames(object), names(grp)==colnames(object))
  }

  # correlation within each group
  ct <- list()
  for (g in unique(sort(grp))){
    object_tmp <- object[, grp==g, drop=FALSE]
    phenotype_tmp <- phenotype[grp==g]

    if (!is.null(prefix)){
      prefix_tmp <- paste(prefix, g, sep='.')
    } else{
      prefix_tmp <- g
    }

    ct[[prefix_tmp]] <- ezlimma::ezcor(object_tmp, phenotype_tmp, method=method, reorder.rows=FALSE, prefix=prefix_tmp,
                                       adjust.method="BH", alternative="two.sided")
  }
  ct <- Reduce(cbind, ct)

  # combine directional p-values
  direction_suffix <- gsub(".*\\.", "", colnames(ct)[1])
  comb <- combine_directional_pvalues(ct, direction.suffix=direction_suffix, prefix=prefix)
  ct <- cbind(ct, comb)
  if (reorder.rows) ct <- ct[order(ct[, grep("Combined.p", colnames(ct))]), ]

  return(ct)
}
