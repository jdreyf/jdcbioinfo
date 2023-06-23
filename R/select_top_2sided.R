#' Select equal numbers of top features from the topTable
#'
#' Select equal numbers of top features from the topTable based on p-values.
#'
#' @param ntop Number of top features
#' @param ntop.offset Offset for ntop
#' @param each Logical TRUE for select equal number of top features from each comparison
#' @inheritParams multi_pval2z
#' @return Vector of the top feature's IDs (the rownames of topTable).
#' @export

select_top <-  function(tab, prefix.v=NULL, p.suffix="p", direction.suffix="logFC", ntop=50, ntop.offset=5, each=TRUE){

  if(nrow(tab)<=ntop) return(rownames(tab))

  mat.z <- multi_pval2z (tab, prefix.v=prefix.v, p.suffix=p.suffix, direction.suffix=direction.suffix)

  if(!each){
    top <- rownames(tab)[order(ezlimma::combine_pvalues(tab[, p.cols]))][1:ntop]
  } else {
    ord <- list()
    ntop.each <- ceiling(ntop/2)
    for (i in 1:ncol(mat.z)) {
      ord[[paste0("up", i)]] <- rownames(mat.z)[order(-mat.z[,i])][1:ntop.each]
      ord[[paste0("dn", i)]] <- rownames(mat.z)[order(mat.z[,i])][1:ntop.each]
    }

    ntop.tmp <- Inf
    while(ntop.tmp-ntop>ntop.offset){
      top <- lapply(ord, function(x) x[1:ntop.each])
      top <- Reduce(union, top)
      ntop.tmp <- length(top)
      ntop.each <- ntop.each - 1
    }
  }
  return(top)
}
