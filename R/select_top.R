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

select_top <-  function(tab, prefix.v=NULL, p.suffix="p", ntop=50, ntop.offset=5, each=TRUE){

  if(nrow(tab)<=ntop) return(rownames(tab))

  if(!is.null(prefix.v)){
    p.cols <- paste(prefix.v, p.suffix, sep = ".")
    stopifnot(p.cols %in% colnames(tab))

  } else {
    p.suffix <- paste0("\\.", p.suffix, "$")
    p.cols <- grep(p.suffix, colnames(tab), value = TRUE)
  }

  if(!each){
    top <- rownames(tab)[order(ezlimma::combine_pvalues(tab[, p.cols]))][1:ntop]
  } else {
    ord <- lapply(p.cols, function(p.col) rownames(tab)[order(tab[, p.col])][1:ntop])

    ntop.each <- ntop
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
