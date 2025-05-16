#' Convert multiple columns of p-values to list of named and sorted Z-scores
#'
#' Convert multiple columns (from a table-like object) of p-values from two-sided tests to list of named and sorted Z-scores.
#' Z-scores that are \code{NA} are removed.
#'
#' @param sort Logical. Should the Z-scores be sorted?
#' @param decreasing Logical. Should the sort be increasing or decreasing?
#' @param nm  Names of the Zscores.
#' @inheritParams multi_pval2z
#' @return List of Z-scores.
#' @export

multi_pval2zlist <- function(tab, prefix.v=NULL, p.suffix="p", direction.suffix="logFC", sort=TRUE, decreasing=TRUE, nm=rownames(tab)){

  zList <- multi_pval2z(tab, prefix.v=prefix.v, p.suffix=p.suffix, direction.suffix=direction.suffix) %>%
    as.data.frame() %>%
    as.list() %>%
    lapply(FUN=\(x){
      y <- stats::setNames(x, nm=nm)
      y[!is.na(y)]
    })

  if(sort) {
    zList <- zList %>%
      lapply(FUN=sort, decreasing=decreasing)
  }
  return(zList)
}
