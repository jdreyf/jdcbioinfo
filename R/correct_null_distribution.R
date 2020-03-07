#' Correct null distribution
#'
#' Correct null distribution using  \pkg{fdrtool}.
#'
#' @param p.suffix The suffix of p-value columns.
#' @param fdr.suffix The suffix of FDR columns.
#' @param verbose Print out status messages.
#' @inheritParams multi_pval2z
#' @return Matrix or data.fame with corrected p-values and the corresponding FDR.
#' @export

correct_null_distribution <- function(tab, prefix.v=NULL, p.suffix="p", fdr.suffix="FDR", direction.suffix="logFC", verbose=TRUE){

  stopifnot(length(p.suffix)==1, length(fdr.suffix)==1, length(direction.suffix)==1)

  mat_z <- multi_pval2z(tab, prefix.v=prefix.v, p.suffix=p.suffix, direction.suffix=direction.suffix)
  res <- apply(mat_z, MARGIN=2, FUN=fdrtool::fdrtool, statistic="normal", plot=FALSE, verbose=verbose)

  for(nm in names(res)){
    tab[, paste(nm,  p.suffix, sep=".")] <- res[[nm]]$pval
    tab[, paste(nm,  fdr.suffix, sep=".")] <- stats::p.adjust(res[[nm]]$pval, method="BH")
  }
  return(tab)
}
