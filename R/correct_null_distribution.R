#' Correct null distribution
#'
#' Correct null distribution using  \pkg{fdrtool}
#'
#' @inheritParams multi_pval2z
#' @param p.suffix The suffix of p-value columns
#' @return Matrix or data.fame with corrected p-values and the corresponding FDR.
#' @export
#'
correct_null_distribution <- function(tab, prefix.v=NULL, p.suffix="p", fdr.suffix="FDR", direction.suffix="logFC"){

  stopifnot(length(p.suffix)==1, length( fdr.suffix)==1, length(direction.suffix)==1)

  mat_z <- multi_pval2z(tab, prefix.v=prefix.v, p.suffix=p.suffix, direction.suffix=direction.suffix)
  res <- apply(mat_z, 2, fdrtool::fdrtool, statistic="normal", plot=FALSE)

  for(prefix in prefix.v){
    tab[, paste(prefix,  p.suffix, sep=".")] <- res[[prefix]]$pval
    tab[, paste(prefix,  fdr.suffix, sep=".")] <- p.adjust(res[[prefix]]$pval, method="BH")
  }
  return(tab)
}
