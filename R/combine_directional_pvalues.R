#' Directionally combine p-values
#'
#' Directionally combine two-sided test p-values.
#'
#' @inheritParams multi_pval2z
#' @inheritParams rankprod
#' @inheritParams ezlimma::combine_pvalues
#' @return Vector of p-values and/or direction and FDR .
#' @export

combine_directional_pvalues <- function(tab, prefix.v=NULL, p.suffix="p", direction.suffix="logFC", only.p=FALSE, prefix=NULL){

  mat.z <- multi_pval2z(tab, prefix.v=prefix.v, p.suffix=p.suffix, direction.suffix=direction.suffix)
  combz.v <- rowSums(mat.z)/sqrt(ncol(mat.z))
  combp.v <- stats::pnorm(combz.v)
  combp.v <- 2*pmin(combp.v, 1-combp.v)

  if(only.p){
    combp.v <- stats::setNames(combp.v, nm=rownames(tab))
    return(combp.v)
  }

  direction <- ifelse(combz.v>0, "Up", "Down")
  fdr <- stats::p.adjust(combp.v, method="BH")
  res <- data.frame(Combined.Direction=direction, Combined.p= combp.v, Combined.FDR=fdr)
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")

  return(res)
}
