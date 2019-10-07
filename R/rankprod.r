#' Rank products of a matrix where larger statistics have stronger rank
#'
#' Rank products of a matrix. Larger statistics have stronger rank.
#'
#' @param mat Numeric matrix with columns holding statistics per comparison & rows are analytes.
#' @param nsim Numeric number of simulations.
#' @param seed Random seed.
#' @inheritParams ezlimma::ezcor
#' @return Data frame with statistics from rank products test.
#' @export

rankprod <- function(mat, nsim=1e7-1, reorder.rows=TRUE, prefix=NULL, seed=100){
  stopifnot(!is.null(colnames(mat)))
  rmat <- apply(mat, 2, rank)
  rmat <- rmat/nrow(mat)

  colnames(rmat) <- paste(gsub("\\..$", "", colnames(mat)), "Rank.Prop", sep=".")

  set.seed(seed)
  rmat.sim <- apply(rmat, 2, function(v, nsim) sample(v, size=nsim, replace=TRUE), nsim)

  rankprod <- -1 * apply(rmat, 1, prod)
  rankprod.sim <- -1 * apply(rmat.sim, 1, prod)

  Fn <- stats::ecdf(c(rankprod.sim, Inf))
  pval <- Fn(rankprod)
  fdr <- stats::p.adjust(pval, method="BH")

  res <- data.frame(rmat, Rankprod.p=pval, Rankprod.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Rankprod.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  return(res)
}
