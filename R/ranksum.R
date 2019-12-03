#' Rank sum of a matrix where larger statistics have stronger rank
#'
#' Rank sum of a matrix. Larger statistics have stronger rank.
#'
#' @inheritParams rankprod
#' @return Data frame with statistics from rank sum test.
#' @export

ranksum <- function(mat, nsim=1e7-1, reorder.rows=TRUE, prefix=NULL, seed=100){

  stopifnot(ncol(mat) > 1,  !is.null(colnames(mat)))
  if(nsim > 1e7-1) stop("nsim too large to have enough precision")

  rmat <- apply(mat, 2, rank)
  rmat <- rmat/nrow(mat)
  colnames(rmat) <- paste(gsub("\\..$", "", colnames(mat)), "Rank.Prop", sep=".")

  set.seed(seed)
  rmat.sim <- apply(rmat, 2, function(v, nsim) sample(v, size=nsim, replace=TRUE), nsim)

  ranksum <- rowSums(rmat)
  ranksum.sim <- rowSums(rmat.sim)

  Fn <- stats::ecdf(c(ranksum.sim, Inf))
  pval <- 1 - Fn(ranksum)
  fdr <- stats::p.adjust(pval, method="BH")

  res <- data.frame(rmat, Ranksum.p=pval, Ranksum.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Ranksum.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")

  return(res)
}
