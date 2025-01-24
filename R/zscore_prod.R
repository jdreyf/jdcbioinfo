#' Z-score products of a two column matrix where larger scores in the same or opposite direction
#'
#' Z-score products of a two column matrix where larger scores in the same or opposite directio.
#'
#' @inheritParams signed_rankprod
#' @return Data frame with statistics from Z-score products test.
#' @export

zscore_prod <- function(mat, nsim=1e7-2, same.dirction=FALSE, reorder.rows=TRUE, prefix=NULL, seed=100){

  stopifnot(ncol(mat)==2,  !is.null(colnames(mat)))
  if(nsim > 1e7-2) stop("nsim too large to have enough precision")

  colnames(mat) <- paste(gsub("\\..$", "", colnames(mat)), "z", sep=".")

  set.seed(seed)
  mat.sim <- apply(mat, 2, function(v, nsim) sample(v, size=nsim, replace=TRUE), nsim)

  rankprod <- apply(mat, 1, prod)
  rankprod.sim <- apply(mat.sim, 1, prod)

  Fn <- stats::ecdf(c(rankprod.sim, Inf, -Inf))
  pval <- Fn(rankprod)

  if(same.dirction) {
    pval <- 1 - pval
  } else {
    pval <- 2 * pmin(pval, 1 - pval)
  }

  fdr <- stats::p.adjust(pval, method="BH")

  direction <- rep("", nrow(mat))
  direction[rankprod < 0] <- "Opposite"
  direction[mat[, 1] > 0 & mat[, 2] > 0] <- "Up"
  direction[mat[, 1] < 0 & mat[, 2] < 0] <- "Down"

  res <- data.frame(mat, Direction=direction, Prod.p=pval, Prod.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Prod.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  return(res)
}
