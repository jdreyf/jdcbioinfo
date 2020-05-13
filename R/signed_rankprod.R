#' Signed rank products of a two column matrix where larger statistics have stronger rank in the same or opposite direction
#'
#' Signed rank products of a two column matrix. Larger statistics have stronger rank in the same or opposite direction.
#'
#' @param mat Numeric matrix with two columns holding statistics per comparison & rows are analytes.
#' @param same.dirction Logical indicates whether the two ranks should be in the same direction.
#' @inheritParams rankprod
#' @return Data frame with statistics from signed rank products test.
#' @export

signed_rankprod <- function(mat, nsim=1e7-2, same.dirction=FALSE, reorder.rows=TRUE, prefix=NULL, seed=100){

  stopifnot(ncol(mat)==2,  !is.null(colnames(mat)))
  if(nsim > 1e7-2) stop("nsim too large to have enough precision")

  rmat <- apply(mat, 2, function(v) {
    r <- rank(abs(v))
    r[v < 0] <- -r[v < 0]
    return(r)
  })

  rmat <- rmat/nrow(mat)
  colnames(rmat) <- paste(gsub("\\..$", "", colnames(mat)), "Signed.Rank.Prop", sep=".")

  set.seed(seed)
  rmat.sim <- apply(rmat, 2, function(v, nsim) sample(v, size=nsim, replace=TRUE), nsim)

  rankprod <- apply(rmat, 1, prod)
  rankprod.sim <- apply(rmat.sim, 1, prod)

  Fn <- stats::ecdf(c(rankprod.sim, Inf, -Inf))
  pval <- Fn(rankprod)

  if(same.dirction) {
    pval <- 1 - pval
  } else {
    pval <- 2 * pmin(pval, 1 - pval)
  }

  fdr <- stats::p.adjust(pval, method="BH")

  direction <- rep("", nrow(rmat))
  direction[rankprod < 0] <- "Opposite"
  direction[rmat[, 1] > 0 & rmat[, 2] > 0] <- "Up"
  direction[rmat[, 1] < 0 & rmat[, 2] < 0] <- "Down"

  res <- data.frame(rmat, Direction=direction, Signed.Rankprod.p=pval, Signed.Rankprod.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Signed.Rankprod.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  return(res)
}
