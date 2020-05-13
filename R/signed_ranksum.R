#' Signed rank sum of a matrix where larger statistics have stronger rank in the same direction
#'
#' Signed rank sum of a matrix. Larger statistics have stronger rank in the same direction.
#'
#' @inheritParams rankprod
#' @return Data frame with statistics from signed rank sum test.
#' @export

signed_ranksum <- function(mat, nsim=1e7-2, reorder.rows=TRUE, prefix=NULL, seed=100){

  stopifnot(ncol(mat) > 1,  !is.null(colnames(mat)))
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

  ranksum <- rowSums(rmat)
  ranksum.sim <- rowSums(rmat.sim)

  Fn <- stats::ecdf(c(ranksum.sim, Inf, -Inf))
  pval <- Fn(ranksum)

  pval <- 2 * pmin(pval, 1 - pval)
  fdr <- stats::p.adjust(pval, method="BH")
  direction <- ifelse(ranksum > 0, "Up", "Down")

  res <- data.frame(rmat, Direction=direction, Signed.Ranksum.p=pval, Signed.Ranksum.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Signed.Ranksum.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  return(res)
}
