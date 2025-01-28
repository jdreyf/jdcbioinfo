#' Z-score products of a two column matrix where larger scores in the same or opposite direction
#'
#' Z-score products of a two column matrix where larger scores in the same or opposite direction.
#'
#' @param direction Character, either "same", "opposite", or default "both"
#' @inheritParams signed_rankprod
#' @return Data frame with statistics from Z-score products test.
#' @export

zscore_prod <- function(mat, nsim=1e7-2, direction = c("both", "same", "opposite"), reorder.rows=TRUE, prefix=NULL, seed=100){

  direction <- match.arg(direction,  c("both", "same", "opposite"))
  stopifnot(ncol(mat)==2,  !is.null(colnames(mat)))
  if(nsim > 1e7-2) stop("nsim too large to have enough precision")

  colnames(mat) <- paste(gsub("\\..$", "", colnames(mat)), "z", sep=".")

  set.seed(seed)
  mat.sim <- apply(mat, 2, function(v, nsim) sample(v, size=nsim, replace=TRUE), nsim)

  prod <- apply(mat, 1, prod)
  prod.sim <- apply(mat.sim, 1, prod)

  Fn <- stats::ecdf(c(prod.sim, Inf, -Inf))
  pval <- Fn(prod)

  if(direction=="same") {
    pval <- 1 - pval
  } else if(direction=="both") {
    pval <- 2 * pmin(pval, 1 - pval)
  }

  fdr <- stats::p.adjust(pval, method="BH")

  sign.mat <- apply(mat, MARGIN = 2, FUN = sign)
  sign.sum <- as.character(rowSums(sign.mat))
  direction.v <- sapply(sign.sum, FUN = switch, "2" = "Up", "-2" = "Down", "0" = "Opposite")

  res <- data.frame(mat, Direction=direction.v, Prod.p=pval, Prod.FDR=fdr, row.names=rownames(mat))
  if(reorder.rows) res <- res[order(res$Prod.p), ]
  if(!is.null(prefix)) colnames(res) <- paste(prefix, colnames(res), sep=".")
  return(res)
}
