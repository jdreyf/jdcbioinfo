#' Post screenDMT
#'
#' Post screenDMT i.e. combining replication & reversed replication test results
#'
#' @import DirectionalMaxPTest
#' @param keep.input Logical, should returned data frame include `tab[, cols]`? If so, `tab` must have column names and these
#' cannot include "chisq", "p", and "FDR" or "FWER" (depending on `p.adj.rate`)
#' @inheritParams DirectionalMaxPTest::screendmt
#' @return Data frame whose rows correspond to the rows of `tab` with the same row names and whose columns include
#' \describe{
#' \item{direction}{diretion of change, up, down, or opposite}
#' \item{chisq}{Chi-square on 1 degree of freedom}
#' \item{p}{P-value}
#' \item{FDR or FWER}{Adjusted p-value}
#' }
#' @details Larger chi-square values are more significant. `tab` must have more than one row.
#' @md
#' @export

post_screendmt <- function(tab, cols=1:4, reorder.rows=FALSE, p.adj.rate=c("FDR", "FWER"), prefix=NULL, keep.input=FALSE){
  p.adj.rate <- match.arg(p.adj.rate, c("FDR", "FWER"))
  prod.sgns=c(same=1, opposite=-1)

  resList <- list()
  for (nm in names(prod.sgns)) {
    resList[[nm]] <- DirectionalMaxPTest::screendmt(tab, cols=cols, prod.sgn=prod.sgns[nm], reorder.rows=FALSE,
                                                    p.adj.rate=p.adj.rate, prefix=NULL, keep.input=FALSE)
  }

  stat.cols <- cols[c(1, 3)]
  sign.mat <- apply(tab[, stat.cols], MARGIN = 2, FUN = sign)
  sign.sum <- apply(sign.mat, MARGIN = 1, FUN = sum)
  direction.v <- switch (sign.sum, `2` = "up", `-2` = "down", `0` = "opposite")

  res <- resList[["same"]]
  res[direction.v == "opposite", ] <- resList[["opposite"]][direction.v == "opposite", ]
  res <- cbind(direction = direction.v, res)

  if (keep.input) res <- data.frame(res, tab[, cols, drop=FALSE])
  if (reorder.rows) res <- res |> dplyr::arrange(dplyr::across(dplyr::matches("F(D|WE)R")), p)
  if (!is.null(prefix)) colnames(res)[1:4] <- paste(prefix, colnames(res)[1:4], sep=".")

  return(res)
}
