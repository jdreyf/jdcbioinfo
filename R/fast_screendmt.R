#' Post screenDMT
#'
#' Post screenDMT i.e. combining replication & reversed replication test results
#'
#' @import DirectionalMaxPTest
#' @param direction Character, either "same", "opposite", or default "both"
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

fast_screendmt <- function(tab, cols=1:4, direction = c("both", "same", "opposite"), reorder.rows=FALSE,
                           p.adj.rate=c("FDR", "FWER"), prefix=NULL, keep.input=FALSE){
  direction <- match.arg(direction,  c("both", "same", "opposite"))
  p.adj.rate <- match.arg(p.adj.rate, c("FDR", "FWER"))
  stopifnot(length(cols) == 4)
  stopifnot(DirectionalMaxPTest::check_tab(tab, num.cols = cols))
  stat.cols <- cols[c(1, 3)]
  p.cols <- cols[c(2, 4)]
  stopifnot(0 <= tab[, p.cols], tab[, p.cols] <= 1, !is.na(tab[, cols]),
            is.logical(reorder.rows), is.logical(keep.input), !is.na(prefix),
            !is.null(rownames(tab)), nrow(tab) > 1)

  sign.mat <- apply(tab[, stat.cols], MARGIN = 2, FUN = sign)
  sign.sum <- as.character(rowSums(sign.mat))
  direction.v <- sapply(sign.sum, FUN = switch, "2" = "Up", "-2" = "Down", "0" = "Opposite")

  pval.tab <- tab[, p.cols]
  maxp <-  pmax(pval.tab[, 1], pval.tab[, 2], na.rm = TRUE)
  maxp[is.infinite(maxp)] <- NA
  minp <-  pmin(pval.tab[, 1], pval.tab[, 2], na.rm = TRUE)
  minp[is.infinite(minp)] <- NA

  direstion.list <- list(both = c("same", "opposite"), same = "same", opposite = "opposite")

  res.list <- list()
  for ( d in direstion.list[[direction]]) {
    pval <- rep(1, length(maxp))

    if (d == "same") {
      direction.bool <- direction.v != "Opposite"
    } else {
      direction.bool <- direction.v == "Opposite"
    }

    pval[direction.bool] <- 0.5 * maxp[direction.bool]
    chisq <- stats::qchisq(pval, df = 1, lower.tail = FALSE)

    res <- data.frame(direction = direction.v, chisq = chisq, p = pval, row.names = rownames(tab)) %>%
      dplyr::mutate(max2 = pmax(pval, minp, na.rm = TRUE),
                    max2 = replace(max2, is.infinite(max2), NA)) %>%
      dplyr::arrange(max2) %>%
      dplyr::mutate(rnk = dplyr::row_number(),
                    adj.num = findInterval(max2, sort(minp, na.last = TRUE), rightmost.closed = TRUE))

    if (p.adj.rate=="FDR") {
      res <- res %>%
        dplyr::mutate(bh.point = adj.num * max2/rnk,
                      FDR = rev(cummin(rev(bh.point))))
    } else {
      res <- res %>%
        dplyr::mutate(bh.point = adj.num * max2,
                      FWER = pmin(rev(cummin(rev(bh.point))), 1))
    }
    res.list[[d]] <- res %>%
      dplyr::select(!(max2:bh.point))
  }

  if(direction != "both") {
    res <- res.list[[direction]]
  } else {
    res <- rbind(res.list[["same"]] %>% dplyr::filter(direction != "Opposite"),
                 res.list[["opposite"]] %>% dplyr::filter(direction == "Opposite"))
  }

  res <- res[rownames(tab), ]

  if (keep.input) res <- data.frame(res, tab[, cols, drop=FALSE])
  if (reorder.rows) res <- res %>% dplyr::arrange(dplyr::across(dplyr::matches("F(D|WE)R")), p)
  if (!is.null(prefix)) colnames(res)[1:4] <- paste(prefix, colnames(res)[1:4], sep=".")

  return(res)
}
