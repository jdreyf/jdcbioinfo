#' Two-proportion z-test (pooled) for non-zero null
#'
#' Two-proportion z-test (pooled) for non-zero null.
#'
#' @param x1 number of successes in group 1.
#' @param x2 number of successes in group 2.
#' @param n1 total number of games in group 1.
#' @param n2 total number of games in group 2.
#' @param delta null hypothesis of the difference in proportion.
#' @inheritParams stats::prop.test
#' @return a number of p-value.
#' @export

z_2prop_pooled <- function(x1, x2, n1, n2, delta=0, alternative=c("two.sided", "less", "greater"), correct=TRUE){

  stopifnot(x1<=n1, x2<=n2, x1>=0, x2>=0, n1>0, n2>0)
  alternative <- match.arg(alternative)

  # get cc
  if (correct) {
    cc <- (1/n1+1/n2)/2
    cc <- ifelse(alternative=="greater", cc, -cc)
  } else {
    cc <- 0
  }

  # get z-score
  p.diff <- (x1/n1)-(x2/n2)-delta+cc
  p.common <- (x1+x2)/(n1+n2)
  se <- sqrt(p.common * (1-p.common) * (1/n1+1/n2))
  if (se==0) {
    z.score <- NaN
  } else {
    z.score <- p.diff/se
  }

  # get p-value
  if (alternative=="two.sided") {
    pval <- stats::pnorm(z.score, lower.tail=TRUE)
    pval <- 2*min(pval, 1-pval)
    if (pval>1) pval <- 1

  } else if (alternative=="less") {
    pval <- stats::pnorm(z.score, lower.tail=TRUE)

  } else if (alternative=="greater") {
    pval <- stats::pnorm(z.score, lower.tail=FALSE)
  }

  return(unname(pval))
}
