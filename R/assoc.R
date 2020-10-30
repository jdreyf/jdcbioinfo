#' Association test between two variables
#'
#' Association test for between variables.
#'
#' @param v1 Variable 1.
#' @param v2 Variable 2.
#' @param seed Random seed.
#' @return p-value.

assoc <- function(v1, v2, seed=100) {
  stopifnot(length(v1) == length(v2))

  keep1 <- !(is.na(v1) | v1 %in% c("", "NA"))
  keep2 <- !(is.na(v2) | v2 %in% c("", "NA"))
  keep <- keep1 & keep2
  v1 <- v1[keep]
  v2 <- v2[keep]

  is.num1 <- is.numeric(v1)
  is.num2 <- is.numeric(v2)
  is.num <- c(is.num1, is.num2)

  if (!any(is.num)) {
    if (length(v1) >= 50) {
      pval <- stats::chisq.test(table(v1, v2))$p.value
    } else {
      set.seed(seed)
      pval <- stats::chisq.test(table(v1, v2), simulate.p.value=TRUE)$p.value
    }
  } else {
    if (all(is.num)) {
      pval <- stats::cor.test(x=v1, y=v2)$p.value
    } else {
      if (is.num1) {
        pval <- summary(stats::aov(v1~v2))[[1]]["v2", "Pr(>F)"]
      } else {
        pval <- summary(stats::aov(v2~v1))[[1]]["v1", "Pr(>F)"]
      }
    }
  }
  return(pval)
}
