#' Convert p-values to Z-scores
#'
#' Counvert two-sided test p-values to Z-scores using \code{qnorm}.
#'
#' @param pval Vector of two-sided test p-values.
#' @param direction Numeric or character vector of same length of p-value vector. Examples of numeric values are log fold-change,
#' regression slope, or correlation coefficient. For character values, only "Up" and "Down" are allowed.
#' @return Vector of Z-scores.
#' @details
#' P-values of zero are changed to half the lowest non-zero p-value and p-values of one are analogously changed, as p-values of zero or one
#' would yield infinite z-scores.
#'
#' @export

pval2z <- function(pval, direction){
  stopifnot(length(pval)==length(direction), names(pval)==names(direction))
  if(any(pval<0 | pval>1, na.rm = TRUE)) stop("p-values should be within [0, 1]")

  if(is.numeric(direction)){
    direction <- sign(direction)
  }else{
    stopifnot(direction %in% c("Up", "Down"))
    direction <- ifelse(direction=="Up", 1, -1)
  }
  # clip 0 and 1
  LB <- min(pval[pval>0])/2
  UB <- 1-(1-max(pval[pval<1]))/2
  pval <- pmax(LB, pmin(UB, pval))

  z <- stats::qnorm(pval/2, lower.tail=FALSE) * direction
  return(z)
}
