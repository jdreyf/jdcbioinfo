#' Convert multiple columns of p-values to Z-scores
#'
#' Convert multiple columns (from a table-like object) of p-values from two-sided tests to Z-scores.
#'
#' @param p.suffix The suffix of p-value columns, i.e. the character string after the period "."
#' @param direction.suffix The suffix of direction columns such as log fold-change, regression slope, correlation
#' coefficient, or the "Up/Down" values
#' @inheritParams ezlimmaplot::ezvenn
#' @return Matrix of Z-scores.
#' @export

multi_pval2z <- function(tab, prefix.v=NULL, p.suffix="p", direction.suffix="logFC"){

  stopifnot(length(p.suffix)==1, length(direction.suffix)==1)

  if (!is.null(prefix.v)){
    p.cols <- paste(prefix.v, p.suffix, sep = ".")
    stopifnot(p.cols %in% colnames(tab))

  } else {
    p.suffix <- paste0("\\.", p.suffix, "$")
    p.cols <- grep(p.suffix, colnames(tab), value = TRUE)
    prefix.v <- gsub(p.suffix, "", p.cols)
  }

  direction.cols <- paste(prefix.v, direction.suffix, sep = ".")
  stopifnot(direction.cols %in% colnames(tab))

  res <- vapply(seq_along(prefix.v), FUN=function(i) {pval2z(tab[,p.cols[i]], direction=tab[,direction.cols[i]])}, FUN.VALUE=numeric(nrow(tab)))
  dimnames(res) <-  list(rownames(tab), prefix.v)

  return(res)
}
