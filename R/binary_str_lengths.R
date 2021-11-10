#' Calculate lengths of each string of ones in a binary vector for CGM analyses
#'
#' Calculate lengths of each string of ones in a binary vector for CGM analyses.
#' If none are found, a numeric vector of length zero is returned.
#'
#' @param vv Numeric vector of ones and zeroes only.
#' @export

binary_str_lengths <- function(vv){
  stopifnot(is.numeric(vv), vv %in% 0:1)
  xx <- paste(vv, collapse = "")
  yy <- unlist(strsplit(xx, split="0"))
  sort(nchar(yy[yy!=""]))
}
