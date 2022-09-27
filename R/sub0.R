#' Sub with empty replacement
#'
#' Sub with empty replacement.
#'
#' @param x character vector where matches are sought or an object which can be coerced by `as.character` to a character vector.
#' @inheritParams base::sub
#' @export

sub0 <- function(pattern, x, replacement="", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE){
  sub(pattern=pattern, replacement=replacement, x=x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}
