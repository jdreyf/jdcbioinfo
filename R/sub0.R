#' Sub with empty replacement
#'
#' Sub with empty replacement.
#'
#' @inheritParams base::grep
#' @export

sub0 <- function(pattern, x, replacement="", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE){
  sub(pattern=pattern, replacement=replacement, x=x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}
