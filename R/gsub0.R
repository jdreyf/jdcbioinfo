#' Gsub with empty replacement
#'
#' Gsub with empty replacement.
#'
#' @inheritParams base::grep
#' @export

gsub0 <- function(pattern, x, replacement="", ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE){
  gsub(pattern=pattern, replacement=replacement, x=x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes)
}
