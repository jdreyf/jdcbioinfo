#' Match names of chemicals
#'
#' Match names of chemicals, which are often common, non-systematic names.
#'
#' @inheritParams base::match
#' @return Vector of the positions of (first) matches of its first argument in its second.
#' @export

match_met_nms <- function(x, table){
  stopifnot(dim(x)==0, dim(table)==0, length(x) > 0, length(table) > 0)

  xx <- tolower(gsub0("-|\\(|\\)", x))
  tble <- tolower(gsub0("-|\\(|\\)", table))

  tble2 <- apply(as.matrix(tble), 1, FUN=function(xx){
    if (grepl("ate$", x = xx)){
      sub("ate$", "ic acid", xx)
    } else if (grepl("ic acid$", x = xx)){
      sub("ic acid$", "ate", xx)
    } else xx
  })

  res <- match(xx, tble)
  res2 <- match(xx, tble2)

  res[is.na(res)] <- res2[is.na(res)]
  res
}
