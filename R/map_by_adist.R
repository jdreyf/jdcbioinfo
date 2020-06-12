#' Map two character vectors by approximate string distances
#'
#' Map two character vectors by generalized Levenshtein (edit) distance i.e. approximate string distances.
#' @param x	a character vector.
#' @param y	a character vector.
#' @param ignore.case	a logical. If TRUE, case is ignored for computing the distances.
#' @param max.norm.dist maximal normalized distance between 0 and 1.
#' @return Two-column Data frame of which each row is the mapped elements in x and y.
#' @export

map_by_adist <- function(x, y, ignore.case=TRUE, max.norm.dist=0.2) {

  l1 <- length(x)
  l2 <- length(y)

  if (l1 > l2) {
    z <- x
    x <- y
    y <- z
    rm(z)
  }

  mat.dist <- utils::adist(x=x, y=y, ignore.case=ignore.case)
  dimnames(mat.dist) <- list(x, y)

  nc1 <- nchar(x)
  nc2 <- nchar(y)

  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      mat.dist[i, j] <- mat.dist[i, j]/max(nc1[i], nc2[j])
    }
  }

  mat.dist <- mat.dist[rowSums(mat.dist<=max.norm.dist)>0, ]

  idxs <- apply(mat.dist, 1, which.min)

  x <- rownames(mat.dist)
  y <- colnames(mat.dist)[idxs]

  if (l1 > l2) {
    tab <- data.frame(x=y, y=x)
  } else {
    tab <- data.frame(x=x, y=y)
  }

  return(tab)
}
