#' Perform limma's moderated F-test for all pairwise comparions, and return a table
#'
#' Make contrasts for all pairwise comparions, and pass them to code{\link[jdcbioinfo]{limma_ftest_contrasts}}.
#'
#' @inheritParams limma_ftest_contrasts
#' @return Data frame.
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means}
#' is \code{FALSE}.
#'
#' @references McCarthy DJ & Smyth GK (2009). Testing significance relative to a fold-change threshold is a TREAT.
#' Bioinformatics 25, 765-771.
#' @export

limma_ftest_pairwise <- function(object, grp, add.means=TRUE, weights=NA, design=NULL, prefix='',
                                 trend=FALSE, block=NULL, correlation=NULL){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp))
  if(!is.null(design)) stopifnot(grp %in% colnames(design))

  # make contrast
  comb <- combn(unique(grp), 2)
  contrast.v <- character(0)
  for(i in 1:ncol(comb)){
    contrast.v[i] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  ttf <- limma_ftest_contrasts(object, grp=grp, contrast.v=contrast.v, add.means=add.means, weights=weights, design=design, prefix=prefix,
                               trend=trend, block=block, correlation=correlation)
  return(ttf)
}
