#' Perform dream's moderated F-test for all pairwise comparisons, and return a table
#'
#' Make contrasts for all pairwise comparisons, and pass them to \code{\link[jdcbioinfo]{dream_ftest_contrasts}}.
#'
#' @inheritParams dream_ftest_contrasts
#' @return Data frame.
#' @export

dream_ftest_pairwise <- function(object, formula, pheno, grp, weights=NA, add.means=TRUE, moderated=TRUE, prefix=""){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp))

  # make contrast
  comb <- utils::combn(unique(grp), 2)
  contrast.v <- character(0)
  for(i in 1:ncol(comb)){
    contrast.v[paste0(comb[2, i], "vs", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  tt <- dream_ftest_contrasts(object=object, formula=formula, pheno=pheno, contrast.v=contrast.v, weights=weights, grp=grp, add.means=add.means, moderated=moderated, prefix=prefix)
  return(tt)
}
