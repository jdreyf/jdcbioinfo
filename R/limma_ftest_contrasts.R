#' Apply limma's lmFit, contrasts.fit, & eBayes to one or more contrasts, perform moderated F-test, and return a table
#'
#' Apply \pkg{limma}'s \code{lmFit}, \code{contrasts.fit}, & \code{eBayes} to one or more contrasts,  perform moderated F-test, and return
#' a table. See examples in vignette.
#'
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams ezlimma::limma_cor
#' @return Data frame.
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means}
#' is \code{FALSE}.
#'
#' @references McCarthy DJ & Smyth GK (2009). Testing significance relative to a fold-change threshold is a TREAT.
#' Bioinformatics 25, 765-771.
#' @seealso \code{\link[limma]{lmFit}}; \code{\link[limma]{eBayes}}.
#' @export

limma_ftest_contrasts <- function(object, grp=NULL, contrast.v, add.means=!is.null(grp), weights=NA, design=NULL, prefix='',
                                 trend=FALSE, block=NULL, correlation=NULL){

  if (is.null(design)|add.means) stopifnot(ncol(object)==length(grp), colnames(object)==names(grp))

  # make model
  if (is.null(design)){
    design <- model.matrix(~0+grp)
    colnames(design) <- sub('grp', '', colnames(design), fixed=TRUE)
  }

  # lmFit
  if (!is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ cat('object$weights are being ignored\n') }
    fit <- lmFit(object, design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- lmFit(object, design, block = block, correlation = correlation)
  }

  # contrast
  contr.mat <- makeContrasts(contrasts=contrast.v, levels=design)

  # contrasts.fit & eBayes
  fit2 <- contrasts.fit(fit, contr.mat)
  fit2 <- eBayes(fit2, trend=trend)

  # topTable
  ttf <- topTable(fit2, number=Inf, adjust.method='BH', coef=contrast.v)

  if("F" %in% colnames(ttf)){
    ttf <- ttf[, c('F', 'P.Value', 'adj.P.Val')]
  }else if("t" %in% colnames(ttf)){
    ttf <- ttf[, c('t', 'P.Value', 'adj.P.Val')]
  }
  colnames(ttf) <- sub('P.Value', 'p', sub('adj.P.Val', 'FDR', colnames(ttf)))
  if (prefix!=''){ colnames(ttf) <- paste(prefix, colnames(ttf), sep='.') }

  #cbind grp means
  if (add.means){
    grp.means <- t(apply(object, 1, FUN=function(v) tapply(v, grp, mean, na.rm=TRUE)))
    colnames(grp.means) <- paste(colnames(grp.means), 'avg', sep='.')
    ttf <- cbind(grp.means[rownames(ttf),], ttf)
  }
  return(ttf)
}
