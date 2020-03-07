#' Apply limma's lmFit, contrasts.fit, & treat to one or more contrasts, and return a table
#'
#' Apply \pkg{limma}'s \code{lmFit}, \code{contrasts.fit}, & \code{treat} to one or more contrasts,  and return
#' a table.
#'
#' @inheritParams rankprod
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

limma_treat_max <- function(object, grp=NULL, contrast.v, treat.lfc=log2(1.2), add.means=!is.null(grp), weights=NA, design=NULL, prefix='',
                                 trend=FALSE, block=NULL, correlation=NULL){

  if (is.null(design)|add.means) stopifnot(ncol(object)==length(grp), colnames(object)==names(grp))

  # make model
  if (is.null(design)){
    design <- stats::model.matrix(~0+grp)
    colnames(design) <- sub('grp', '', colnames(design), fixed=TRUE)
  }

  # lmFit
  if (!is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ cat('object$weights are being ignored\n') }
    fit <- limma::lmFit(object, design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design, block = block, correlation = correlation)
  }

  # contrast
  contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=design)

  # contrasts.fit & treat
  fit2 <- limma::contrasts.fit(fit, contr.mat)
  fit2 <- limma::treat(fit2, lfc=treat.lfc, trend=trend)

  # topTable
  if (length(contrast.v) == 1) {
    ttf <- limma::topTreat(fit2, number=Inf, adjust.method='BH', coef=contrast.v, sort.by="p")
    ttf <- ttf[, c('P.Value', 'adj.P.Val')]
    colnames(ttf) <- sub('P.Value', 'p', sub('adj.P.Val', 'FDR', colnames(ttf)))
  } else {
    p.min <- apply(fit2$p.value, MARGIN=1, FUN=min)
    fdr <- stats::p.adjust(p.min, method="BH")
    ttf <- data.frame(p = p.min, FDR = fdr)
    ttf <- ttf[order(ttf$p), ]
  }

  if (prefix!=''){ colnames(ttf) <- paste(prefix, colnames(ttf), sep='.') }

  #cbind grp means
  if (add.means){
    grp.means <- t(apply(object, 1, FUN=function(v) tapply(v, grp, mean, na.rm=TRUE)))
    colnames(grp.means) <- paste(colnames(grp.means), 'avg', sep='.')
    ttf <- data.frame(grp.means[rownames(ttf),], ttf)
  }
  return(ttf)
}
