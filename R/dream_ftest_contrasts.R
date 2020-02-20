#' Apply variancePartition's dream to one or more contrasts, perform moderated F-test, and return a table using limma's topTable
#'
#' Apply \pkg{variancePartition}'s \code{dream}, to one or more contrasts,  perform perform moderated F-test, and return
#' a table using \pkg{limma}'s \code{topTable}.
#'
#' @param pheno data.frame with columns corresponding to formula
#' @inheritParams limma_ftest_contrasts
#' @inheritParams dream_contrasts
#' @return Data frame.
#' @export

dream_ftest_contrasts <- function(object, formula, pheno, contrast.v, weights=NA, grp=NULL, add.means=!is.null(grp), prefix=""){

  stopifnot(is.na(weights) || is.null(weights) || dim(weights)==dim(object) || length(weights)==nrow(object) ||
              length(weights)==ncol(object))
  stopifnot(!(is.null(grp) & add.means))

  if (is.vector(object)) stop("'object' must be a matrix-like object; you can coerce it to one with 'as.matrix()'")
  if (any(duplicated(rownames(object)))) stop("object cannot have duplicated rownames.")
  if (any(rownames(object)=="")) stop("object cannot have an empty rowname ''.")

  contrast.v <- strsplit(contrast.v, split=" *- *")
  stopifnot(sapply(contrast.v, length) %in% c(1,2))
  L <- sapply(contrast.v, FUN=function(contr) {
    variancePartition::getContrast(object, formula=formula, data=pheno, coefficient=contr)})
  if (is.vector(L)) {
    L <- as.matrix(L)
    colnames(L) <- names(contrast.v)
  }

  # can't set weights=NULL in lmFit when using voom, since lmFit only assigns
  # weights "if (missing(weights) && !is.null(y$weights))"
  # can't make this into separate function, since then !missing(weights)
  # length(NULL)=0; other weights should have length > 1

if (length(weights)!=1 || !is.na(weights)){
  if (!is.matrix(object) && !is.null(object$weights)){ warning("object$weights are being ignored") }
    fit <- variancePartition::dream(object, formula=formula, data=pheno, L=L, weights=weights)
  } else {
    fit <- variancePartition::dream(object, formula=formula, data=pheno, L=L)
  }

  tt <- limma::topTable(fit, coef=colnames(L), number=Inf , sort.by="none")

  if("F" %in% colnames(tt)){
    tt <- tt[, c("F", "P.Value", "adj.P.Val")]
  } else if("t" %in% colnames(tt)){
    tt <- tt[, c("t", "P.Value", "adj.P.Val")]
  }

  colnames(tt) <- gsub("P\\.Value", "p", gsub("adj\\.P\\.Val", "FDR", colnames(tt)))
  tt <- tt[order(tt$p), ]
  if (prefix!=""){ colnames(tt) <- paste(prefix, colnames(tt), sep=".") }

  if (add.means) {
    mat <- as.matrix(object)
    grps <- unique(sort(grp))
    mat.avg <- sapply(grps, FUN=function(g) rowMeans(mat[, grp==g], na.rm=TRUE))
    colnames(mat.avg) <- paste(grps, "avg", sep=".")
    tt <- cbind(mat.avg[rownames(tt), ], tt)
  }
  return(tt)
}
