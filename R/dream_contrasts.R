#' Apply variancePartition's dream to one or more contrasts, perform moderated t-test, and return a table using limma's topTable
#'
#' Apply \pkg{variancePartition}'s \code{dream}, to one or more contrasts,  perform moderated t-test, and return
#' a table using \pkg{limma}'s \code{topTable}.
#'
#' @param pheno data.frame with columns corresponding to formula
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams variancePartition::dream
#' @return Data frame.
#' @export

dream_contrasts <- function(object, formula, pheno, contrast.v, weights=NA, grp=NULL, add.means=!is.null(grp), cols=c("P.Value", "adj.P.Val", "logFC")){

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

  mtt <- list()
  for (contr.nm in colnames(L)) {
   tt <- limma::topTable(fit, coef=contr.nm, number=Inf , sort.by="none")[, cols]
   tt$FC <- sign(tt$logFC)*2^(abs(tt$logFC))
   colnames(tt) <- gsub("P.Value", "p", gsub("adj.P.Val", "FDR", colnames(tt)))
   colnames(tt) <- paste(contr.nm, colnames(tt), sep=".")
   mtt[[contr.nm]] <- tt
  }
  mtt <- Reduce(cbind, mtt)
  mtt <- mtt[order(ezlimma::combine_pvalues(mtt)), ]

  if (add.means) {
    mat.avg <- sapply(unique(sort(grp)), FUN=function(g) rowMeans(object[, grp==g]))
    mtt <- cbind(mat.avg[rownames(mtt), ], mtt)
  }
  return(mtt)
}
