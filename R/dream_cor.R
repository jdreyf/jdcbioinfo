#' Apply variancePartition's dream to one or more contrasts, perform moderated t-test, and return a table using limma's topTable
#'
#' Apply \pkg{variancePartition}'s \code{dream}, to one or more contrasts,  perform moderated t-test, and return
#' a table using \pkg{limma}'s \code{topTable}.
#'
#' @param pheno data.frame with columns corresponding to formula
#' @param coef column name of the mixed model
#' @param moderated Logical; should variancePartition::eBayes be used?
#' @inheritParams ezlimma::limma_cor
#' @inheritParams variancePartition::dream
#' @return Data frame.
#' @export

dream_cor <- function(object, formula, pheno, weights=NA, coef="", moderated=TRUE, cols=c("P.Value", "adj.P.Val", "logFC")){

  if (!requireNamespace("BiocParallel", quietly = TRUE)) stop("Package \"BiocParallel\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("variancePartition", quietly = TRUE)) stop("Package \"variancePartition\" must be installed to use this function.", call. = FALSE)

  stopifnot(is.na(weights) || is.null(weights) || dim(weights)==dim(object) || length(weights)==nrow(object) || length(weights)==ncol(object))

  if (is.vector(object)) stop("'object' must be a matrix-like object; you can coerce it to one with 'as.matrix()'")
  if (any(duplicated(rownames(object)))) stop("object cannot have duplicated rownames.")
  if (any(rownames(object)=="")) stop("object cannot have an empty rowname ''.")

  # can't set weights=NULL in lmFit when using voom, since lmFit only assigns
  # weights "if (missing(weights) && !is.null(y$weights))"
  # can't make this into separate function, since then !missing(weights)
  # length(NULL)=0; other weights should have length > 1

  bp <- BiocParallel::SerialParam(progressbar=TRUE)
  BiocParallel::register(BiocParallel::bpstart(bp))

  if (length(weights)!=1 || !is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ warning("object$weights are being ignored") }
      fit <- variancePartition::dream(object, formula=formula, data=pheno, weights=weights, BPPARAM=bp)
    } else {
      fit <- variancePartition::dream(object, formula=formula, data=pheno, BPPARAM=bp)
    }
  if (moderated) {
    fit <- variancePartition::eBayes(fit)
  }
  BiocParallel::bpstop(bp)

  tt <- variancePartition::topTable(fit, coef=coef, number=Inf , sort.by="p")[, cols]
  colnames(tt) <- gsub("P\\.Value", "p", gsub("adj\\.P\\.Val", "FDR", gsub("logFC", "slope", colnames(tt))))
  colnames(tt) <- paste(coef, colnames(tt), sep=".")

  return(tt)
}
