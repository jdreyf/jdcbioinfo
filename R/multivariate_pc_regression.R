#' Multivariate regression for meta variables and principal components
#'
#' Multivariate regression for meta variables and principal components
#'
#' @inheritParams pc_assoc
#' @return Matrix of p-values.
#' @export

multivariate_pc_regression <- function(object, pheno.df, metavars, num.pc=10, plot=TRUE, name="assoc_with_pc_pvals", width=8, height=7, xrot=0, yrot=0) {
  stopifnot(limma::isNumeric(object), nrow(object[rowSums(is.na(object))==0,]) > 0, ncol(object) > 0,
            !is.null(colnames(object)), is.logical(plot), colnames(object)==rownames(pheno.df), metavars %in% colnames(pheno.df))
  fm <- stats::as.formula(paste0("~", paste(metavars, collapse="+")))
  des <- stats::model.matrix(fm, pheno.df)
  if (ncol(des) > ncol(object)) {
    stop("Too many covariates")}

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale.=FALSE)
  num.pc <- min(num.pc, ncol(object))
  pcs <- paste0("PC", 1:num.pc)
  dat <- cbind(pheno.df, pca$x[rownames(pheno.df), pcs, drop=FALSE])

  pvals <- sapply(pcs, FUN=function(pc) {
    fm <- stats::as.formula(paste0(pc, "~", paste(metavars, collapse="+")))
    lmfit <- stats::lm(fm, data = dat)
    av <- stats::aov(lmfit)
    avRes <- as.data.frame(summary(av)[[1]])
    avRes <- avRes[metavars, ] %>%
      dplyr::pull(`Pr(>F)`)
  })
  colnames(pval) <- pcs

  if (plot) {plot_pvals(pvals, name=name, width=width, height=height, xrot=xrot, yrot=yrot)}

  return(pvals)
}
