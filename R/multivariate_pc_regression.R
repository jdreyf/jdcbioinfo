#' Multivariate regression for meta variables and principal components
#'
#' Multivariate regression for meta variables and principal components
#'
#' @inheritParams pc_assoc
#' @param metavars Character vector of column names in pheno.df to include as independent variables in the regression.
#' @param formula Optional. A formula object specifying the regression model. If provided, metavars is ignored.
#' @param num.pc Number of principal components to analyze. Default is 10.
#' @param scale Logical indicating whether to scale the data before PCA. Default is FALSE.
#' @param plot Logical indicating whether to plot the p-values heatmap. Default is TRUE.
#' @param name Name prefix for the output plot file. Default is "assoc_with_pc_pvals".
#' @param width Width of the output plot in inches. Default is 8.
#' @param height Height of the output plot in inches. Default is 7.
#' @param xrot Rotation angle for x-axis labels in the plot. Default is 0.
#' @param yrot Rotation angle for y-axis labels in the plot. Default is 0.
#' @return Matrix of p-values.
#' @export

multivariate_pc_regression <- function(object, pheno.df, metavars=NULL, formula=NULL, num.pc=10, scale=FALSE, plot=TRUE, name="assoc_with_pc_pvals", width=8, height=7, xrot=0, yrot=0) {
  stopifnot(limma::isNumeric(object), nrow(object[rowSums(is.na(object))==0,]) > 0, ncol(object) > 0,
            !is.null(colnames(object)), is.logical(plot), colnames(object)==rownames(pheno.df))
  stopifnot(!is.null(metavars) | !is.null(formula))
  stopifnot(is.null(metavars) | is.null(formula))

  if (!is.null(formula)) {
    fm <- formula
    metavars <- all.vars(formula)
  } else {
    fm <- stats::as.formula(paste0("~", paste(metavars, collapse="+")))
  }
  stopifnot(metvars %in% colnames(pheno.df))

  des <- stats::model.matrix(fm, pheno.df)
  if (ncol(des) > ncol(object)) {
    stop("Too many covariates")}

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale.=scale)
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
  colnames(pvals) <- pcs
  rownames(pvals) <- metavars
  pvals[pvals < 0] <- 0
  pvals[pvals > 1] <- 1
  if (plot) {plot_pvals(pvals, name=name, width=width, height=height, xrot=xrot, yrot=yrot)}

  return(pvals)
}
