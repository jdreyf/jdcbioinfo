#' Association tests for meta variables and principal components
#'
#' Association tests for meta variables and principal components.
#'
#' @param object Matrix-like object with features (e.g. genes) as rows and samples as columns.
#' @param pheno.df Data frame with rows as samples and columns as phenotypes.
#' @param metavars A vector of column names in pheno.df.
#' @param num.pc Number of principal components
#' @param seed Random seed.
#' @param plot Logical for plotting the p-values.
#' @param name Name of the plot file.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @param xrot rot angle (in degrees) by which the x- axis labels are to be rotated
#' @param yrot rot angle (in degrees) by which the y- axis labels are to be rotated
#' @return matrix of p-values.
#' @export

pc_assoc <- function(object, pheno.df, metavars, num.pc=10, seed=100, plot=TRUE, name="assoc_with_pc_pvals", width=8, height=7, xrot=0, yrot=0) {
  stopifnot(limma::isNumeric(object), nrow(object[rowSums(is.na(object))==0,]) > 0, ncol(object) > 0,
            !is.null(colnames(object)), is.logical(plot), colnames(object)==rownames(pheno.df), metavars %in% colnames(pheno.df))

  pca <- stats::prcomp(t(object[rowSums(is.na(object))==0,]), scale.=FALSE)
  num.pc <- min(num.pc, ncol(object))
  obj1 <- pheno.df[, metavars, drop=FALSE]
  obj2 <- pca$x[rownames(pheno.df), 1:num.pc, drop=FALSE]
  pvals <-  multi_assoc(object1=obj1, object2=obj2, seed=seed, plot=plot, name=name, width=width, height=height, xrot=xrot, yrot=yrot)
  return(pvals)
}
