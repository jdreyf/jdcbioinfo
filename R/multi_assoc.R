#' Association tests for two groups of variables
#'
#' Association tests for two groups of variables
#'
#' @param object1 Vector, matrix or data.frame of group 1 variables, each column is 1 variable.
#' @param object2 Vector, matrix or data.frame of group 2 variables, each column is 1 variable.
#' @param seed Random seed.
#' @param plot Logical for plotting the p-values.
#' @param name Name of the plot file.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @return Matrix of p-values.
#' @export

multi_assoc <- function(object1, object2, seed=100, plot=TRUE, name="assoc_pvals", width=8, height=7) {
  if (is.vector(object1)) {
    nm1 <- deparse(substitute(object1))
    object1 <- as.matrix(object1)
    colnames(object1) <- nm1
  }
  if (is.vector(object2)) {
    nm2 <- deparse(substitute(object2))
    object2 <- as.matrix(object2)
    colnames(object2) <- nm2
  }
  stopifnot(nrow(object1)==nrow(object2))

  pvals <- matrix(NA, nrow=ncol(object1), ncol=ncol(object2), dimnames=list(colnames(object1), colnames(object2)))
  for (i in 1:ncol(object1)) {
    for (j in 1:ncol(object2)) {
      pvals[i, j] <- assoc(v1=object1[,i], v2=object2[,j], seed=seed)
    }
  }

  if (plot) {plot_pvals(pvals, name=name, width=width, height=height)}

  return(pvals)
}
