#' Fisher exact test for 2-way Venn
#'
#' Fisher exact test for enrichment of overlap in 2-way Venn.
#'
#' @param mat.sig Numeric matrix with two columns whose elements are 1, -1, or 0 for significantly up, down, or
#' insignificant genes (rows) in a test (column) if `directional=TRUE`.
#' If `directional=FALSE`, we only consider significant and non significant genes, so elements are 0 or 1.
#' @param directional Logical; are we considering whether significant genes are up or down?
#' Otherwise, we only consider if genes are significant or not.
#' @return Data frame with statistics from Fisher exact test.
#' @export

venn_fisher_test <- function(mat.sig, directional=TRUE){

  if (!is.matrix(mat.sig)) mat.sig <- as.matrix(mat.sig)
  stopifnot(!all(mat.sig == 0), ncol(mat.sig)==2, mat.sig %in% c(-1, 0, 1))
  if (!directional) stopifnot(mat.sig %in% c(0, 1))

  mat.sig.up <- mat.sig
  mat.sig.up[mat.sig.up != 1] <- 2
  mat.sig.up <- as.data.frame(mat.sig.up)
  mat.sig.up[,1] <- factor(mat.sig.up[,1], levels=c(1,2))
  mat.sig.up[,2] <- factor(mat.sig.up[,2], levels=c(1,2))
  tab.up <- table(mat.sig.up[,1], mat.sig.up[,2])
  res.up <- stats::fisher.test(tab.up, alternative="greater")
  vec.up <- c(res.up$p.value, res.up$estimate, res.up$conf.int)

  if (directional){
    mat.sig.down <- mat.sig
    mat.sig.down[mat.sig.down != -1] <- 0
    mat.sig.down <- as.data.frame(mat.sig.down)
    mat.sig.down[,1] <- factor(mat.sig.down[,1], levels=c(-1,0))
    mat.sig.down[,2] <- factor(mat.sig.down[,2], levels=c(-1,0))
    tab.down <- table(mat.sig.down[,1], mat.sig.down[,2])
    res.down <- stats::fisher.test(tab.down, alternative="greater")
    vec.down <- c(res.down$p.value, res.down$estimate, res.down$conf.int)

    res <- data.frame(vec.up, vec.down)
    dimnames(res) <- list(c("p-Value", "Odds Ratio", "Lower Bound 95% CI for Odds Ratio", "Upper Bound 95% CI for Odds Ratio"),
                          c("Overlapping of Up Regulation", "Overlapping of Down Regulation"))
  } else {
    res <- data.frame(vec.up)
    dimnames(res) <- list(c("p-Value", "Odds Ratio", "Lower Bound 95% CI for Odds Ratio", "Upper Bound 95% CI for Odds Ratio"),
                          "Overlapping of Enrichment")
  }
  return(res)
}
