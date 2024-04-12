#' Make cnets plot for a list of GSEA results
#'
#' Make cnet plots for the community-representing pathways
#'
#' @param gseaRes a named list of GSEA object
#' @param fc a list of Z-scores with the same names as gseaRes
#' @return a list of ggplot2 objects.
#' @export

make_gsea_cnet <- function(gseaRes, fc) {
  gseaRes <- lapply(gseaRes, FUN = pairwise_termsim, method = "JC")

  pltList <- list()
  for(nm in names(gseaRes)) {
    gsea <- gseaRes[[nm]]
    mat <- gsea@termsim
    mat[is.na(mat)] <- 0
    net <- graph_from_adjacency_matrix(mat, mode = "max", weighted = TRUE)

    cluster.out <- cluster_louvain(net)

    pwys <- sapply(unique(cluster.out$membership), FUN = function(i) {
      x <- cluster.out$names[cluster.out$membership == i]

      gsea@result %>%
        dplyr::filter(Description %in% x) %>%
        top_n(1, wt = -1*pvalue) %>%
        pull(Description)})

    gsea2 <- gsea
    gsea2@result <- gsea2@result %>%
      dplyr::filter(Description %in% pwys)

    set.seed(100)
    plt <- cnetplot(gsea2, categorySize = "geneNum",  showCategory = length(pwys), node_label = "category",
                    color.params = list(foldChange = fc[[nm]], edge = TRUE, category = "orange"), cex.params = list(category_node = 1.5))
    plt <- plt + scale_colour_gradient2(name = "Z-scores", low = "blue", mid = "white", high = "red") + ggtitle(nm)
    plot(plt)

    pltList[[nm]] <- plt
  }
  return(invisible(pltList))
}
