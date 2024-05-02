#' Make cnets plot for a list of GSEA results
#'
#' Make cnet plots for the community-representing pathways
#'
#' @param gseaRes a named list of GSEA object
#' @param fc a list of Z-scores with the same names as gseaRes
#' @resolution Optional resolution parameter
#' @label label for color bar
#' @return a list of ggplot2 objects.
#' @export

make_gsea_cnet <- function(gseaRes, fc, resolution = 1, label = "Z-scores") {
  gseaRes <- lapply(gseaRes, FUN = enrichplot::pairwise_termsim, method = "JC")

  pltList <- list()
  for(nm in names(gseaRes)) {
    gsea <- gseaRes[[nm]]
    mat <- gsea@termsim
    mat[is.na(mat)] <- 0
    net <- igraph::graph_from_adjacency_matrix(mat, mode = "max", weighted = TRUE)

    cluster.out <- igraph::cluster_louvain(net, resolution = resolution)

    pwys <- sapply(unique(cluster.out$membership), FUN = function(i) {
      x <- cluster.out$names[cluster.out$membership == i]

      gsea@result %>%
        dplyr::filter(Description %in% x) %>%
        dplyr::top_n(1, wt = -1*pvalue) %>%
        dplyr::pull(Description)})

    gsea2 <- gsea
    gsea2@result <- gsea2@result %>%
      dplyr::filter(Description %in% pwys)

    set.seed(100)
    plt <- enrichplot::cnetplot(gsea2, categorySize = "geneNum",  showCategory = length(pwys), node_label = "category",
                    color.params = list(foldChange = fc[[nm]], edge = TRUE, category = "orange"), cex.params = list(category_node = 1.5))
    plt <- plt + ggplot2::scale_colour_gradient2(name = label, low = "blue", mid = "white", high = "red") + ggplot2::ggtitle(nm)
    graphics::plot(plt)

    pltList[[nm]] <- plt
  }
  return(invisible(pltList))
}
