#' Make cnets plot for a list of GSEA results
#'
#' Make cnet plots for the community-representing pathways
#'
#' @param gseaRes a named list of GSEA object
#' @param fc a list of Z-scores with the same names as gseaRes
#' @param resolution optional resolution parameter
#' @param algorithm either "louvain" or "leiden"
#' @param label label for color bar
#' @return a list of ggplot2 objects.
#' @export

make_gsea_cnet <- function(gseaRes, fc, resolution = 1, algorithm = "louvain", label = "Z-scores") {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package \"igraph\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("enrichplot", quietly = TRUE)) stop("Package \"enrichplot\" must be installed to use this function.", call. = FALSE)
  gseaRes <- lapply(gseaRes, FUN = enrichplot::pairwise_termsim, method = "JC")

  pltList <- list()
  for(nm in names(gseaRes)) {
    gsea <- gseaRes[[nm]]
    mat <- gsea@termsim
    mat[is.na(mat)] <- 0
    net <- igraph::graph_from_adjacency_matrix(mat, mode = "max", weighted = TRUE)

    if (algorithm == "louvain") {
      cluster.out <- igraph::cluster_louvain(net, resolution = resolution)
    } else if (algorithm == "leiden") {
      cluster.out <- igraph::cluster_leiden(net, resolution_parameter = resolution, objective_function = "modularity")
    }

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
    plt <- enrichplot::cnetplot(gsea2, foldChange = fc[[nm]], categorySize = ~itemNum,  showCategory = length(pwys), node_label = "category",
                    color_category = "orange", size_category = 1.5)
    plt <- plt + ggplot2::scale_colour_gradient2(name = label, low = "blue", mid = "white", high = "red")
    plt <- plt + ggplot2::ggtitle(nm)
    plt <- plt + ggplot2::theme(
      # 1. Center the title, increase its font size, and add a bottom margin
      plot.title = element_text(
        hjust = 0.5,                   # Centers the title text (0 = left, 1 = right)
        size = 18,                     # Increases title font size
        face = "bold",                 # Makes title bold
        margin = margin(b = 20)        # Adds 20pt margin between title and plot surface
      ),
      # 2. Add an explicit margin around the legend container
      legend.margin = margin(t = 10, r = 15, b = 10, l = 15),
      # 3. Add space explicitly between the legend body and the main plot panel
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 20)
    )
    graphics::plot(plt)

    pltList[[nm]] <- plt
  }
  return(invisible(pltList))
}
