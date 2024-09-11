#' Make a barplot for GSEA results
#'
#' Make a barplot showing the NES of top pathways of GSEA results
#'
#' @param gseaObj GSEA object
#' @param title plot title
#' @param ntop number of up and down top pathways to shown
#' @return ggplot2 object.
#' @export

make_gsea_barplot <- function(gseaObj, title = "GSEA", ntop = 10) {
  ggp <- gseaObj@result %>%
    dplyr::mutate(Direction = ifelse(NES > 0, "Positive", "Negative")) %>%
    dplyr::group_by(Direction) %>%
    dplyr::top_n(ntop, wt = abs(NES)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Description = factor(Description, levels = unique(Description))) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = Description, y = NES, fill = Direction)) + ggplot2::geom_col() + ggplot2::coord_flip()
  ggp <- ggp + ggplot2::ggtitle(title) + ggplot2::scale_fill_manual(values = c(Negative = "blue", Positive = "red"))
  ggp <- ggp + ggplot2::scale_x_discrete(labels = scales::label_wrap(60))
  graphics::plot(ggp)

  return(invisible(ggp))
}
