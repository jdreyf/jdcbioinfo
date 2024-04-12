#' Make a barplot for GSEA results
#'
#' Make a barplot showing the NES of top pathways of GSEA results
#'
#' @param gseaObj GSEA object
#' @param ntop number of up and down top patwys to shown
#' @return ggplot2 object.
#' @export

make_gsea_barplot <- function(gseaObj, ntop = 10) {
  ggp <- gseaObj@result %>%
    mutate(Direction = ifelse(NES > 0, "Positive", "Negative")) %>%
    group_by(Direction) %>%
    top_n(ntop, wt = abs(NES)) %>%
    ungroup() %>%
    arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description))) %>%
    ggplot(mapping = aes(x = Description, y = NES, fill = Direction)) + geom_col() + coord_flip()
  ggp <- ggp + ggtitle(nm) + scale_fill_manual(values = c(Negative = "blue", Positive = "red"))
  ggp <- ggp + scale_x_discrete(labels = scales::label_wrap(60))
  plot(ggp)

  return(invisible(ggp))
}
