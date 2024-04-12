#' Make a pie chart
#'
#' Make a pie chart with labels of counts
#'
#' @param tab a table of counts and category
#' @param Count column name of counts in tab
#' @param Category column name of categories in tab
#' @param title title of the chart
#' @return ggplot2 object.
#' @export

make_piechart <- function(tab, Count = "Count", Category = "Category", title = NULL) {
  tab <- tab %>%
    dplyr::mutate(csum = rev(cumsum(rev(!!rlang::sym(Count)))),
           pos = !!sym(Count)/2 + lead(csum, 1),
           pos = if_else(is.na(pos), !!sym(Count)/2, pos))

  plt <- ggplot2::ggplot(tab, ggplot2::aes(x = "" , y = !!rlang::sym(Count), fill = forcats::fct_inorder(!!rlang::sym(Category)))) +
    ggplot2::geom_col(width = 1, color = 1) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_brewer(palette = "Set3") +
    ggrepel::geom_label_repel(ggplot2::aes(y = pos, label = !!rlang::sym(Count)),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    ggplot2::guides(fill = guide_legend(title = Category)) +
    ggplot2::theme_void()

  if (!is.null(title)) {
    plt <- plt + ggplot2::ggtitle(title)
  }
  graphics::plot(plt)

  return(invisible(plt))
}
