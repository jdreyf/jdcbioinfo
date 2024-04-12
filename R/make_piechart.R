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
    mutate(csum = rev(cumsum(rev(!!sym(Count)))),
           pos = !!sym(Count)/2 + lead(csum, 1),
           pos = if_else(is.na(pos), !!sym(Count)/2, pos))

  plt <- ggplot(tab, aes(x = "" , y = !!sym(Count), fill = fct_inorder(!!sym(Category)))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set3") +
    geom_label_repel(aes(y = pos, label = !!sym(Count)),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = Category)) +
    theme_void()

  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  plot(plt)

  return(invisible(plt))
}
