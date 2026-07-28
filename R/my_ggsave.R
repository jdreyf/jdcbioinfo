#' Wrapper for ggplot2::ggsave
#'
#' Call ggsave twice to save the plot as both PNG and PDF
#'
#' @param name The name of plots (without the format extension)
#' @inheritParams ggplot2::ggsave
#' @return A character vector of file names of the PNG and PDF
#' @export

my_ggsave <- function(name,
                      plot = get_last_plot(),
                      width = NA,
                      height = NA,
                      units = "in",
                      dpi = 300,
                      limitsize = FALSE,
                      bg = "white",
                      ...) {
  formats <- c("png", "pdf")
  fileanmes <- paste(name, formats, sep = ".")
  for (fileanme in fileanmes) {
    ggsave(filename = fileanme,
           plot = plot,
           width = width,
           height = height,
           units = "in",
           dpi = dpi,
           limitsize = limitsize,
           bg = bg,
           ...
    )
  }
  return(invisible(fileanmes))
}
