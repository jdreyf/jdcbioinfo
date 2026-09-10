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
  for (format in formats) {
    fileanme <- paste(name, format, sep = ".")
    device <- ifelse(format == "png", grDevices::png, grDevices::cairo_pdf)
    ggsave(filename = fileanme,
           plot = plot,
           device = device,
           width = width,
           height = height,
           units = "in",
           dpi = dpi,
           limitsize = limitsize,
           bg = bg,
           ...
    )
  }
  fileanmes <- paste(name, formats, sep = ".")
  return(invisible(fileanmes))
}
