#' Wrapper for knitr::kable
#'
#' Call knitr::kable to print out the table in HTML format. This is useful for R Markdown reports.
#'
#' @param x The table to be printed
#' @param caption A caption for the table
#' @param digits The number of digits to display
#' @return NULL
#' @export

show_table <- function(x, caption, digits = 3) {
  cat('<div class="report-table">\n')
  print(knitr::kable(x, format = "html", row.names = FALSE,
                     caption = caption, digits = digits, escape = TRUE))
  cat('</div>\n\n')
}
