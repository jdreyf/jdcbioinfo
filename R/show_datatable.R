#' Display an interactive HTML table
#'
#' @param x A data frame or matrix.
#' @param caption Optional table caption.
#' @param digits Decimal places to display in numeric columns.
#' @param highlight_col Name of the numeric column used for highlighting.
#' @param cutoff Highlight rows whose selected value is <= this cutoff.
#' @param highlight_color Background color for highlighted rows.
#' @param options Options passed to DT::datatable().
#' @return A DT HTML widget.
#' @export
show_datatable <- function(
    x, caption = NULL, digits = 3,
    highlight_col = NULL, cutoff = NULL,
    highlight_color = "#fff3cd",
    options = list(pageLength = 10)) {

  x <- as.data.frame(x)

  if (xor(is.null(highlight_col), is.null(cutoff))) {
    stop("Supply both `highlight_col` and `cutoff`, or neither.",
         call. = FALSE)
  }

  if (!is.null(highlight_col)) {
    if (!is.character(highlight_col) ||
        length(highlight_col) != 1L ||
        is.na(highlight_col) ||
        sum(names(x) == highlight_col) != 1L) {
      stop("`highlight_col` must identify one existing column.",
           call. = FALSE)
    }

    if (!is.numeric(x[[highlight_col]])) {
      stop("The highlighting column must be numeric.", call. = FALSE)
    }

    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        !is.finite(cutoff)) {
      stop("`cutoff` must be one finite numeric value.", call. = FALSE)
    }
  }

  tbl <- DT::datatable(
    x,
    rownames = FALSE,
    escape = TRUE,
    caption = if (!is.null(caption)) htmltools::tags$caption(caption),
    options = options
  )

  numeric_cols <- which(vapply(x, is.numeric, logical(1)))
  if (length(numeric_cols)) {
    tbl <- DT::formatRound(tbl, numeric_cols, digits = digits)
  }

  if (!is.null(highlight_col)) {
    tbl <- DT::formatStyle(
      tbl,
      columns = highlight_col,
      target = "row",
      backgroundColor = DT::styleInterval(
        cutoff, c(highlight_color, "")
      )
    )
  }

  tbl
}
